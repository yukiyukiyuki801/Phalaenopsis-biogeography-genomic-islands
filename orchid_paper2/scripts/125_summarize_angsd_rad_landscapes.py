#!/usr/bin/env python3
"""Aggregate independent ANGSD estimates over operational RAD loci."""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

GROUP_N = {"PAM": 88, "PAP": 78, "PSA": 7}
PAIRS = (("PAM", "PAP"), ("PAM", "PSA"), ("PAP", "PSA"))


def exp_log(value: str) -> float:
    if value in ("-inf", "nan", "NA"):
        return 0.0
    return math.exp(float(value))


def tajima_from_thetas(theta_pi: float, theta_w: float, chromosomes: int) -> float:
    """Plug-in Tajima D used for locus aggregation of ANGSD posterior thetas."""
    if chromosomes < 2:
        return math.nan
    a1 = sum(1.0 / i for i in range(1, chromosomes))
    a2 = sum(1.0 / (i * i) for i in range(1, chromosomes))
    b1 = (chromosomes + 1.0) / (3.0 * (chromosomes - 1.0))
    b2 = 2.0 * (chromosomes * chromosomes + chromosomes + 3.0) / (
        9.0 * chromosomes * (chromosomes - 1.0)
    )
    c1 = b1 - 1.0 / a1
    c2 = b2 - (chromosomes + 2.0) / (a1 * chromosomes) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    expected_s = theta_w * a1
    variance = e1 * expected_s + e2 * expected_s * max(0.0, expected_s - 1.0)
    return (theta_pi - theta_w) / math.sqrt(variance) if variance > 0 else math.nan


def quantile(values, probability):
    values = sorted(x for x in values if math.isfinite(x))
    if not values:
        return math.nan
    index = max(0, min(len(values) - 1, math.ceil(probability * len(values)) - 1))
    return values[index]


def fmt(value):
    return "NA" if isinstance(value, float) and not math.isfinite(value) else f"{value:.10g}" if isinstance(value, float) else value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--min-bases", type=int, default=100)
    args = parser.parse_args()
    work = args.work
    results = work / "results"
    results.mkdir(exist_ok=True)

    site_locus = {}
    with (work / "inputs" / "site_to_rad_locus.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            site_locus[(row["chrom"], int(row["pos"]))] = row["locus_id"]

    locus_info = {}
    with (work / "inputs" / "rad_locus_target_summary.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locus_info[row["locus_id"]] = row

    snp_counts = {}
    with (work / "inputs" / "pair_locus_informative_snps.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            snp_counts[(row["pair"], row["locus_id"])] = int(row["informative_snps"])

    # Posterior per-site theta estimates are log-scaled in thetaStat output.
    pop_stats = {}
    for group in GROUP_N:
        accum = defaultdict(lambda: [0, 0.0, 0.0])  # sites, thetaW, thetaPi
        with gzip.open(work / "theta" / f"{group}.thetas.tsv.gz", "rt") as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                fields = raw.split()
                key = (fields[0], int(fields[1]))
                locus = site_locus.get(key)
                if locus is None:
                    continue
                record = accum[locus]
                record[0] += 1
                record[1] += exp_log(fields[2])
                record[2] += exp_log(fields[3])
        for locus, (sites, theta_w, theta_pi) in accum.items():
            pop_stats[(group, locus)] = {
                "bases": sites,
                "pi": theta_pi / sites if sites else math.nan,
                "theta_w": theta_w / sites if sites else math.nan,
                "tajima_d": tajima_from_thetas(theta_pi, theta_w, 2 * GROUP_N[group]),
            }

    # Fixed common REF/ALT coding makes frequencies directly comparable.
    frequencies = {}
    for group in GROUP_N:
        with gzip.open(work / "saf" / f"{group}.mafs.gz", "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                frequencies[(group, row["chromo"], int(row["position"]))] = float(row["knownEM"])

    # realSFS fst print gives chromosome, position, numerator (a), denominator (b).
    fst_accum = {}
    for group1, group2 in PAIRS:
        pair = f"{group1}-{group2}"
        accum = defaultdict(lambda: [0, 0.0, 0.0])
        with gzip.open(work / "fst" / f"{pair}.per_site.tsv.gz", "rt") as handle:
            for raw in handle:
                fields = raw.split()
                if len(fields) < 4 or fields[0].lower().startswith(("chr", "region")):
                    continue
                try:
                    pos, numerator, denominator = int(fields[1]), float(fields[2]), float(fields[3])
                except ValueError:
                    continue
                locus = site_locus.get((fields[0], pos))
                if locus is None or not (math.isfinite(numerator) and math.isfinite(denominator)):
                    continue
                record = accum[locus]
                record[0] += 1
                record[1] += numerator
                record[2] += denominator
        fst_accum[pair] = accum

    all_rows = []
    for group1, group2 in PAIRS:
        pair = f"{group1}-{group2}"
        dxy_accum = defaultdict(lambda: [0, 0.0])
        for (chrom, pos), locus in site_locus.items():
            p1 = frequencies.get((group1, chrom, pos))
            p2 = frequencies.get((group2, chrom, pos))
            if p1 is None or p2 is None:
                continue
            dxy_accum[locus][0] += 1
            dxy_accum[locus][1] += p1 * (1.0 - p2) + p2 * (1.0 - p1)
        for locus, (dxy_sites, dxy_sum) in dxy_accum.items():
            info = locus_info[locus]
            first = pop_stats.get((group1, locus), {})
            second = pop_stats.get((group2, locus), {})
            fst_sites, fst_num, fst_den = fst_accum[pair].get(locus, (0, 0.0, 0.0))
            all_rows.append({
                "pair": pair, "group1": group1, "group2": group2,
                "locus_id": locus, "chrom": info["chrom"],
                "start": int(info["first_target_pos"]), "end": int(info["last_target_pos"]),
                "target_sites": int(info["target_sites"]),
                "informative_snps": snp_counts.get((pair, locus), 0),
                "dxy_callable_bases": dxy_sites,
                "fst_callable_bases": fst_sites,
                "pi_bases_group1": first.get("bases", 0),
                "pi_bases_group2": second.get("bases", 0),
                "dxy": dxy_sum / dxy_sites if dxy_sites else math.nan,
                "fst_angsd": fst_num / fst_den if fst_den > 0 else math.nan,
                "pi_group1": first.get("pi", math.nan),
                "pi_group2": second.get("pi", math.nan),
                "pi_pair_mean": (first.get("pi", math.nan) + second.get("pi", math.nan)) / 2.0,
                "tajima_d_group1": first.get("tajima_d", math.nan),
                "tajima_d_group2": second.get("tajima_d", math.nan),
            })

    fields = list(all_rows[0]) + ["shang_model_5pct"]
    for min_snps in (0, 1, 2, 3):
        retained = [row for row in all_rows if
                    row["informative_snps"] >= min_snps and
                    row["dxy_callable_bases"] >= args.min_bases and
                    row["fst_callable_bases"] >= args.min_bases and
                    row["pi_bases_group1"] >= args.min_bases and
                    row["pi_bases_group2"] >= args.min_bases and
                    math.isfinite(row["fst_angsd"]) and math.isfinite(row["dxy"])]
        cutoffs = {}
        for pair in ("PAM-PAP", "PAM-PSA", "PAP-PSA"):
            pair_rows = [row for row in retained if row["pair"] == pair]
            cutoffs[pair] = {
                "fst05": quantile([r["fst_angsd"] for r in pair_rows], 0.05),
                "fst95": quantile([r["fst_angsd"] for r in pair_rows], 0.95),
                "dxy05": quantile([r["dxy"] for r in pair_rows], 0.05),
                "dxy95": quantile([r["dxy"] for r in pair_rows], 0.95),
            }
        path = results / f"angsd_rad_locus_statistics_min{min_snps}snp.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in sorted(retained, key=lambda x: (x["pair"], x["chrom"], x["start"])):
                q = cutoffs[row["pair"]]
                fst, dxy = row["fst_angsd"], row["dxy"]
                model = "unclassified"
                if fst >= q["fst95"] and dxy >= q["dxy95"]:
                    model = "divergence_with_gene_flow"
                elif fst >= q["fst95"] and dxy <= q["dxy05"]:
                    model = "recurrent_selection"
                elif fst >= q["fst95"] and q["dxy05"] < dxy < q["dxy95"]:
                    model = "allopatric_selection"
                elif fst <= q["fst05"] and dxy >= q["dxy95"]:
                    model = "balancing_selection"
                output = {key: fmt(value) for key, value in row.items()}
                output["shang_model_5pct"] = model
                writer.writerow(output)

    with (results / "estimator_notes.tsv").open("w") as handle:
        handle.write("quantity\testimator\n")
        handle.write("pi\tFolded 1D-SFS posterior pairwise theta divided by retained sites\n")
        handle.write("tajima_d\tPlug-in Tajima D from summed ANGSD posterior theta_pi and theta_W\n")
        handle.write("fst\trealSFS two-population folded 2D-SFS ratio of summed numerator and denominator\n")
        handle.write("dxy\tANGSD genotype-likelihood allele frequencies with common fixed REF/ALT coding, including invariant sites\n")


if __name__ == "__main__":
    main()

