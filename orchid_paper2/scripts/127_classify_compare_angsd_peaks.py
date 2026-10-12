#!/usr/bin/env python3
"""Compare old VCF and new ANGSD RAD-locus peak classifications.

Two deliberately separate classifiers are implemented:
  1. Shang-paper FST/dXY classifier: 5% tails, no pi requirement.
  2. Orchid extended classifier: FST/dXY/pi quartiles, species-specific pi,
     with the final central-25% dXY definition for allopatric selection.
"""
from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "21_angsd_rad_landscapes" / "peak_classification"
OUT.mkdir(parents=True, exist_ok=True)

INPUTS = {
    "Old VCF": ROOT / "results" / "rad_locus_statistics" / "rad_locus_statistics.tsv",
    "New ANGSD": ROOT / "21_angsd_rad_landscapes" / "results" / "angsd_rad_locus_statistics_min3snp.tsv",
}
PAIRS = ("PAM-PAP", "PAM-PSA", "PAP-PSA")
MODELS = (
    "divergence_with_gene_flow",
    "allopatric_selection",
    "recurrent_selection",
    "balancing_selection",
    "unclassified",
)


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def quantile_nearest(values, probability):
    values = sorted(v for v in values if math.isfinite(v))
    return values[max(0, math.ceil(probability * len(values)) - 1)]


def paper_class(fst, dxy, q):
    if fst >= q["fst95"] and dxy >= q["dxy95"]:
        return "divergence_with_gene_flow"
    if fst >= q["fst95"] and dxy <= q["dxy05"]:
        return "recurrent_selection"
    if fst >= q["fst95"] and q["dxy05"] < dxy < q["dxy95"]:
        return "allopatric_selection"
    if fst <= q["fst05"] and dxy >= q["dxy95"]:
        return "balancing_selection"
    return "unclassified"


def extended_class(fst, dxy, pi, q, pi25, pi75):
    if fst >= q["fst75"] and dxy >= q["dxy75"] and pi <= pi25:
        return "divergence_with_gene_flow"
    if fst >= q["fst75"] and q["dxy375"] < dxy < q["dxy625"] and pi <= pi25:
        return "allopatric_selection"
    if fst >= q["fst75"] and dxy <= q["dxy25"] and pi <= pi25:
        return "recurrent_selection"
    if fst <= q["fst25"] and dxy >= q["dxy75"] and pi >= pi75:
        return "balancing_selection"
    return "unclassified"


all_rows = []
threshold_rows = []
for dataset, path in INPUTS.items():
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fst_field = "hudson_fst" if dataset == "Old VCF" else "fst_angsd"
    by_pair = defaultdict(list)
    for row in rows:
        by_pair[row["pair"]].append(row)

    for pair, pair_rows in by_pair.items():
        fst_values = [fnum(r[fst_field]) for r in pair_rows]
        dxy_values = [fnum(r["dxy"]) for r in pair_rows]
        q = {
            "fst05": quantile_nearest(fst_values, 0.05),
            "fst25": quantile_nearest(fst_values, 0.25),
            "fst75": quantile_nearest(fst_values, 0.75),
            "fst95": quantile_nearest(fst_values, 0.95),
            "dxy05": quantile_nearest(dxy_values, 0.05),
            "dxy25": quantile_nearest(dxy_values, 0.25),
            "dxy375": quantile_nearest(dxy_values, 0.375),
            "dxy625": quantile_nearest(dxy_values, 0.625),
            "dxy75": quantile_nearest(dxy_values, 0.75),
            "dxy95": quantile_nearest(dxy_values, 0.95),
        }
        group1, group2 = pair.split("-")
        pi_thresholds = {}
        for group, field in ((group1, "pi_group1"), (group2, "pi_group2")):
            values = [fnum(r[field]) for r in pair_rows]
            pi_thresholds[group] = (
                quantile_nearest(values, 0.25), quantile_nearest(values, 0.75)
            )
        threshold_rows.append({"dataset": dataset, "pair": pair, **q,
                               f"pi25_{group1}": pi_thresholds[group1][0],
                               f"pi75_{group1}": pi_thresholds[group1][1],
                               f"pi25_{group2}": pi_thresholds[group2][0],
                               f"pi75_{group2}": pi_thresholds[group2][1]})

        for row in pair_rows:
            fst, dxy = fnum(row[fst_field]), fnum(row["dxy"])
            species_models = []
            for group, field in ((group1, "pi_group1"), (group2, "pi_group2")):
                p25, p75 = pi_thresholds[group]
                species_models.append((group, extended_class(
                    fst, dxy, fnum(row[field]), q, p25, p75
                )))
            supported = [(g, model) for g, model in species_models if model != "unclassified"]
            # FST/dXY determine the eligible model, so supported species cannot conflict.
            extended = supported[0][1] if supported else "unclassified"
            all_rows.append({
                "dataset": dataset,
                "pair": pair,
                "locus_id": row["locus_id"],
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "informative_snps": row["informative_snps"],
                "fst": fst,
                "dxy": dxy,
                "pi_group1": fnum(row["pi_group1"]),
                "pi_group2": fnum(row["pi_group2"]),
                "paper_fst_dxy_5pct": paper_class(fst, dxy, q),
                "extended_fst_dxy_pi_central25": extended,
                "extended_model_group1": species_models[0][1],
                "extended_model_group2": species_models[1][1],
                "supporting_species": ",".join(g for g, _ in supported) or "none",
                "n_supporting_species": len(supported),
            })

with (OUT / "old_vcf_new_angsd_peak_assignments.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(all_rows)

threshold_fields = sorted({key for row in threshold_rows for key in row},
                          key=lambda x: (x not in ("dataset", "pair"), x))
with (OUT / "old_vcf_new_angsd_peak_thresholds.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=threshold_fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(threshold_rows)

classifiers = {
    "Paper: FST+dXY 5%": "paper_fst_dxy_5pct",
    "Orchid: FST+dXY+pi": "extended_fst_dxy_pi_central25",
}
summary = []
for dataset in INPUTS:
    for pair in PAIRS:
        selected = [r for r in all_rows if r["dataset"] == dataset and r["pair"] == pair]
        for classifier, field in classifiers.items():
            counts = Counter(r[field] for r in selected)
            n_classified = len(selected) - counts["unclassified"]
            for model in MODELS:
                summary.append({
                    "dataset": dataset, "pair": pair, "classifier": classifier,
                    "model": model, "n": counts[model], "n_all_loci": len(selected),
                    "n_classified": n_classified,
                    "fraction_all_loci": counts[model] / len(selected),
                    "fraction_classified_loci": (counts[model] / n_classified
                                                   if model != "unclassified" and n_classified else math.nan),
                })

with (OUT / "old_vcf_new_angsd_peak_ratio_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

# Percentage-point changes and shared-locus transitions are diagnostic only;
# the primary ratios use every retained locus in each estimator.
delta_rows = []
for pair in PAIRS:
    for classifier in classifiers:
        for model in MODELS:
            old = next(r for r in summary if r["dataset"] == "Old VCF" and
                       r["pair"] == pair and r["classifier"] == classifier and r["model"] == model)
            new = next(r for r in summary if r["dataset"] == "New ANGSD" and
                       r["pair"] == pair and r["classifier"] == classifier and r["model"] == model)
            delta_rows.append({
                "pair": pair, "classifier": classifier, "model": model,
                "old_fraction_all": old["fraction_all_loci"],
                "new_fraction_all": new["fraction_all_loci"],
                "change_all_percentage_points": 100 * (new["fraction_all_loci"] - old["fraction_all_loci"]),
                "old_fraction_classified": old["fraction_classified_loci"],
                "new_fraction_classified": new["fraction_classified_loci"],
                "change_classified_percentage_points": (
                    100 * (new["fraction_classified_loci"] - old["fraction_classified_loci"])
                    if math.isfinite(old["fraction_classified_loci"]) and math.isfinite(new["fraction_classified_loci"])
                    else math.nan),
            })
with (OUT / "old_to_new_peak_ratio_changes.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(delta_rows[0]), delimiter="\t")
    writer.writeheader(); writer.writerows(delta_rows)

transition_rows = []
for pair in PAIRS:
    old_map = {(r["pair"], r["locus_id"]): r for r in all_rows if r["dataset"] == "Old VCF"}
    new_map = {(r["pair"], r["locus_id"]): r for r in all_rows if r["dataset"] == "New ANGSD"}
    shared = sorted(set(old_map) & set(new_map))
    for classifier, field in classifiers.items():
        counts = Counter((old_map[key][field], new_map[key][field]) for key in shared if key[0] == pair)
        for (old_model, new_model), n in counts.items():
            transition_rows.append({"pair": pair, "classifier": classifier,
                                    "old_model": old_model, "new_model": new_model, "n": n})
with (OUT / "shared_locus_peak_transitions.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(transition_rows[0]), delimiter="\t")
    writer.writeheader(); writer.writerows(transition_rows)

