#!/usr/bin/env python3
"""Recalculate species-only pi/Tajima's D and plot corrected RAD outlier Venns.

FST and dXY sets are defined for the three pairwise comparisons.  Pi and
Tajima's D sets are defined independently for PAM, PAP, and PSA, using each
species' callable sites rather than the pair-specific callable intersections.
"""
from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/results/rad_locus_statistics"
VCF = ROOT / "orchid_paper2/05_all_sites_vcf/filtered/orchid_173.paper1_like.biallelic_snps_invariant.vcf.gz"
BED = ROOT / "orchid_paper2/05_all_sites_vcf/callable/cohort_callable.bed"
PAIRS = ["PAM-PAP", "PAM-PSA", "PAP-PSA"]
SPECIES = ["PAM", "PAP", "PSA"]
GROUP_N = {"PAM": 88, "PAP": 78, "PSA": 7}
COLORS = ["#4C78A8", "#F58518", "#54A24B"]
CALL_RATE = 0.70


def load_bed(path, gap=500):
    by_chrom, bounds = defaultdict(list), {}
    cluster = 0
    last_chrom, last_end = None, None
    with open(path) as handle:
        for line in handle:
            chrom, start, end = line.rstrip().split("\t")[:3]
            start, end = int(start), int(end)
            if chrom != last_chrom or last_end is None or start - last_end > gap:
                cluster += 1
            by_chrom[chrom].append((start + 1, end, cluster))
            if cluster not in bounds:
                bounds[cluster] = [chrom, start + 1, end]
            else:
                bounds[cluster][1] = min(bounds[cluster][1], start + 1)
                bounds[cluster][2] = max(bounds[cluster][2], end)
            last_chrom, last_end = chrom, end
    return by_chrom, bounds


def parse_info(value):
    return dict(x.split("=", 1) for x in value.split(";") if "=" in x)


def scalar(info, key):
    value = info.get(key, ".").split(",")[0]
    return float(value) if value not in ("", ".") else 0.0


def pi_site(ac, an):
    p = ac / an
    return 2.0 * p * (1.0 - p) * an / (an - 1.0)


def tajima_d(theta_pi, segregating, n):
    n = int(n)
    if n < 2 or segregating < 2:
        return math.nan
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n*n + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2.0) / (a1*n) + a2 / (a1*a1)
    e1 = c1 / a1
    e2 = c2 / (a1*a1 + a2)
    denominator = math.sqrt(e1*segregating + e2*segregating*(segregating-1.0))
    return (theta_pi - segregating/a1) / denominator if denominator > 0 else math.nan


def calculate_species_statistics(common_loci):
    intervals, bounds = load_bed(BED)
    cursor = defaultdict(int)
    dat = defaultdict(lambda: {"callable_bases": 0, "pi_sum": 0.0,
                               "segregating_sites": 0, "min_an": 10**9})
    with gzip.open(VCF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t", 9)
            chrom, pos, alt = fields[0], int(fields[1]), fields[4]
            iv = intervals[chrom]
            j = cursor[chrom]
            while j < len(iv) and pos > iv[j][1]:
                j += 1
            cursor[chrom] = j
            if j == len(iv) or pos < iv[j][0]:
                raise RuntimeError(f"Position outside callable BED: {chrom}:{pos}")
            cluster = iv[j][2]
            locus = f"RAD{cluster:05d}"
            if locus not in common_loci:
                continue
            info = parse_info(fields[7])
            for species in SPECIES:
                an = scalar(info, f"AN_{species}")
                ac = scalar(info, f"AC_{species}")
                if an < 2 * GROUP_N[species] * CALL_RATE:
                    continue
                x = dat[(species, cluster)]
                x["callable_bases"] += 1
                x["pi_sum"] += pi_site(ac, an)
                x["min_an"] = min(x["min_an"], int(an))
                if alt != "." and 0 < ac < an:
                    x["segregating_sites"] += 1

    rows = []
    for locus in sorted(common_loci):
        cluster = int(locus[3:])
        chrom, start, end = bounds[cluster]
        for species in SPECIES:
            x = dat[(species, cluster)]
            bases = x["callable_bases"]
            pi = x["pi_sum"] / bases if bases else math.nan
            td = tajima_d(x["pi_sum"], x["segregating_sites"], x["min_an"])
            rows.append({"species": species, "locus_id": locus, "chrom": chrom,
                         "start": start, "end": end, "callable_bases": bases,
                         "segregating_sites": x["segregating_sites"],
                         "min_an": x["min_an"] if bases else np.nan,
                         "pi": pi, "tajima_d": td})
    return pd.DataFrame(rows)


def venn_regions(sets, labels):
    a, b, c = (sets[x] for x in labels)
    return {
        f"{labels[0]} only": a-b-c,
        f"{labels[1]} only": b-a-c,
        f"{labels[2]} only": c-a-b,
        f"{labels[0]} & {labels[1]} only": (a & b)-c,
        f"{labels[0]} & {labels[2]} only": (a & c)-b,
        f"{labels[1]} & {labels[2]} only": (b & c)-a,
        "All three": a & b & c,
    }


def main():
    pair = pd.read_csv(BASE / "rad_locus_statistics.tsv", sep="\t", na_values="NA")
    pair = pair[pair.chrom != "NC_017609.1"].copy()
    common = set.intersection(*(set(pair.loc[pair.pair == p, "locus_id"]) for p in PAIRS))
    pair = pair[pair.locus_id.isin(common)].copy()

    species = calculate_species_statistics(common)
    species.to_csv(BASE / "rad_species_locus_statistics.tsv", sep="\t", index=False,
                   na_rep="NA", float_format="%.10g")

    definitions = []
    thresholds = []
    for title, metric, direction in (("High FST", "hudson_fst", "upper"),
                                     ("High dXY", "dxy", "upper")):
        sets = {}
        for label in PAIRS:
            z = pair[pair.pair == label]
            threshold = z[metric].quantile(.95)
            sets[label] = set(z.loc[z[metric] >= threshold, "locus_id"])
            thresholds.append([title, "pair", label, metric, direction, .95, threshold,
                               len(sets[label])])
        definitions.append((f"{title} (upper 5%)", PAIRS, sets))
    for title, metric, direction, quantile in (
            ("Low pi", "pi", "lower", .05),
            ("Low Tajima's D", "tajima_d", "lower", .05),
            ("High Tajima's D", "tajima_d", "upper", .95)):
        sets = {}
        for label in SPECIES:
            z = species[species.species == label]
            threshold = z[metric].dropna().quantile(quantile)
            hit = z[metric] <= threshold if direction == "lower" else z[metric] >= threshold
            sets[label] = set(z.loc[hit, "locus_id"])
            thresholds.append([title, "species", label, metric, direction, quantile, threshold,
                               len(sets[label])])
        definitions.append((f"{title} ({direction} 5%)", SPECIES, sets))

    centers = [(.40, .59), (.60, .59), (.50, .41)]
    label_xy = [(.22, .88), (.78, .88), (.50, .08)]
    number_xy = [(.25,.63),(.75,.63),(.50,.25),(.50,.67),(.36,.44),(.64,.44),(.50,.50)]
    fig, axes = plt.subplots(2, 3, figsize=(17, 11))
    summary, membership = [], []
    for ax, (title, labels, sets) in zip(axes.flat, definitions):
        for center, color in zip(centers, COLORS):
            ax.add_patch(Circle(center, .265, facecolor=color, edgecolor=color, alpha=.28, lw=2))
        for xy, label, color in zip(label_xy, labels, COLORS):
            ax.text(*xy, f"{label}\n(n={len(sets[label])})", ha="center", va="center",
                    color=color, fontweight="bold", fontsize=10)
        regions = venn_regions(sets, labels)
        for xy, (_, loci) in zip(number_xy, regions.items()):
            ax.text(*xy, str(len(loci)), ha="center", va="center", fontsize=12,
                    fontweight="bold")
        for region, loci in regions.items():
            summary.append([title, " / ".join(labels), region, len(loci)])
            membership.extend([title, " / ".join(labels), region, x] for x in sorted(loci))
        shared = set().union(*(sets[a] & sets[b] for i, a in enumerate(labels)
                               for b in labels[i+1:]))
        unit = "pairwise comparisons" if labels == PAIRS else "species"
        ax.set_title(f"{title}\nshared by >=2 {unit}: {len(shared)}; all 3: {len(regions['All three'])}",
                     fontsize=11)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off")
    axes.flat[-1].axis("off")
    fig.suptitle("Overlap of 5% RAD-locus outliers", fontsize=16, y=.985)
    fig.text(.5, .012,
             f"FST/dXY: pairwise sets; pi/Tajima's D: species sets. Common nuclear loci: n={len(common)}. Areas are schematic.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=[0, .035, 1, .96], w_pad=2.5)
    stem = BASE / "plots/rad_locus_outlier_venn_pair_and_species_top5"
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    pd.DataFrame(thresholds, columns=["panel", "set_unit", "set_label", "metric",
                 "tail", "quantile", "threshold", "n_outlier_loci"]).to_csv(
                     BASE / "pair_species_outlier_venn_top5_thresholds.tsv", sep="\t", index=False)
    pd.DataFrame(summary, columns=["panel", "set_labels", "venn_region", "n_loci"]).to_csv(
        BASE / "pair_species_outlier_venn_top5_summary.tsv", sep="\t", index=False)
    pd.DataFrame(membership, columns=["panel", "set_labels", "venn_region", "locus_id"]).to_csv(
        BASE / "pair_species_outlier_venn_top5_membership.tsv", sep="\t", index=False)
    (BASE / "pair_species_outlier_venn_top5_common_loci.txt").write_text(
        "\n".join(sorted(common)) + "\n")


if __name__ == "__main__":
    main()
