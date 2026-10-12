#!/usr/bin/env python3
"""Test overlap of pairwise FST/dXY and species pi/Tajima-D outliers."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/results/rad_locus_statistics"
FIG = BASE / "plots"
PAIRS = ["PAM-PAP", "PAM-PSA", "PAP-PSA"]
SPECIES = ["PAM", "PAP", "PSA"]
N_PERM = 49999
SEED = 20261002
BATCH = 500


def permute_vector(vector, strata, rng, batch):
    out = np.empty((batch, len(vector)), dtype=bool)
    for ids in strata:
        order = np.argsort(rng.random((batch, len(ids))), axis=1)
        out[:, ids] = vector[ids][order]
    return out


def overlap_counts(x):
    return np.column_stack(((x[:, 0] & x[:, 1]).sum(1),
                            (x[:, 0] & x[:, 2]).sum(1),
                            (x[:, 1] & x[:, 2]).sum(1),
                            (x[:, 0] & x[:, 1] & x[:, 2]).sum(1)))


def add_bh(frame, pcol, qcol):
    p = frame[pcol].to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(ranked) /
                               np.arange(1, len(ranked)+1))[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1)
    frame[qcol] = restored
    return frame


pair = pd.read_csv(BASE / "rad_locus_statistics.tsv", sep="\t", na_values="NA")
pair = pair[pair.chrom != "NC_017609.1"].copy()
common = set.intersection(*(set(pair.loc[pair.pair == p, "locus_id"]) for p in PAIRS))
loci = np.array(sorted(common))
pair = pair[pair.locus_id.isin(common)].copy()
species = pd.read_csv(BASE / "rad_species_locus_statistics.tsv", sep="\t", na_values="NA")

# Two ascertainment strata appropriate to the set unit. Pair statistics use
# mean pairwise callable length/informative SNP count; single-species statistics
# use mean species callable length/segregating-site count. Each is divided into
# 4 x 4 rank-based bins, and permutations retain exact set sizes within bins.
pair_cov = pair.groupby("locus_id")[["callable_bases", "informative_snps"]].mean().loc[loci]
species_cov = species.groupby("locus_id")[["callable_bases", "segregating_sites"]].mean().loc[loci]


def make_strata(frame, columns):
    a = pd.qcut(frame[columns[0]].rank(method="first"), 4, labels=False).astype(int)
    b = pd.qcut(frame[columns[1]].rank(method="first"), 4, labels=False).astype(int)
    code = (a * 4 + b).to_numpy()
    return [np.flatnonzero(code == value) for value in np.unique(code)]


pair_strata = make_strata(pair_cov, ["callable_bases", "informative_snps"])
species_strata = make_strata(species_cov, ["callable_bases", "segregating_sites"])

analyses = []
for title, metric, labels, tail, strata in (
        ("High FST", "hudson_fst", PAIRS, "upper", pair_strata),
        ("High dXY", "dxy", PAIRS, "upper", pair_strata),
        ("Low pi", "pi", SPECIES, "lower", species_strata),
        ("Low Tajima's D", "tajima_d", SPECIES, "lower", species_strata),
        ("High Tajima's D", "tajima_d", SPECIES, "upper", species_strata)):
    vectors = []
    for label in labels:
        source = pair[pair.pair == label].set_index("locus_id").loc[loci] if labels == PAIRS \
            else species[species.species == label].set_index("locus_id").loc[loci]
        threshold = source[metric].dropna().quantile(.95 if tail == "upper" else .05)
        vector = (source[metric] >= threshold if tail == "upper" else
                  source[metric] <= threshold).fillna(False).to_numpy()
        vectors.append(vector)
    analyses.append((title, labels, vectors, strata))

rows = []
for analysis_index, (metric, labels, vectors, strata) in enumerate(analyses):
    observed = overlap_counts(np.stack(vectors)[None, :, :])[0]
    null = np.empty((N_PERM, 4), dtype=np.int16)
    rng = np.random.default_rng(SEED + 1000 * analysis_index)
    start = 0
    while start < N_PERM:
        batch = min(BATCH, N_PERM - start)
        simulation = np.stack([permute_vector(v, strata, rng, batch)
                               for v in vectors], axis=1)
        null[start:start+batch] = overlap_counts(simulation)
        start += batch
    tests = [f"{labels[0]} & {labels[1]}", f"{labels[0]} & {labels[2]}",
             f"{labels[1]} & {labels[2]}", "All three"]
    for j, test in enumerate(tests):
        mean = float(null[:, j].mean())
        lower, upper = np.quantile(null[:, j], [.025, .975])
        p = (1 + np.sum(null[:, j] >= observed[j])) / (N_PERM + 1)
        rows.append({"metric": metric, "set_unit": "pair" if labels == PAIRS else "species",
                     "overlap": test, "set_1": labels[0], "set_2": labels[1],
                     "set_3": labels[2], "set_size_1": int(vectors[0].sum()),
                     "set_size_2": int(vectors[1].sum()), "set_size_3": int(vectors[2].sum()),
                     "locus_universe": len(loci), "observed_loci": int(observed[j]),
                     "null_mean_loci": mean, "null_sd_loci": float(null[:, j].std(ddof=1)),
                     "null_q025_loci": float(lower), "null_q975_loci": float(upper),
                     "above_95pct_null_interval": bool(observed[j] > upper),
                     "enrichment_ratio": observed[j] / mean if mean else np.inf,
                     "p_enriched_one_sided": p, "n_permutations": N_PERM})

result = add_bh(pd.DataFrame(rows), "p_enriched_one_sided", "q_BH_all_20_tests")
result.to_csv(BASE / "pair_species_outlier_overlap_permutation_summary.tsv",
              sep="\t", index=False)

# Four compact forest panels: null mean and 95% permutation interval versus
# observed count, with the one-sided empirical enrichment P value.
fig, axes = plt.subplots(2, 3, figsize=(18, 9))
for ax, (metric, _, _, _) in zip(axes.flat, analyses):
    z = result[result.metric == metric].reset_index(drop=True)
    y = np.arange(len(z))[::-1]
    low = z.null_mean_loci - z.null_q025_loci
    high = z.null_q975_loci - z.null_mean_loci
    ax.errorbar(z.null_mean_loci, y, xerr=np.vstack([low, high]), fmt="o",
                color="#555555", ecolor="#777777", capsize=5, lw=2,
                label="Null mean and 95% interval")
    colors = ["#D62728" if sig else "#4C78A8" for sig in z.above_95pct_null_interval]
    ax.scatter(z.observed_loci, y, marker="D", s=55, c=colors, zorder=4,
               label="Observed")
    xmax = max(z.observed_loci.max(), z.null_q975_loci.max())
    for yi, row in zip(y, z.itertuples()):
        ptext = f"P <= {1/(N_PERM+1):.5f}" if row.p_enriched_one_sided <= 1/(N_PERM+1) \
            else f"P={row.p_enriched_one_sided:.4g}"
        ax.text(xmax * 1.03 + .25, yi, ptext, va="center", fontsize=9)
    ax.set_yticks(y, z.overlap)
    ax.set_xlim(left=0, right=xmax * 1.38 + 1)
    ax.set_xlabel("Number of overlapping loci")
    ax.set_title(metric, fontweight="bold")
    ax.grid(axis="x", alpha=.2)
axes.flat[-1].axis("off")
legend_handles = [
    Line2D([0], [0], marker="D", color="none", markerfacecolor="#4C78A8",
           markeredgecolor="#4C78A8", markersize=7, label="Observed"),
    Line2D([0], [0], marker="o", color="#777777", markerfacecolor="#555555",
           markersize=6, lw=2, label="Null mean and 95% interval"),
]
fig.legend(handles=legend_handles, frameon=False, fontsize=9, loc="upper center",
           bbox_to_anchor=(.5, .925), ncol=2)
fig.text(.5, .975, "Observed RAD-outlier overlap versus permutation expectation",
         ha="center", va="top", fontsize=15, fontweight="bold")
fig.text(.5, .015,
         "49,999 ascertainment-stratified permutations; red observed diamond = above the 95% null interval; one-sided empirical P values",
         ha="center", fontsize=9.5)
fig.tight_layout(rect=[0, .04, 1, .86], w_pad=5)
stem = FIG / "rad_pair_species_outlier_overlap_permutation_95interval"
fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)

# Slide-friendly bar version. Each test has its own discrete 97.5th-percentile
# null threshold, so the red dashed threshold is drawn separately on each bar
# rather than as one misleading common horizontal line.
fig, axes = plt.subplots(2, 3, figsize=(18, 9))
for ax, (metric, _, _, _) in zip(axes.flat, analyses):
    z = result[result.metric == metric].reset_index(drop=True)
    x = np.arange(len(z))
    colors = ["#D95F02" if sig else "#4C78A8" for sig in z.above_95pct_null_interval]
    bars = ax.bar(x, z.observed_loci, color=colors, alpha=.86, width=.68,
                  edgecolor="white", linewidth=.8)
    ax.scatter(x, z.null_mean_loci, color="#444444", s=32, zorder=4)
    for xi, row, bar in zip(x, z.itertuples(), bars):
        ax.hlines(row.null_q975_loci, xi-.34, xi+.34, color="#D62728",
                  linestyle="--", linewidth=2.2, zorder=5)
        ptext = f"P <= {1/(N_PERM+1):.5f}" if row.p_enriched_one_sided <= 1/(N_PERM+1) \
            else f"P={row.p_enriched_one_sided:.4g}"
        offset = max(0.22, z.observed_loci.max() * .025)
        ax.text(xi, bar.get_height()+offset, f"{int(row.observed_loci)}\n{ptext}",
                ha="center", va="bottom", fontsize=8.5)
    short_labels = [s.replace(" & ", "\n&\n").replace("All three", "All three")
                    for s in z.overlap]
    ax.set_xticks(x, short_labels, fontsize=8.5)
    ymax = max(z.observed_loci.max(), z.null_q975_loci.max())
    ax.set_ylim(0, ymax*1.28+1)
    ax.set_ylabel("Observed overlapping loci")
    ax.set_title(metric, fontweight="bold")
    ax.grid(axis="y", alpha=.2)
axes.flat[-1].axis("off")
bar_legend = [
    Line2D([0], [0], color="#D62728", linestyle="--", lw=2.2,
           label="Upper 95% permutation threshold"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#444444",
           markeredgecolor="#444444", markersize=6, label="Null mean"),
]
fig.legend(handles=bar_legend, frameon=False, fontsize=9, loc="upper center",
           bbox_to_anchor=(.5, .925), ncol=2)
fig.text(.5, .975, "Observed overlap and the 95% permutation threshold",
         ha="center", va="top", fontsize=15, fontweight="bold")
fig.text(.5, .015,
         "Numbers above bars are observed overlapping loci; orange bars exceed their test-specific 95% null threshold",
         ha="center", fontsize=9.5)
fig.tight_layout(rect=[0, .04, 1, .86], w_pad=4)
stem = FIG / "rad_pair_species_outlier_overlap_permutation_barplot_95threshold"
fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)
