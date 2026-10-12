#!/usr/bin/env python3
"""Plot combined-chain MCMCTree node ages/95% CI against geological time."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("orchid_paper2")
DATING = ROOT / "15_ortholog_phylogenomics/single_copy_dating/dating"
OUT = ROOT / "15_ortholog_phylogenomics/single_copy_dating/plots"
OUT.mkdir(parents=True, exist_ok=True)

TIP_LABEL = {
    "RICE": "Oryza sativa", "ASPA": "Asparagus officinalis",
    "APO": "Apostasia shenzhenica", "VPLA": "Vanilla planifolia",
    "PGUA": "Platanthera guangdongensis", "DCAT": "Dendrobium catenatum",
    "DNOB": "Dendrobium nobile", "PEQ": "Phalaenopsis equestris",
    "PAP": "Phalaenopsis aphrodite", "PAM": "Phalaenopsis amabilis",
}
TIPS = ["RICE", "ASPA", "APO", "VPLA", "PGUA", "DCAT", "DNOB", "PEQ", "PAP", "PAM"]
# MCMCTree node IDs follow the labeled fixed input topology.
CHILDREN = {
    "n11": ["RICE", "n12"], "n12": ["ASPA", "n13"],
    "n13": ["APO", "n14"], "n14": ["VPLA", "n15"],
    "n15": ["PGUA", "n16"], "n16": ["n17", "n18"],
    "n17": ["DCAT", "DNOB"], "n18": ["n19", "PEQ"],
    "n19": ["PAP", "PAM"],
}
PARAM = {f"n{i}": f"t_n{i}" for i in range(11, 20)}

PERIODS = [
    ("Quaternary", 0, 2.58, "#F6E8A6"),
    ("Neogene", 2.58, 23.03, "#F2B56B"),
    ("Paleogene", 23.03, 66.0, "#E99D8F"),
    ("Cretaceous", 66.0, 145.0, "#8FCB82"),
    ("Jurassic", 145.0, 201.4, "#6FC5C8"),
]

def summaries(scenario):
    frames = []
    for chain in ("chain1", "chain2"):
        frames.append(pd.read_csv(DATING/scenario/chain/"mcmc.txt", sep="\t"))
    x = pd.concat(frames, ignore_index=True)
    rows = []
    for node, parameter in PARAM.items():
        z = x[parameter].to_numpy() * 100.0  # time unit in ctl is 100 Ma
        rows.append({"node": node, "parameter": parameter,
                     "median_ma": np.median(z),
                     "ci_low_ma": np.quantile(z, .025),
                     "ci_high_ma": np.quantile(z, .975),
                     "mean_ma": np.mean(z)})
    return pd.DataFrame(rows).set_index("node")

def descendant_tips(node):
    if node in TIPS:
        return [node]
    return sum((descendant_tips(c) for c in CHILDREN[node]), [])

def draw(scenario, title):
    s = summaries(scenario)
    y = {tip: len(TIPS)-i for i, tip in enumerate(TIPS)}
    def set_y(node):
        if node in y: return y[node]
        yy = [set_y(c) for c in CHILDREN[node]]
        y[node] = sum(yy)/len(yy)
        return y[node]
    set_y("n11")
    age = {tip: 0.0 for tip in TIPS}
    age.update({node: s.loc[node, "median_ma"] for node in CHILDREN})
    max_age = max(205.0, s.loc["n11", "ci_high_ma"] * 1.04)

    fig, ax = plt.subplots(figsize=(15, 8.7))
    band_y, band_h = len(TIPS)+0.70, 0.55
    for label, start, end, color in PERIODS:
        if start >= max_age: continue
        end2 = min(end, max_age)
        ax.add_patch(Rectangle((start, band_y), end2-start, band_h,
                               facecolor=color, edgecolor="white", lw=1, zorder=0))
        if end2-start > 8:
            ax.text((start+end2)/2, band_y+band_h/2, label,
                    ha="center", va="center", fontsize=10, weight="bold")

    for node, children in CHILDREN.items():
        ys = [y[c] for c in children]
        ax.plot([age[node], age[node]], [min(ys), max(ys)], color="black", lw=1.15)
        for child in children:
            ax.plot([age[node], age[child]], [y[child], y[child]], color="black", lw=1.15)

    # Posterior 95% intervals and age labels for every internal node.
    for node in CHILDREN:
        row = s.loc[node]
        ax.plot([row.ci_low_ma, row.ci_high_ma], [y[node], y[node]],
                color="#4C9CC4", lw=7, alpha=.38, solid_capstyle="butt", zorder=3)
        ax.plot(row.median_ma, y[node], "o", ms=3.6, color="#155D7A", zorder=4)
        ax.text(row.ci_high_ma + max_age*.008, y[node]+.10,
                f"{row.median_ma:.1f} [{row.ci_low_ma:.1f}–{row.ci_high_ma:.1f}]",
                fontsize=7.5, va="bottom")

    for tip in TIPS:
        ax.text(-max_age*.018, y[tip], TIP_LABEL[tip], ha="left", va="center",
                fontsize=10.5, fontstyle="italic")

    ax.set_xlim(max_age, -max_age*.27)
    ax.set_ylim(.35, band_y+band_h+.55)
    ax.set_yticks([])
    ax.set_xlabel("Time before present (Ma)", fontsize=12)
    ax.set_title(title + "\nmedian node age [95% posterior credible interval]", fontsize=14)
    ax.grid(axis="x", color="#D9D9D9", lw=.6, zorder=-2)
    for side in ("left", "right", "top"): ax.spines[side].set_visible(False)
    fig.tight_layout()
    stem = OUT / f"single_copy_mcmctree_{scenario}_geological_time"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    table = s.reset_index()
    table["descendant_tips"] = table.node.map(lambda n: ",".join(TIP_LABEL[t] for t in descendant_tips(n)))
    table.to_csv(OUT/f"single_copy_mcmctree_{scenario}_node_ages.tsv", sep="\t", index=False)

for scenario, title in [
    ("zhang_multipoint", "Single-copy nuclear chronogram — Zhang calibration set"),
    ("givnish_multipoint", "Single-copy nuclear chronogram — Givnish calibration set"),
]:
    draw(scenario, title)
