#!/usr/bin/env python3
"""Replot unconstrained RAD-BEAST trees with the finalized old-tree code style."""

from pathlib import Path
import io
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from Bio import Phylo
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/16_rad_beast/final_unconstrained_60m_summary"
OUT = BASE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PAM = {"AMB", "PAL", "PNG", "PON", "QL", "SAB", "SK", "SUL", "SUR", "TUA", "WJ"}
PSA = {"MDA"}
COLORS = {"PAM": "#D95F02", "PAP": "#1B9E77", "PSA": "#7570B3", "PEQ": "#4D4D4D"}
EPOCHS = [
    ("Holocene", 0.0, 0.0117, "#FFF3B0"),
    ("Pleistocene", 0.0117, 2.58, "#CFE8F3"),
    ("Pliocene", 2.58, 5.333, "#F6C28B"),
    ("Miocene", 5.333, 23.03, "#F4A6A6"),
]


def species(name):
    if name == "PEQ":
        return "PEQ"
    if name in PAM:
        return "PAM"
    if name in PSA:
        return "PSA"
    return "PAP"


def translation(path):
    out, active = {}, False
    for line in open(path):
        if line.strip().lower() == "translate":
            active = True
            continue
        if active:
            match = re.match(r"\s*(\d+)\s+([^,;]+)", line)
            if match:
                out[match.group(1)] = match.group(2).strip()
            if ";" in line:
                break
    return out


def read_annotated_mcc(scenario):
    path = BASE / f"{scenario}_unconstrained_final_mcc.tree"
    trans = translation(path)
    raw = next(line for line in open(path) if line.startswith("tree "))
    newick = re.sub(r"^\[&R\]\s*", "", raw.split("=", 1)[1].strip())
    tree = Phylo.read(io.StringIO(newick), "newick")
    for tip in tree.get_terminals():
        tip.name = trans.get(tip.name, tip.name)
    return tree


def layout(tree):
    tips = tree.get_terminals()
    order = [tip.name for tip in tips]
    y = {name: len(order) - i for i, name in enumerate(order)}

    def ypos(clade):
        if clade.is_terminal():
            return y[clade.name]
        return sum(ypos(child) for child in clade.clades) / len(clade.clades)

    depths = tree.depths()
    height = max(depths[tip] for tip in tips)
    age = {clade: height - depths[clade] for clade in depths}
    return order, y, ypos, age, height


def key(clade):
    return frozenset(t.name for t in clade.get_terminals())


def annotation_value(comment, names, pair=False):
    comment = comment or ""
    for name in names:
        if pair:
            hit = re.search(rf"(?:^|,){re.escape(name)}=\{{([^,}}]+),([^}}]+)\}}", comment)
            if hit:
                return float(hit.group(1)), float(hit.group(2))
        else:
            hit = re.search(rf"(?:^|,){re.escape(name)}=([0-9.eE+-]+)", comment)
            if hit:
                return float(hit.group(1))
    return (float("nan"), float("nan")) if pair else float("nan")


def node_rows(tree):
    _, _, _, ages, _ = layout(tree)
    rows = []
    for i, clade in enumerate(tree.get_nonterminals(order="preorder")):
        median = annotation_value(clade.comment,
            ("CAheight_median", "height_median"))
        lo, hi = annotation_value(clade.comment,
            ("CAheight_95%_HPD", "height_95%_HPD"), pair=True)
        posterior = annotation_value(clade.comment, ("posterior",))
        if pd.isna(median):
            median = ages[clade]
        if pd.isna(lo) or pd.isna(hi):
            lo = hi = ages[clade]
        if pd.isna(posterior):
            posterior = 1.0
        descendants = sorted(t.name for t in clade.get_terminals())
        rows.append({"node": f"N{i+1}", "n_descendants": len(descendants),
            "descendants": ",".join(descendants), "mcc_age": ages[clade],
            "posterior": posterior, "median_age": median,
            "q025_age": lo, "q975_age": hi})
    return pd.DataFrame(rows)


def draw_base(ax, tree, max_age):
    order, y, ypos, age, _ = layout(tree)
    band_y, band_h = len(order) + 0.70, 0.62
    for label, start, end, color in EPOCHS:
        if start >= max_age:
            continue
        end = min(end, max_age)
        ax.add_patch(Rectangle((start, band_y), end - start, band_h,
                               facecolor=color, edgecolor="white", lw=1, zorder=0))
        if end - start >= 0.35:
            ax.text((start + end) / 2, band_y + band_h / 2, label,
                    ha="center", va="center", fontsize=9, weight="bold")

    def recurse(clade):
        if not clade.clades:
            return
        ys = [ypos(child) for child in clade.clades]
        ax.plot([age[clade], age[clade]], [min(ys), max(ys)], color="#222222", lw=1.05)
        for child in clade.clades:
            ax.plot([age[clade], age[child]], [ypos(child), ypos(child)],
                    color="#222222", lw=1.05)
            recurse(child)

    recurse(tree.root)
    for tip in tree.get_terminals():
        ax.text(-max_age * 0.020, y[tip.name], tip.name, ha="left", va="center",
                fontsize=9.2, fontweight="bold", color=COLORS[species(tip.name)])
    ax.set_xlim(max_age * 1.04, -max_age * 0.16)
    ax.set_ylim(0.25, band_y + band_h + 0.55)
    ax.set_yticks([])
    ax.set_xlabel("Time before present (Ma)", fontsize=12)
    ax.grid(axis="x", color="#D9D9D9", lw=0.6, zorder=-2)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.legend(handles=[Patch(facecolor=COLORS[x], label=x)
              for x in ("PAM", "PAP", "PSA", "PEQ")],
              loc="lower left", frameon=False, ncol=4, fontsize=8)
    return ypos, age


def save(fig, stem):
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=350, bbox_inches="tight")
    plt.close(fig)


def plot_ages(scenario, tree, rows):
    lookup = {frozenset(str(r.descendants).split(",")): r for _, r in rows.iterrows()}
    max_age = max(rows.q975_age) * 1.035
    fig, ax = plt.subplots(figsize=(15.5, 13.0))
    ypos, _ = draw_base(ax, tree, max_age)
    for number, clade in enumerate(tree.get_nonterminals(order="preorder")):
        row = lookup[key(clade)]
        yy = ypos(clade)
        ax.plot([row.q025_age, row.q975_age], [yy, yy], color="#4C9CC4", lw=6,
                alpha=0.42, solid_capstyle="butt", zorder=4)
        ax.plot(row.median_age, yy, "o", ms=3.3, color="#155D7A", zorder=5)
        dy = 0.10 if number % 2 == 0 else -0.23
        ax.text(row.q975_age + max_age * 0.007, yy + dy,
                f"{row.median_age:.2f} [{row.q025_age:.2f}–{row.q975_age:.2f}]",
                fontsize=6.3, color="#155D7A", va="center")
    ax.set_title(f"Lineage-unconstrained RAD-seq MCC chronogram — {scenario.capitalize()} calibration\n"
                 "median node age [95% posterior credible interval]", fontsize=14)
    save(fig, f"{scenario}_unconstrained_final_mcc_geological_node_ages_CI")


def plot_support(scenario, tree, rows):
    lookup = {frozenset(str(r.descendants).split(",")): r for _, r in rows.iterrows()}
    max_age = max(rows.q975_age) * 1.035
    fig, ax = plt.subplots(figsize=(14.0, 12.5))
    ypos, age = draw_base(ax, tree, max_age)
    all_tips = key(tree.root)
    for clade in tree.get_nonterminals(order="preorder"):
        descendants = key(clade)
        if descendants == all_tips:
            continue
        ax.text(age[clade] + max_age * 0.007, ypos(clade) + 0.10,
                f"{lookup[descendants].posterior:.2f}", fontsize=7.2,
                color="#7A1F5C", ha="left", va="bottom")
    ax.set_title(f"Lineage-unconstrained RAD-seq MCC topology — {scenario.capitalize()} calibration\n"
                 "node labels: posterior clade support", fontsize=14)
    save(fig, f"{scenario}_unconstrained_final_mcc_geological_node_support")


def cloud_coords(tree, order):
    y = {name: i for i, name in enumerate(order)}
    def ypos(clade):
        if clade.is_terminal():
            return y[clade.name]
        return sum(ypos(child) for child in clade.clades) / len(clade.clades)
    depth = tree.depths()
    height = max(depth[t] for t in tree.get_terminals())
    age = {c: height - depth[c] for c in depth}
    return y, ypos, age, height


def draw_cloud_tree(ax, tree, order, labels=False, color="#367FA8", alpha=.10, lw=1.05):
    y, ypos, age, height = cloud_coords(tree, order)
    def recurse(clade):
        if not clade.clades:
            return
        ys = [ypos(x) for x in clade.clades]
        ax.plot([age[clade], age[clade]], [min(ys), max(ys)], color=color, alpha=alpha, lw=lw)
        for child in clade.clades:
            ax.plot([age[child], age[clade]], [ypos(child), ypos(child)],
                    color=color, alpha=alpha, lw=lw)
            recurse(child)
    recurse(tree.root)
    if labels:
        for name, yy in y.items():
            ax.text(-height*.025, yy, name, color=COLORS[species(name)],
                    va="center", ha="left", fontsize=8, fontweight="bold")
    return height


def plot_cloud(scenario, mcc_tree):
    trees = list(Phylo.parse(BASE / f"{scenario}_unconstrained_100_posterior_full.newick", "newick"))
    order = [q.name for q in mcc_tree.get_terminals()]
    fig, ax = plt.subplots(figsize=(12, 10))
    max_height = 0
    for tree in trees:
        max_height = max(max_height, draw_cloud_tree(ax, tree, order))
    max_height = max(max_height,
        draw_cloud_tree(ax, mcc_tree, order, labels=True, color="black", alpha=1, lw=1.85))
    ax.set_xlim(max_height*1.05, -max_height*.12)
    ax.set_ylim(-1, len(order))
    ax.set_yticks([])
    ax.set_xlabel("Time before present (Ma)")
    ax.set_title(f"Lineage-unconstrained posterior tree cloud — {scenario.capitalize()} calibration")
    ax.grid(axis="x", color="#ddd", lw=.5)
    save(fig, f"{scenario}_unconstrained_posterior_tree_cloud")


for scenario in ("zhang", "givnish"):
    mcc_tree = read_annotated_mcc(scenario)
    rows = node_rows(mcc_tree)
    rows.to_csv(BASE / f"{scenario}_unconstrained_node_ages_support_and_CI.tsv",
                sep="\t", index=False)
    plot_ages(scenario, mcc_tree, rows)
    plot_support(scenario, mcc_tree, rows)
    plot_cloud(scenario, mcc_tree)

