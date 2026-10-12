#!/usr/bin/env python3
"""Plot honest phased RAD genealogies with identical haplotypes collapsed."""
from collections import defaultdict
from pathlib import Path

from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/19_local_genealogy"
RESULTS, PLOTS = BASE / "results", BASE / "plots"
TREE_DIR = BASE / "trees_retained501_phased_ml"
COLLAPSED_DIR = BASE / "representative_collapsed_phased_trees"
COLLAPSED_DIR.mkdir(exist_ok=True)

ORDER = ["balancing_selection", "recurrent_selection", "allopatric_selection",
         "divergence_with_gene_flow", "unclassified"]
LABEL = {"balancing_selection":"Balancing selection", "recurrent_selection":"Recurrent selection",
         "allopatric_selection":"Allopatric selection",
         "divergence_with_gene_flow":"Divergence with gene flow",
         "unclassified":"Unclassified background"}
CLASS_COLOR = {"balancing_selection":"#7B3294", "recurrent_selection":"#D73027",
               "allopatric_selection":"#4575B4", "divergence_with_gene_flow":"#1B9E77",
               "unclassified":"#777777"}
TIP_COLOR = {"PAP":"#E66101", "PSA":"#1B9E77", "shared":"#7B3294", "REF":"#202020"}
METRICS = ["relative_mean_depth_median", "nearest_neighbor_mixing_median",
           "species_separation_median"]


def read_fasta(path):
    result, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]; result[name] = ""
        else:
            result[name] += line.strip().upper()
    return result


def species(name):
    return "PAP" if name.startswith("PAP") else ("PSA" if name.startswith("PSA") else "REF")


# Reuse the five representatives actually illustrated in the manuscript;
# do not rerun or export the unreported phased-subset sensitivity analysis.
selected = pd.read_csv(RESULTS / "representative_rad_tree_loci.tsv", sep="\t")


def collapse_identical(row):
    locus = row.locus_id
    seqs = read_fasta(BASE / f"alignments/{locus}.fa")
    groups = defaultdict(list)
    for name, seq in seqs.items():
        if name != "REF_P_equestris":
            groups[seq].append(name)
    tree = Phylo.read(TREE_DIR / f"{locus}.treefile", "newick")
    ref = next(x for x in tree.get_terminals() if x.name == "REF_P_equestris")
    tree.root_with_outgroup(ref)
    group_rows, rename = [], {}
    ordered_groups = sorted(groups.values(), key=lambda names: sorted(names)[0])
    for index, names in enumerate(ordered_groups, 1):
        names = sorted(names); keep = names[0]
        for drop in names[1:]:
            tip = next((x for x in tree.get_terminals() if x.name == drop), None)
            if tip is not None:
                tree.prune(tip)
        pap = sum(x.startswith("PAP") for x in names); psa = sum(x.startswith("PSA") for x in names)
        state = "shared" if pap and psa else ("PAP" if pap else "PSA")
        node = f"H{index:03d}"
        rename[keep] = node
        group_rows.append({"locus_id": locus, "model_class": row.model_class,
                           "node": node, "representative_tip": keep,
                           "n_haplotypes": len(names), "n_PAP": pap, "n_PSA": psa,
                           "composition": state, "members": ",".join(names)})
    for tip in tree.get_terminals():
        if tip.name in rename:
            tip.name = rename[tip.name]
    # IQ-TREE uses ~1e-6 as its minimum numerical branch. Treat internal
    # branches at that floor as unresolved and display them as polytomies.
    tree.collapse_all(lambda c: (not c.is_terminal()) and
                      c.branch_length is not None and c.branch_length <= 1.01e-6)
    ref = next(x for x in tree.get_terminals() if x.name == "REF_P_equestris")
    tree.root_with_outgroup(ref)
    Phylo.write(tree, COLLAPSED_DIR / f"{locus}.collapsed_identical.treefile", "newick")
    return tree, pd.DataFrame(group_rows)


def coords(tree, normalizer):
    tips = tree.get_terminals(); y = {tip:i for i, tip in enumerate(tips)}; x = {tree.root:0.0}
    def walk_x(clade):
        for child in clade.clades:
            x[child] = x[clade] + (child.branch_length or 0)/normalizer
            walk_x(child)
    def walk_y(clade):
        if clade in y: return y[clade]
        values = [walk_y(c) for c in clade.clades]
        y[clade] = (min(values)+max(values))/2
        return y[clade]
    walk_x(tree.root); walk_y(tree.root)
    return x, y, tips


panels, membership = [], []
for row in selected.itertuples(index=False):
    tree, groups = collapse_identical(row)
    membership.append(groups)
    # The reference has already oriented the tree. Remove it from the display
    # so the plotted genealogy represents PAP--PSA only.
    ref = next(x for x in tree.get_terminals() if x.name == "REF_P_equestris")
    tree.prune(ref)
    while len(tree.root.clades) == 1:
        tree.root = tree.root.clades[0]
        tree.root.branch_length = 0.0
    x, y, tips = coords(tree, 1.0)
    metadata = groups.set_index("node").to_dict("index")
    panels.append((row, tree, x, y, tips, metadata))

members = pd.concat(membership, ignore_index=True)
members.to_csv(RESULTS / "representative_collapsed_haplotype_nodes.tsv", sep="\t", index=False)
selection_rows = []
for row, tree, x, y, tips, metadata in panels:
    selection_rows.append({"locus_id":row.locus_id, "model_class":row.model_class,
                           "representative_score":row.representative_score,
                           "alignment_sites":row.alignment_sites,
                           "informative_snps":row.informative_snps,
                           "relative_mean_depth_median":row.relative_mean_depth_median,
                           "nearest_neighbor_mixing_median":row.nearest_neighbor_mixing_median,
                           "species_separation_median":row.species_separation_median,
                           "phased_haplotypes":sum(v["n_haplotypes"] for v in metadata.values()),
                           "unique_haplotype_nodes":len(metadata),
                           "shared_exact_haplotype_nodes":sum(v["composition"]=="shared" for v in metadata.values())})
pd.DataFrame(selection_rows).to_csv(RESULTS / "representative_rad_tree_loci.tsv", sep="\t", index=False)

fig, axes = plt.subplots(5, 1, figsize=(15, 17), sharex=True)
fig.subplots_adjust(hspace=.38)
depth_limit = selected.relative_mean_depth_median.max()*1.08

for ax, (row, tree, x, y, tips, metadata) in zip(axes, panels):
    # Coalescent-style aligned-tip projection: every sampled terminal is at
    # the present (depth=0). Root depth is the independently calculated,
    # reference-normalized mean pairwise divergence, allowing roots to differ
    # among loci on one common scale. Internal positions follow their fraction
    # of maximum remaining descendant ML depth.
    remaining = {}
    def descendant_height(clade):
        if clade.is_terminal():
            remaining[clade] = 0.0
        else:
            remaining[clade] = max((child.branch_length or 0.0) + descendant_height(child)
                                   for child in clade.clades)
        return remaining[clade]
    root_height = descendant_height(tree.root)
    locus_depth = row.relative_mean_depth_median
    sx = {clade: locus_depth*remaining[clade]/root_height for clade in tree.find_clades()}
    for parent in tree.find_clades(order="preorder"):
        if not parent.clades: continue
        cy = [y[c] for c in parent.clades]
        ax.plot([sx[parent],sx[parent]],[min(cy),max(cy)],color="#555555",lw=.65)
        for child in parent.clades:
            ax.plot([sx[parent],sx[child]],[y[child],y[child]],color="#555555",lw=.65)
    for tip in tips:
        state, count = metadata[tip.name]["composition"], metadata[tip.name]["n_haplotypes"]
        ax.scatter(sx[tip], y[tip], s=13+8*np.sqrt(count), c=TIP_COLOR[state],
                   edgecolor="white", linewidth=.25, zorder=4)
        if count >= 10:
            ax.text(sx[tip]+.008, y[tip], f"{count}", fontsize=7, va="center")
    ax.set_xlim(depth_limit, -.015); ax.set_ylim(-1,len(tips)); ax.invert_yaxis()
    ax.set_yticks([]); ax.set_xticks([]); ax.set_xlabel("")
    ax.set_title(f"{LABEL[row.model_class]} — {row.locus_id}  |  "
                 f"depth={locus_depth:.3f}; {len(metadata)} unique nodes from "
                 f"{sum(v['n_haplotypes'] for v in metadata.values())} phased haplotypes",
                 loc="left",fontsize=12,fontweight="bold")
    for side in ("top","right","left","bottom"): ax.spines[side].set_visible(False)

legend = [Line2D([0],[0],marker="o",color="none",markerfacecolor=TIP_COLOR[x],markeredgecolor="none",
                 markersize=8,label={"PAP":"PAP-only haplotype","PSA":"PSA-only haplotype",
                                     "shared":"identical PAP + PSA haplotype"}[x])
          for x in ("PAP","PSA","shared")]
fig.legend(handles=legend,loc="lower left",bbox_to_anchor=(.055,.045),ncol=4,frameon=False,fontsize=10.5)
fig.text(.99,.012,"Circle area indicates the number of identical phased haplotypes; labels show groups with n>=10.\n"
         "All PAP/PSA tips align at the present (0); roots vary by reference-normalized depth on one shared scale. "
         "The reference was used to root the trees and then removed from display.",ha="right",va="bottom",fontsize=9.1)
fig.suptitle("Representative PAP–PSA phased RAD genealogies",
             fontsize=17,fontweight="bold",y=.995)
fig.subplots_adjust(top=.965,bottom=.115,left=.04,right=.985)
stem=PLOTS/"pap_psa_representative_rad_trees_all_classes"
fig.savefig(f"{stem}.png",dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig(f"{stem}.pdf",bbox_inches="tight",facecolor="white")
plt.close(fig)

# Alternative diagonal/fan rendering: the data, coordinates, and node symbols
# are identical to the rectangular figure, but parent--child connections are
# single straight lines rather than right-angle elbows.
fig, axes = plt.subplots(5, 1, figsize=(15, 17), sharex=True)
fig.subplots_adjust(hspace=.38)
for ax, (row, tree, x, y, tips, metadata) in zip(axes, panels):
    remaining = {}
    def descendant_height_diagonal(clade):
        if clade.is_terminal():
            remaining[clade] = 0.0
        else:
            remaining[clade] = max((child.branch_length or 0.0) +
                                   descendant_height_diagonal(child)
                                   for child in clade.clades)
        return remaining[clade]
    root_height = descendant_height_diagonal(tree.root)
    locus_depth = row.relative_mean_depth_median
    sx = {clade: locus_depth*remaining[clade]/root_height for clade in tree.find_clades()}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            ax.plot([sx[parent],sx[child]],[y[parent],y[child]],color="#555555",lw=.65)
    for tip in tips:
        state, count = metadata[tip.name]["composition"], metadata[tip.name]["n_haplotypes"]
        ax.scatter(sx[tip],y[tip],s=13+8*np.sqrt(count),c=TIP_COLOR[state],
                   edgecolor="white",linewidth=.25,zorder=4)
        if count >= 10:
            ax.text(sx[tip]+.008,y[tip],f"{count}",fontsize=7,va="center")
    ax.set_xlim(depth_limit,-.015); ax.set_ylim(-1,len(tips)); ax.invert_yaxis()
    ax.set_yticks([]); ax.set_xticks([]); ax.set_xlabel("")
    ax.set_title(f"{LABEL[row.model_class]} — {row.locus_id}  |  "
                 f"depth={locus_depth:.3f}; {len(metadata)} unique nodes from "
                 f"{sum(v['n_haplotypes'] for v in metadata.values())} phased haplotypes",
                 loc="left",fontsize=12,fontweight="bold")
    for side in ("top","right","left","bottom"):
        ax.spines[side].set_visible(False)
fig.legend(handles=legend,loc="lower left",bbox_to_anchor=(.055,.045),ncol=3,
           frameon=False,fontsize=10.5)
fig.text(.99,.012,"Circle area indicates the number of identical phased haplotypes; labels show groups with n>=10.\n"
         "Straight parent--child lines provide the diagonal/triangular rendering; tree coordinates match the rectangular version.",
         ha="right",va="bottom",fontsize=9.1)
fig.suptitle("Representative PAP–PSA phased RAD genealogies — diagonal branches",
             fontsize=17,fontweight="bold",y=.995)
fig.subplots_adjust(top=.965,bottom=.115,left=.04,right=.985)
diagonal_stem=PLOTS/"pap_psa_representative_rad_trees_all_classes_diagonal"
fig.savefig(f"{diagonal_stem}.png",dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig(f"{diagonal_stem}.pdf",bbox_inches="tight",facecolor="white")
plt.close(fig)
