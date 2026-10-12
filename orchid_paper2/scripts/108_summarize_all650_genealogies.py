#!/usr/bin/env python3
"""Summarize complete nuclear RAD genealogies and original 12-contrast permutation family."""
from copy import deepcopy
from pathlib import Path
import math

from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/19_local_genealogy"
RNG = np.random.default_rng(20261002)
N_PERM = 49999


def species(name):
    if name.startswith("PAP"):
        return "PAP"
    if name.startswith("PSA"):
        return "PSA"
    return "REF"


def fitch_score(tree):
    def visit(clade):
        if clade.is_terminal():
            return {species(clade.name)}, 0
        children = [visit(x) for x in clade.clades]
        states, score = children[0]
        score += sum(x[1] for x in children[1:])
        for nxt, _ in children[1:]:
            overlap = states & nxt
            if overlap:
                states = overlap
            else:
                states |= nxt
                score += 1
        return states, score
    return visit(tree.root)[1]


manifest = pd.read_csv(BASE / "inputs/iqtree_all650_diploid_manifest.tsv", sep="\t")
tree_rows = []
for row in manifest.itertuples(index=False):
    path = BASE / f"trees_all650_diploid_ml/{row.locus_id}.treefile"
    if not path.exists():
        continue
    try:
        tree = Phylo.read(path, "newick")
        ref = next(x for x in tree.get_terminals() if x.name == "REF_P_equestris")
    except Exception:
        continue
    tree.root_with_outgroup(ref)
    ingroup = [x for x in tree.get_terminals() if species(x.name) in ("PAP", "PSA")]
    pap = [x for x in ingroup if species(x.name) == "PAP"]
    psa = [x for x in ingroup if species(x.name) == "PSA"]
    if not pap or not psa:
        continue
    mrca = tree.common_ancestor(ingroup)
    height = max(tree.distance(mrca, x) for x in ingroup)
    ref_distance = np.mean([tree.distance(ref, x) for x in ingroup])
    pap_mrca, psa_mrca = tree.common_ancestor(pap), tree.common_ancestor(psa)
    pap_mono = {species(x.name) for x in pap_mrca.get_terminals()} == {"PAP"}
    psa_mono = {species(x.name) for x in psa_mrca.get_terminals()} == {"PSA"}
    copy = deepcopy(tree)
    copy.prune(next(x for x in copy.get_terminals() if x.name == "REF_P_equestris"))
    steps = fitch_score(copy)
    tree_rows.append({"locus_id": row.locus_id, "model_class": row.model_class,
                      "n_ingroup_individuals": len(ingroup),
                      "ml_ingroup_height": height,
                      "ml_reference_distance": ref_distance,
                      "ml_relative_ingroup_height": height/ref_distance if ref_distance > 0 else math.nan,
                      "pap_monophyletic": pap_mono, "psa_monophyletic": psa_mono,
                      "both_species_monophyletic": pap_mono and psa_mono,
                      "species_fitch_steps": steps,
                      "species_fitch_steps_per_tip": steps/len(ingroup)})

tree_metrics = pd.DataFrame(tree_rows)
tree_metrics.to_csv(BASE / "results/all650_diploid_tree_topology_metrics.tsv", sep="\t", index=False, na_rep="NA")
if len(tree_metrics):
    tree_metrics.groupby("model_class").agg(
        n_trees=("locus_id", "size"),
        median_relative_height=("ml_relative_ingroup_height", "median"),
        median_fitch_steps_per_tip=("species_fitch_steps_per_tip", "median"),
        fraction_PAP_monophyletic=("pap_monophyletic", "mean"),
        fraction_PSA_monophyletic=("psa_monophyletic", "mean"),
        fraction_both_species_monophyletic=("both_species_monophyletic", "mean"),
    ).reset_index().to_csv(BASE / "results/all650_diploid_tree_class_summary.tsv", sep="\t", index=False)


def bh(values):
    p = np.asarray(values); order = np.argsort(p); ranked = p[order]
    q = np.minimum.accumulate((ranked*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty_like(q); result[order] = np.minimum(q, 1)
    return result


def compare(data, analysis_name, metric_specs):
    summary, tests = [], []
    classes = ["balancing_selection", "recurrent_selection", "allopatric_selection",
               "divergence_with_gene_flow", "unclassified"]
    for model in classes:
        subset = data[data.model_class == model]
        for metric, _, direction in metric_specs:
            values = subset[metric].dropna().to_numpy()
            if len(values):
                summary.append({"analysis": analysis_name, "model_class": model,
                                "metric": metric, "n_loci": len(values),
                                "median": np.median(values), "mean": np.mean(values)})
    balancing = data[data.model_class == "balancing_selection"]
    for metric, _, direction in metric_specs:
        a = balancing[metric].dropna().to_numpy()
        for comparator in ["unclassified", "recurrent_selection", "allopatric_selection"]:
            b = data.loc[data.model_class == comparator, metric].dropna().to_numpy()
            if not len(a) or not len(b):
                continue
            observed = np.median(a)-np.median(b)
            combined = np.r_[a, b]
            exceed = 0
            for _ in range(N_PERM):
                perm = RNG.permutation(combined)
                delta = np.median(perm[:len(a)])-np.median(perm[len(a):])
                exceed += delta >= observed if direction == "greater" else delta <= observed
            tests.append({"analysis": analysis_name, "metric": metric,
                          "alternative": f"balancing_{direction}",
                          "comparison": f"balancing_vs_{comparator}",
                          "n_balancing": len(a), "n_comparator": len(b),
                          "median_balancing": np.median(a), "median_comparator": np.median(b),
                          "median_difference": observed,
                          "p_one_sided": (exceed+1)/(N_PERM+1),
                          "n_permutations": N_PERM})
    tests = pd.DataFrame(tests)
    if len(tests):
        tests["q_BH_within_analysis"] = bh(tests.p_one_sided)
    return pd.DataFrame(summary), tests


all_metrics = pd.read_csv(BASE / "results/all650_diploid_distance_metrics.tsv", sep="\t", na_values="NA")
all_specs = [("relative_mean_depth_median", "Relative mean depth", "greater"),
             ("relative_max_depth_median", "Relative maximum depth", "greater"),
             ("nearest_neighbor_mixing_median", "PAP–PSA nearest-neighbor mixing", "greater"),
             ("species_separation_median", "Between-minus-within distance", "less")]
all_summary, all_tests = compare(all_metrics, "all_650_phase_free", all_specs)


all_summary.to_csv(BASE / "results/genealogy_class_summary.tsv",sep="\t",index=False)
all_tests.to_csv(BASE / "results/genealogy_class_permutation_tests.tsv",sep="\t",index=False)
