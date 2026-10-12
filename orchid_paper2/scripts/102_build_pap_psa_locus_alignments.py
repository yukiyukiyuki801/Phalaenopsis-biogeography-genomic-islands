#!/usr/bin/env python3
"""Build phased RAD-locus haplotype alignments and balanced distance metrics."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math

import numpy as np
import pandas as pd
import pysam
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/19_local_genealogy"
ALL_VCF = ROOT / "orchid_paper2/05_all_sites_vcf/filtered/orchid_173.paper1_like.biallelic_snps_invariant.vcf.gz"
LOCI = BASE / "inputs/pap_psa_nuclear_loci.tsv"
SAMPLES = [x.strip() for x in (BASE / "inputs/pap_psa_samples.txt").read_text().splitlines() if x.strip()]
RNG = np.random.default_rng(20261002)
N_RESAMPLE = 500
MAX_MISSING = 0.30
MIN_SPECIES_INDIVIDUALS = 5


def allele_for(record, allele_index):
    if allele_index is None:
        return "N"
    alleles = record.alleles
    if allele_index < 0 or allele_index >= len(alleles) or alleles[allele_index] is None:
        return "N"
    allele = alleles[allele_index]
    return allele if len(allele) == 1 else "N"


def pdistance(a, b, minimum):
    comparable = [(x, y) for x, y in zip(a, b) if x != "N" and y != "N"]
    if len(comparable) < minimum:
        return math.nan
    return sum(x != y for x, y in comparable) / len(comparable)


def mean_finite(values):
    values = [x for x in values if math.isfinite(x)]
    return float(np.mean(values)) if values else math.nan


def distance_matrices(seqs, ref):
    names = sorted(seqs)
    encode = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
    array = np.asarray([[encode.get(base.upper(), 4) for base in seqs[name]] for name in names],
                       dtype=np.uint8)
    reference = np.asarray([encode.get(base.upper(), 4) for base in ref], dtype=np.uint8)
    comparable = (array[:, None, :] != 4) & (array[None, :, :] != 4)
    n_comp = comparable.sum(axis=2)
    differences = ((array[:, None, :] != array[None, :, :]) & comparable).sum(axis=2)
    distance = np.divide(differences, n_comp, out=np.full_like(differences, np.nan, dtype=float),
                         where=n_comp >= max(20, math.ceil(.70 * len(ref))))
    ref_comp = (array != 4) & (reference[None, :] != 4)
    ref_n = ref_comp.sum(axis=1)
    ref_diff = ((array != reference[None, :]) & ref_comp).sum(axis=1)
    ref_distance = np.divide(ref_diff, ref_n, out=np.full(len(names), np.nan),
                             where=ref_n >= max(20, math.ceil(.70 * len(ref))))
    name_index = {name: i for i, name in enumerate(names)}
    person_indices = {person: [name_index[f"{person}_h1"], name_index[f"{person}_h2"]]
                      for person in {name.rsplit("_h", 1)[0] for name in names}}
    return distance, ref_distance, person_indices


def one_resample(distance, ref_distance, person_indices, pap_people, psa_people, n):
    chosen_pap = RNG.choice(pap_people, size=n, replace=False)
    chosen_psa = np.asarray(psa_people)
    if len(chosen_psa) > n:
        chosen_psa = RNG.choice(chosen_psa, size=n, replace=False)
    labels, indices = [], []
    for person in list(chosen_pap) + list(chosen_psa):
        species = person[:3]
        labels.extend([species, species])
        indices.extend(person_indices[person])
    dist = distance[np.ix_(indices, indices)]
    within_pap, within_psa, between = [], [], []
    for i in range(len(indices)):
        for j in range(i + 1, len(indices)):
            if not math.isfinite(dist[i, j]):
                continue
            if labels[i] == labels[j] == "PAP":
                within_pap.append(dist[i, j])
            elif labels[i] == labels[j] == "PSA":
                within_psa.append(dist[i, j])
            else:
                between.append(dist[i, j])
    all_dist = [dist[i, j] for i in range(len(indices)) for j in range(i + 1, len(indices))
                if math.isfinite(dist[i, j])]
    ref_dist = [ref_distance[i] for i in indices]
    ref_mean = mean_finite(ref_dist)
    mean_depth = mean_finite(all_dist)
    max_depth = max(all_dist) if all_dist else math.nan

    # Nearest-neighbor mixing: ties contribute their opposite-species fraction.
    mixing = []
    for i in range(len(indices)):
        finite = [(j, dist[i, j]) for j in range(len(indices)) if j != i and math.isfinite(dist[i, j])]
        if not finite:
            continue
        nearest = min(x[1] for x in finite)
        tied = [j for j, value in finite if abs(value - nearest) < 1e-12]
        mixing.append(sum(labels[j] != labels[i] for j in tied) / len(tied))
    cross_species_zero = []
    for i in range(len(indices)):
        opposite = [dist[i, j] for j in range(len(indices))
                    if labels[j] != labels[i] and math.isfinite(dist[i, j])]
        cross_species_zero.append(bool(opposite) and min(opposite) < 1e-12)

    wp, ws, bt = mean_finite(within_pap), mean_finite(within_psa), mean_finite(between)
    within = mean_finite([wp, ws])
    return {
        "n_balanced_individuals_per_species": n,
        "mean_ingroup_pdistance": mean_depth,
        "max_ingroup_pdistance": max_depth,
        "mean_reference_pdistance": ref_mean,
        "relative_mean_depth": mean_depth / ref_mean if ref_mean > 0 else math.nan,
        "relative_max_depth": max_depth / ref_mean if ref_mean > 0 else math.nan,
        "mean_within_pap": wp,
        "mean_within_psa": ws,
        "mean_between": bt,
        "species_separation": bt - within if math.isfinite(bt) and math.isfinite(within) else math.nan,
        "nearest_neighbor_mixing": mean_finite(mixing),
        "cross_species_zero_distance_fraction": mean_finite(cross_species_zero),
    }


def summarize_resamples(rows):
    output = {}
    for key in rows[0]:
        values = np.asarray([x[key] for x in rows], dtype=float)
        values = values[np.isfinite(values)]
        output[f"{key}_median"] = float(np.median(values)) if len(values) else math.nan
        output[f"{key}_q025"] = float(np.quantile(values, .025)) if len(values) else math.nan
        output[f"{key}_q975"] = float(np.quantile(values, .975)) if len(values) else math.nan
    return output


loci = pd.read_csv(LOCI, sep="\t", na_values="NA")
all_vcf = pysam.VariantFile(str(ALL_VCF))
qc_rows, metric_rows, eligible = [], [], []

for index, locus_row in loci.iterrows():
    locus = locus_row.locus_id
    phase_path = BASE / f"phased_loci/{locus}.phased.vcf.gz"
    if not phase_path.exists():
        continue
    phased_vcf = pysam.VariantFile(str(phase_path))
    phased_calls = {}
    for record in phased_vcf:
        for sample in SAMPLES:
            call = record.samples[sample]
            phased_calls[(record.pos, sample)] = (call.get("GT"), call.phased, call.get("PS"))
    phased_vcf.close()

    records = list(all_vcf.fetch(locus_row.chrom, int(locus_row.start)-1, int(locus_row.end)))
    if not records:
        continue
    ref = "".join(record.ref if len(record.ref) == 1 else "N" for record in records)
    sequences, valid_people = {}, defaultdict(list)
    for sample in SAMPLES:
        h1, h2, phase_sets = [], [], []
        heterozygous = phased_heterozygous = unphased_heterozygous = 0
        missing = 0
        for record in records:
            original = record.samples[sample]
            gt = original.get("GT")
            phased = False
            phase_set = None
            if record.alts and record.alts != (".",) and (record.pos, sample) in phased_calls:
                pgt, phased, phase_set = phased_calls[(record.pos, sample)]
                if pgt is not None:
                    gt = pgt
            if gt is None or len(gt) != 2 or None in gt:
                h1.append("N"); h2.append("N"); missing += 1
                continue
            if gt[0] != gt[1]:
                heterozygous += 1
                if phased:
                    phased_heterozygous += 1
                    if phase_set is not None:
                        phase_sets.append(str(phase_set))
                else:
                    unphased_heterozygous += 1
            h1.append(allele_for(record, gt[0])); h2.append(allele_for(record, gt[1]))
        phase_resolved = heterozygous <= 1 or (
            unphased_heterozygous == 0 and len(set(phase_sets)) <= 1)
        missing_fraction = missing / len(records)
        valid = phase_resolved and missing_fraction <= MAX_MISSING
        qc_rows.append({"locus_id": locus, "sample": sample, "species": sample[:3],
                        "alignment_sites": len(records), "missing_sites": missing,
                        "missing_fraction": missing_fraction,
                        "heterozygous_sites": heterozygous,
                        "phased_heterozygous_sites": phased_heterozygous,
                        "unphased_heterozygous_sites": unphased_heterozygous,
                        "phase_sets": len(set(phase_sets)), "phase_resolved": phase_resolved,
                        "valid_for_haplotype_analysis": valid})
        if valid:
            sequences[f"{sample}_h1"] = "".join(h1)
            sequences[f"{sample}_h2"] = "".join(h2)
            valid_people[sample[:3]].append(sample)

    fasta = BASE / f"alignments/{locus}.fa"
    if len(valid_people["PAP"]) >= MIN_SPECIES_INDIVIDUALS and len(valid_people["PSA"]) >= MIN_SPECIES_INDIVIDUALS:
        with open(fasta, "w") as handle:
            handle.write(f">REF_P_equestris\n{ref}\n")
            for name in sorted(sequences):
                handle.write(f">{name}\n{sequences[name]}\n")
        n = min(7, len(valid_people["PAP"]), len(valid_people["PSA"]))
        distance, ref_distance, person_indices = distance_matrices(sequences, ref)
        resamples = [one_resample(distance, ref_distance, person_indices,
                                  valid_people["PAP"], valid_people["PSA"], n)
                     for _ in range(N_RESAMPLE)]
        summary = summarize_resamples(resamples)
        row = {"locus_id": locus, "model_class": locus_row.four_model_class_central25_dxy,
               "chrom": locus_row.chrom, "start": locus_row.start, "end": locus_row.end,
               "alignment_sites": len(records), "informative_snps": locus_row.informative_snps,
               "callable_bases_pair": locus_row.callable_bases,
               "valid_pap_individuals": len(valid_people["PAP"]),
               "valid_psa_individuals": len(valid_people["PSA"]),
               "n_resamples": N_RESAMPLE, **summary}
        metric_rows.append(row); eligible.append(locus)

qc = pd.DataFrame(qc_rows)
metrics = pd.DataFrame(metric_rows)
qc.to_csv(BASE / "results/sample_locus_phase_qc.tsv.gz", sep="\t", index=False,
          compression="gzip")
metrics.to_csv(BASE / "results/locus_balanced_distance_metrics.tsv", sep="\t", index=False,
               na_rep="NA", float_format="%.10g")

# Choose three quality-matched unclassified controls per balancing locus.
eligible_metrics = metrics.dropna(subset=["relative_mean_depth_median"]).copy()
bal = eligible_metrics[eligible_metrics.model_class == "balancing_selection"].copy()
pool = eligible_metrics[eligible_metrics.model_class == "unclassified"].copy()
features = ["alignment_sites", "informative_snps", "valid_pap_individuals", "valid_psa_individuals"]
controls = []
if len(bal) and len(pool):
    combined = pd.concat([bal[features], pool[features]])
    means, stds = combined.mean(), combined.std().replace(0, 1)
    a = (bal[features] - means) / stds
    b = (pool[features] - means) / stds
    distance = cdist(a.to_numpy(), b.to_numpy())
    used = set()
    for i, locus in enumerate(bal.locus_id):
        for j in np.argsort(distance[i]):
            candidate = pool.iloc[j].locus_id
            if candidate not in used:
                controls.append({"balancing_locus": locus, "control_locus": candidate,
                                 "match_distance": distance[i, j]})
                used.add(candidate)
                if sum(x["balancing_locus"] == locus for x in controls) == 3:
                    break
matches = pd.DataFrame(controls)
matches.to_csv(BASE / "results/balancing_matched_unclassified_controls.tsv", sep="\t", index=False)

tree_loci = set(eligible_metrics.loc[eligible_metrics.model_class != "unclassified", "locus_id"])
tree_loci.update(matches.control_locus if len(matches) else [])
manifest = eligible_metrics[eligible_metrics.locus_id.isin(tree_loci)][
    ["locus_id", "model_class", "alignment_sites", "informative_snps",
     "valid_pap_individuals", "valid_psa_individuals"]].sort_values(["model_class", "locus_id"])
manifest["alignment"] = manifest.locus_id.map(lambda x: f"orchid_paper2/19_local_genealogy/alignments/{x}.fa")
manifest.to_csv(BASE / "inputs/iqtree_locus_manifest.tsv", sep="\t", index=False)

with open(BASE / "results/alignment_build_summary.tsv", "w") as handle:
    handle.write("quantity\tvalue\n")
    handle.write(f"phased_locus_vcfs\t{len(set(qc.locus_id))}\n")
    handle.write(f"genealogy_eligible_loci\t{len(metrics)}\n")
    handle.write(f"eligible_balancing_loci\t{sum(metrics.model_class == 'balancing_selection')}\n")
    handle.write(f"iqtree_manifest_loci\t{len(manifest)}\n")
