#!/usr/bin/env python3
"""Build individual-level, ambiguity-coded alignments for every PAP--PSA RAD locus.

No locus is filtered here. Heterozygous diploid genotypes are represented by IUPAC
codes, missing genotypes by N, and the P. equestris reference sequence is included
as an explicit outgroup. Expected allele-mismatch distances provide a phase-free
genealogy-depth sensitivity analysis for all loci.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math

import numpy as np
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/19_local_genealogy"
VCF_PATH = ROOT / "orchid_paper2/05_all_sites_vcf/filtered/orchid_173.paper1_like.biallelic_snps_invariant.vcf.gz"
LOCI_PATH = BASE / "inputs/pap_psa_nuclear_loci.tsv"
SAMPLE_PATH = BASE / "inputs/pap_psa_samples.txt"
OUTDIR = BASE / "alignments_all650_diploid"
RNG = np.random.default_rng(20261002)
N_RESAMPLE = 500

IUPAC = {
    frozenset(("A", "G")): "R", frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S", frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K", frozenset(("A", "C")): "M",
}


def base(record, index):
    if index is None or index < 0 or index >= len(record.alleles):
        return None
    value = record.alleles[index]
    return value.upper() if value is not None and len(value) == 1 else None


def consensus(record, gt):
    if gt is None or len(gt) != 2 or None in gt:
        return "N"
    alleles = [base(record, x) for x in gt]
    if None in alleles:
        return "N"
    if alleles[0] == alleles[1]:
        return alleles[0]
    return IUPAC.get(frozenset(alleles), "N")


def mean_finite(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if len(values) else math.nan


def summarize_locus(genotypes, species, reference):
    """Repeated equal-n distance summaries using expected allele mismatch.

    genotypes is sample x site alternate-allele dosage (0, .5, 1), with NaN for
    missing. The expected difference between random alleles is
    p_i(1-p_j)+(1-p_i)p_j. Invariant sites are retained.
    """
    names = list(genotypes)
    matrix = np.asarray([genotypes[x] for x in names], dtype=float)
    n = len(names)
    dist = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            ok = np.isfinite(matrix[i]) & np.isfinite(matrix[j])
            if ok.sum() < max(20, math.ceil(.70 * matrix.shape[1])):
                continue
            pi, pj = matrix[i, ok], matrix[j, ok]
            d = np.mean(pi * (1 - pj) + (1 - pi) * pj)
            dist[i, j] = dist[j, i] = d
    refdist = np.full(n, np.nan)
    for i in range(n):
        ok = np.isfinite(matrix[i]) & np.isfinite(reference)
        if ok.sum() >= max(20, math.ceil(.70 * matrix.shape[1])):
            # reference is encoded as 0/1 relative to each record's REF/ALT state
            refdist[i] = np.mean(matrix[i, ok] * (1-reference[ok]) +
                                 (1-matrix[i, ok]) * reference[ok])
    index = {x: i for i, x in enumerate(names)}
    pap = [x for x in names if species[x] == "PAP" and np.isfinite(dist[index[x]]).sum()]
    psa = [x for x in names if species[x] == "PSA" and np.isfinite(dist[index[x]]).sum()]
    if not pap or not psa:
        return {}
    take = min(7, len(pap), len(psa))
    rows = []
    for _ in range(N_RESAMPLE):
        a = list(RNG.choice(pap, take, replace=False))
        b = list(RNG.choice(psa, take, replace=False))
        chosen = a + b
        ids = [index[x] for x in chosen]
        within_a, within_b, between, all_d = [], [], [], []
        for ii in range(len(ids)):
            for jj in range(ii + 1, len(ids)):
                value = dist[ids[ii], ids[jj]]
                if not np.isfinite(value):
                    continue
                all_d.append(value)
                if ii < take and jj < take:
                    within_a.append(value)
                elif ii >= take and jj >= take:
                    within_b.append(value)
                else:
                    between.append(value)
        nearest_mix = []
        for ii, idx in enumerate(ids):
            choices = [(jj, dist[idx, jdx]) for jj, jdx in enumerate(ids)
                       if jj != ii and np.isfinite(dist[idx, jdx])]
            if not choices:
                continue
            minimum = min(x[1] for x in choices)
            ties = [jj for jj, value in choices if abs(value-minimum) < 1e-12]
            nearest_mix.append(np.mean([(jj < take) != (ii < take) for jj in ties]))
        wp, ws, bt = mean_finite(within_a), mean_finite(within_b), mean_finite(between)
        depth, rdepth = mean_finite(all_d), mean_finite(refdist[ids])
        rows.append({
            "mean_ingroup_distance": depth,
            "max_ingroup_distance": max(all_d) if all_d else math.nan,
            "mean_reference_distance": rdepth,
            "relative_mean_depth": depth/rdepth if rdepth > 0 else math.nan,
            "relative_max_depth": max(all_d)/rdepth if all_d and rdepth > 0 else math.nan,
            "mean_within_pap": wp, "mean_within_psa": ws, "mean_between": bt,
            "species_separation": bt-mean_finite([wp, ws]) if np.isfinite(bt) else math.nan,
            "nearest_neighbor_mixing": mean_finite(nearest_mix),
        })
    output = {"balanced_individuals_per_species": take, "n_resamples": N_RESAMPLE}
    for key in rows[0]:
        values = np.asarray([x[key] for x in rows], dtype=float)
        values = values[np.isfinite(values)]
        output[f"{key}_median"] = float(np.median(values)) if len(values) else math.nan
        output[f"{key}_q025"] = float(np.quantile(values, .025)) if len(values) else math.nan
        output[f"{key}_q975"] = float(np.quantile(values, .975)) if len(values) else math.nan
    return output


OUTDIR.mkdir(parents=True, exist_ok=True)
(BASE / "results").mkdir(exist_ok=True)
loci = pd.read_csv(LOCI_PATH, sep="\t", na_values="NA")
samples = [x.strip() for x in SAMPLE_PATH.read_text().splitlines() if x.strip()]
species = {x: x[:3] for x in samples}
vcf = pysam.VariantFile(str(VCF_PATH))
manifest, qc, metrics = [], [], []

for row in loci.itertuples(index=False):
    records = list(vcf.fetch(row.chrom, int(row.start)-1, int(row.end)))
    sequences = {x: [] for x in samples}
    dosage = {x: [] for x in samples}
    reference_dosage = []
    variable_sites = 0
    for record in records:
        ref = record.ref.upper() if len(record.ref) == 1 else "N"
        is_snp = bool(record.alts) and record.alts != (".",) and len(record.alts) == 1 and len(record.alts[0]) == 1
        variable_sites += int(is_snp)
        reference_dosage.append(0.0 if ref != "N" else np.nan)
        for sample in samples:
            gt = record.samples[sample].get("GT")
            sequences[sample].append(consensus(record, gt))
            if gt is None or len(gt) != 2 or None in gt or not is_snp or any(x not in (0, 1) for x in gt):
                # invariant callable records are homozygous reference by construction
                if not is_snp and gt is not None and len(gt) == 2 and None not in gt:
                    dosage[sample].append(0.0)
                else:
                    dosage[sample].append(np.nan)
            else:
                dosage[sample].append(sum(gt)/2)
    fasta = OUTDIR / f"{row.locus_id}.fa"
    with fasta.open("w") as handle:
        handle.write(">REF_P_equestris\n" + "".join(r.ref.upper() if len(r.ref) == 1 else "N" for r in records) + "\n")
        for sample in samples:
            handle.write(f">{sample}\n{''.join(sequences[sample])}\n")
            missing = sequences[sample].count("N")
            qc.append({"locus_id": row.locus_id, "sample": sample, "species": species[sample],
                       "alignment_sites": len(records), "missing_sites": missing,
                       "missing_fraction": missing/len(records) if records else 1.0})
    manifest.append({"locus_id": row.locus_id,
                     "model_class": row.four_model_class_central25_dxy,
                     "alignment_sites": len(records), "variable_sites": variable_sites,
                     "informative_snps": row.informative_snps,
                     "alignment": str(fasta.relative_to(ROOT))})
    summary = summarize_locus(dosage, species, np.asarray(reference_dosage, dtype=float)) if records else {}
    metrics.append({"locus_id": row.locus_id,
                    "model_class": row.four_model_class_central25_dxy,
                    "chrom": row.chrom, "start": row.start, "end": row.end,
                    "alignment_sites": len(records), "variable_sites": variable_sites,
                    "informative_snps": row.informative_snps, **summary})

pd.DataFrame(manifest).to_csv(BASE / "inputs/iqtree_all650_diploid_manifest.tsv", sep="\t", index=False)
pd.DataFrame(qc).to_csv(BASE / "results/all650_diploid_sample_locus_qc.tsv.gz", sep="\t", index=False, compression="gzip")
pd.DataFrame(metrics).to_csv(BASE / "results/all650_diploid_distance_metrics.tsv", sep="\t", index=False, na_rep="NA")
with (BASE / "results/all650_alignment_build_summary.tsv").open("w") as handle:
    handle.write("quantity\tvalue\n")
    handle.write(f"input_loci\t{len(loci)}\n")
    handle.write(f"alignments_written\t{len(manifest)}\n")
    handle.write(f"PAP_samples\t{sum(x.startswith('PAP') for x in samples)}\n")
    handle.write(f"PSA_samples\t{sum(x.startswith('PSA') for x in samples)}\n")
    handle.write("outgroup\tREF_P_equestris\n")

