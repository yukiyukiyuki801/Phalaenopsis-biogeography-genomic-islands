#!/usr/bin/env python3
"""Create population-sorted BA3-SNPs immanc files from the Paper 1 VCF."""

from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_mapping(path: Path):
    original_to_final = {}
    with path.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["primary_ba3_population"]:
                raise ValueError(f"Unassigned locality: {row['original_population']}")
            original_to_final[row["original_population"]] = row["primary_ba3_population"]
    return original_to_final


def read_metadata(path: Path, original_to_final):
    sample_to_pop = {}
    with path.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            sid, locality = row["id"], row["pop"]
            if locality not in original_to_final:
                raise ValueError(f"No BA3 mapping for {sid}: {locality}")
            sample_to_pop[sid] = original_to_final[locality]
    return sample_to_pop


def vcf_keys(path: Path):
    keys = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            keys.add((f[0], f[1], f[3], f[4]))
    return keys


def parse_gt(sample_field, gt_index):
    fields = sample_field.split(":")
    if gt_index >= len(fields):
        return None
    gt = fields[gt_index].replace("|", "/")
    if gt in {"./.", ".", ".//."}:
        return None
    alleles = gt.split("/")
    if len(alleles) != 2 or any(a not in {"0", "1"} for a in alleles):
        return None
    return tuple(map(int, alleles))


def convert(vcf, neutral_vcf, mapping, metadata, outdir, callrate, maf):
    original_to_final = read_mapping(mapping)
    sample_to_pop = read_metadata(metadata, original_to_final)
    neutral = set()

    with vcf.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                break
        else:
            raise ValueError("VCF has no #CHROM header")

    if len(samples) != 170 or len(set(samples)) != 170:
        raise ValueError(f"Expected 170 unique VCF samples, found {len(samples)}")
    missing_meta = sorted(set(samples) - set(sample_to_pop))
    extra_meta = sorted(set(sample_to_pop) - set(samples))
    if missing_meta or extra_meta:
        raise ValueError(f"VCF/metadata mismatch: missing={missing_meta}, extra={extra_meta}")

    ordered_samples = sorted(samples, key=lambda s: (sample_to_pop[s], s))
    vcf_index = {s: i for i, s in enumerate(samples)}
    population_counts = Counter(sample_to_pop[s] for s in ordered_samples)
    if len(population_counts) != 18 or sum(population_counts.values()) != 170:
        raise ValueError(f"Expected 18 populations/170 samples, got {population_counts}")

    datasets = {"all_ld_pruned": []}
    rejected = Counter()
    with vcf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                rejected["not_biallelic_snp"] += 1
                continue
            fmt = f[8].split(":")
            if "GT" not in fmt:
                rejected["no_GT"] += 1
                continue
            gi = fmt.index("GT")
            gts = [parse_gt(f[9 + vcf_index[s]], gi) for s in ordered_samples]
            called = [a for gt in gts if gt is not None for a in gt]
            cr = len(called) / (2 * len(ordered_samples))
            if cr < callrate:
                rejected["callrate"] += 1
                continue
            alt_count = sum(called)
            site_maf = min(alt_count / len(called), 1 - alt_count / len(called))
            if site_maf < maf:
                rejected["maf"] += 1
                continue
            record = ((chrom, pos, ref, alt), gts, cr, site_maf)
            datasets["all_ld_pruned"].append(record)

    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "sample_to_ba3_population.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "ba3_population"])
        w.writerows((s, sample_to_pop[s]) for s in ordered_samples)
    with (outdir / "population_sample_counts.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["ba3_population", "n_individuals"])
        w.writerows(sorted(population_counts.items()))

    summary_rows = []
    for name, records in datasets.items():
        immanc = outdir / f"orchid_18pop_{name}.immanc"
        sites = outdir / f"orchid_18pop_{name}.sites.tsv"
        with immanc.open("w") as out:
            for sample_i, sample in enumerate(ordered_samples):
                pop = sample_to_pop[sample]
                for locus_i, (key, gts, _, _) in enumerate(records, 1):
                    gt = gts[sample_i]
                    a1, a2 = (0, 0) if gt is None else (gt[0] + 1, gt[1] + 1)
                    out.write(f"{sample} {pop} locus{locus_i} {a1} {a2}\n")
        with sites.open("w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["locus", "chrom", "pos", "ref", "alt", "call_rate", "maf"])
            for i, (key, _, cr, site_maf) in enumerate(records, 1):
                w.writerow([f"locus{i}", *key, f"{cr:.6f}", f"{site_maf:.6f}"])
        summary_rows.append((name, len(records), len(ordered_samples), len(population_counts), immanc.stat().st_size))

    with (outdir / "input_preparation_summary.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["dataset", "loci", "individuals", "populations", "immanc_bytes"])
        w.writerows(summary_rows)
    with (outdir / "site_filter_rejections.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["reason", "n_sites"])
        w.writerows(sorted(rejected.items()))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", type=Path, required=True)
    p.add_argument("--neutral-vcf", type=Path, default=None, help=argparse.SUPPRESS)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    # The Paper 1 LD-pruned VCF is already the final quality-filtered SNP set.
    # Defaults therefore add no new call-rate or frequency ascertainment.
    p.add_argument("--callrate", type=float, default=0.0)
    p.add_argument("--maf", type=float, default=0.0)
    a = p.parse_args()
    convert(a.vcf, a.neutral_vcf, a.mapping, a.metadata, a.outdir, a.callrate, a.maf)
