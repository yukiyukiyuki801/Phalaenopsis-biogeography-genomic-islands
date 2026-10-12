#!/usr/bin/env python3
"""Prepare common BAM lists, callable positions, and RAD-locus mapping for ANGSD.

The callable target is deliberately defined by the already quality-filtered
all-sites VCF.  This keeps variant and invariant positions, while excluding
reference positions that failed the original cohort calling/filtering.
"""
from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path


GROUPS = ("PAM", "PAP", "PSA")
GROUP_N = {"PAM": 88, "PAP": 78, "PSA": 7}
PAIRS = (("PAM", "PAP"), ("PAM", "PSA"), ("PAP", "PSA"))


def read_rad_intervals(path: Path, max_gap: int = 500):
    by_contig = defaultdict(list)
    locus = 0
    previous_contig = None
    previous_end = None
    with path.open() as handle:
        for raw in handle:
            contig, start, end = raw.rstrip().split("\t")[:3]
            start, end = int(start), int(end)
            if contig != previous_contig or previous_end is None or start - previous_end > max_gap:
                locus += 1
            by_contig[contig].append((start + 1, end, locus))
            previous_contig, previous_end = contig, end
    return by_contig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    out = root / "21_angsd_rad_landscapes" / "inputs"
    out.mkdir(parents=True, exist_ok=True)

    manifest = root / "00_metadata" / "run_manifest.tsv"
    bam_dir = root / "04_bam"
    samples = []
    with manifest.open() as handle:
        header = handle.readline().rstrip().split("\t")
        sample_col = header.index("sample_id")
        for raw in handle:
            fields = raw.rstrip().split("\t")
            if fields:
                samples.append(fields[sample_col])

    observed = []
    for group in GROUPS:
        group_samples = [sample for sample in samples if sample.startswith(group)]
        with (out / f"{group}.bamlist").open("w") as handle:
            for sample in group_samples:
                bam = bam_dir / f"{sample}.sorted.bam"
                bai = Path(str(bam) + ".bai")
                if not bam.is_file() or not bai.is_file():
                    raise FileNotFoundError(f"Missing BAM or index for {sample}: {bam}")
                handle.write(f"{bam}\n")
                observed.append(str(bam))
    with (out / "ALL.bamlist").open("w") as handle:
        handle.write("\n".join(observed) + "\n")

    intervals = read_rad_intervals(root / "05_all_sites_vcf" / "callable" / "cohort_callable.bed")
    cursors = defaultdict(int)
    vcf = root / "05_all_sites_vcf" / "filtered" / "orchid_173.paper1_like.biallelic_snps_invariant.vcf.gz"
    n_sites = n_variant = n_invariant = 0
    locus_counts = defaultdict(lambda: [None, None, None, 0, 0])
    pair_snp_counts = defaultdict(int)
    with gzip.open(vcf, "rt") as source, \
            (out / "callable_sites.txt").open("w") as sites, \
            (out / "callable_sites_alleles.txt").open("w") as allele_sites, \
            (out / "site_to_rad_locus.tsv").open("w") as mapping:
        mapping.write("chrom\tpos\tlocus_id\tvcf_variant\n")
        for raw in source:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip().split("\t", 8)
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            # Reference-block records can carry a multi-base REF with ALT='.'.
            # They are not single nucleotide sites and cannot be represented
            # in ANGSD's forced major/minor site format.
            if ref not in "ACGT" or len(ref) != 1:
                continue
            if alt != "." and (alt not in "ACGT" or len(alt) != 1):
                continue
            iv = intervals[chrom]
            cursor = cursors[chrom]
            while cursor < len(iv) and pos > iv[cursor][1]:
                cursor += 1
            cursors[chrom] = cursor
            if cursor == len(iv) or pos < iv[cursor][0]:
                raise RuntimeError(f"VCF position lies outside callable BED: {chrom}:{pos}")
            _, _, locus_num = iv[cursor]
            locus_id = f"RAD{locus_num:05d}"
            is_variant = int(alt != ".")
            # dXY requires the same allele coding in every population.  Use
            # the cohort VCF's biallelic REF/ALT definition.  At invariant
            # sites, assign a deterministic non-reference minor allele; its
            # estimated frequency should be approximately zero.
            forced_minor = alt if is_variant else next(base for base in "ACGT" if base != ref)
            sites.write(f"{chrom}\t{pos}\n")
            allele_sites.write(f"{chrom}\t{pos}\t{ref}\t{forced_minor}\n")
            mapping.write(f"{chrom}\t{pos}\t{locus_id}\t{is_variant}\n")
            record = locus_counts[locus_id]
            record[0] = chrom
            record[1] = pos if record[1] is None else min(record[1], pos)
            record[2] = pos if record[2] is None else max(record[2], pos)
            record[3] += 1
            record[4] += is_variant
            if is_variant:
                info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
                for group1, group2 in PAIRS:
                    an1 = int(info.get(f"AN_{group1}", 0))
                    an2 = int(info.get(f"AN_{group2}", 0))
                    ac1 = int(info.get(f"AC_{group1}", 0).split(",")[0])
                    ac2 = int(info.get(f"AC_{group2}", 0).split(",")[0])
                    if (an1 >= 2 * GROUP_N[group1] * 0.70 and
                            an2 >= 2 * GROUP_N[group2] * 0.70 and
                            0 < ac1 + ac2 < an1 + an2):
                        pair_snp_counts[(f"{group1}-{group2}", locus_id)] += 1
            n_sites += 1
            n_variant += is_variant
            n_invariant += 1 - is_variant

    with (out / "rad_locus_target_summary.tsv").open("w") as handle:
        handle.write("locus_id\tchrom\tfirst_target_pos\tlast_target_pos\ttarget_sites\tvcf_variant_sites\n")
        for locus_id in sorted(locus_counts, key=lambda x: int(x[3:])):
            chrom, first, last, target_n, variant_n = locus_counts[locus_id]
            handle.write(f"{locus_id}\t{chrom}\t{first}\t{last}\t{target_n}\t{variant_n}\n")

    with (out / "target_rad_regions.txt").open("w") as handle:
        for locus_id in sorted(locus_counts, key=lambda x: int(x[3:])):
            chrom, first, last, _, _ = locus_counts[locus_id]
            handle.write(f"{chrom}:{first}-{last}\n")

    with (out / "pair_locus_informative_snps.tsv").open("w") as handle:
        handle.write("pair\tlocus_id\tinformative_snps\n")
        for group1, group2 in PAIRS:
            pair = f"{group1}-{group2}"
            for locus_id in sorted(locus_counts, key=lambda x: int(x[3:])):
                handle.write(f"{pair}\t{locus_id}\t{pair_snp_counts[(pair, locus_id)]}\n")

    # ANGSD otherwise scans every scaffold in each BAM, including scaffolds
    # absent from the RAD target.  Restrict BAM traversal to target contigs.
    with (out / "target_contigs.txt").open("w") as handle:
        for chrom in intervals:
            if any(record[0] == chrom for record in locus_counts.values()):
                handle.write(f"{chrom}\n")

    with (out / "input_summary.tsv").open("w") as handle:
        handle.write("item\tvalue\n")
        for group in GROUPS:
            handle.write(f"{group}_individuals\t{sum(s.startswith(group) for s in samples)}\n")
        handle.write(f"target_sites\t{n_sites}\n")
        handle.write(f"variant_target_sites\t{n_variant}\n")
        handle.write(f"invariant_target_sites\t{n_invariant}\n")
        handle.write(f"operational_rad_loci\t{len(locus_counts)}\n")
        handle.write("rad_locus_definition\tCallable BED intervals joined when separated by <=500 bp\n")


if __name__ == "__main__":
    main()
