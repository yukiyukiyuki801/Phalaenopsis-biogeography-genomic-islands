#!/usr/bin/env python3
"""Summarize record types and sample missingness for filtered RAD VCFs."""

from pathlib import Path
import csv, gzip, sys

outdir = Path(sys.argv[1])
rows = []
for profile in ("paper1_like", "primary", "strict"):
    path = outdir / f"orchid_173.{profile}.biallelic_snps_invariant.vcf.gz"
    counts = dict(records=0, invariant=0, variant=0)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            counts["records"] += 1
            if line.split("\t", 5)[4] == ".":
                counts["invariant"] += 1
            else:
                counts["variant"] += 1
    global_path = outdir / f"orchid_173.{profile}.global_call30.vcf.gz"
    global_n = 0
    with gzip.open(global_path, "rt") as fh:
        global_n = sum(1 for line in fh if not line.startswith("#"))
    rows.append(dict(profile=profile, global_call30_records=global_n, **counts))

with (outdir / "filter_record_summary.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
    w.writeheader(); w.writerows(rows)
