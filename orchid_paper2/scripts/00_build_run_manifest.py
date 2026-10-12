#!/usr/bin/env python3
"""Build Paper 2 run manifests from the NCBI SRA RunInfo table."""

import csv
from pathlib import Path

root = Path("orchid_paper2/00_metadata")
runinfo = root / "sra_runinfo.csv"
manifest = root / "run_manifest.tsv"
old_samples = root / "paper1_final_vcf.samples.txt"
released_not_final = root / "released_not_in_paper1_final_vcf.samples.txt"
vcf = Path("phalaenopsis_mac20260730/Phalaenopsis_eq_final_snps.recode.vcf")

with runinfo.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

if not rows:
    raise SystemExit(f"No SRA rows in {runinfo}")

with vcf.open() as handle:
    header = next(line for line in handle if line.startswith("#CHROM"))
vcf_samples = header.rstrip("\n").split("\t")[9:]

seen_runs = set()
seen_samples = set()
for row in rows:
    run = row["Run"]
    sample = row["SampleName"]
    if row["LibraryLayout"] != "PAIRED":
        raise SystemExit(f"Non-paired run: {run}")
    if run in seen_runs or sample in seen_samples:
        raise SystemExit(f"Duplicate run or sample: {run} {sample}")
    seen_runs.add(run)
    seen_samples.add(sample)

missing = sorted(set(vcf_samples) - seen_samples)
if missing:
    raise SystemExit(f"Final-VCF samples absent from SRA: {missing}")

with manifest.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["run_accession", "sample_id", "library_id", "biosample"])
    for row in sorted(rows, key=lambda x: x["SampleName"]):
        writer.writerow(
            [row["Run"], row["SampleName"], row["LibraryName"], row["BioSample"]]
        )

old_samples.write_text("\n".join(vcf_samples) + "\n")
extra = sorted(seen_samples - set(vcf_samples))
released_not_final.write_text("\n".join(extra) + "\n")

print(f"manifest_runs={len(rows)} final_vcf_samples={len(vcf_samples)} extras={len(extra)}")
