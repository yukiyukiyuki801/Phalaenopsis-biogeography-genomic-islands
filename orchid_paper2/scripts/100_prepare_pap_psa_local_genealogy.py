#!/usr/bin/env python3
"""Prepare PAP--PSA RAD-locus and sample inputs for local genealogy analysis."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "orchid_paper2/19_local_genealogy"
INPUT = ROOT / "orchid_paper2/results/rad_locus_statistics/rad_locus_statistics.four_models_central25_dxy.tsv"

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "inputs").mkdir(exist_ok=True)
(OUT / "phased_samples").mkdir(exist_ok=True)
(OUT / "logs").mkdir(exist_ok=True)
(OUT / "alignments").mkdir(exist_ok=True)
(OUT / "trees").mkdir(exist_ok=True)
(OUT / "results").mkdir(exist_ok=True)
(OUT / "plots").mkdir(exist_ok=True)

d = pd.read_csv(INPUT, sep="\t", na_values="NA")
d = d[(d["pair"] == "PAP-PSA") & (d["chrom"] != "NC_017609.1")].copy()
d = d.sort_values(["chrom", "start", "end"])

keep = ["locus_id", "chrom", "start", "end", "callable_bases",
        "informative_snps", "pi_group1", "pi_group2", "dxy", "hudson_fst",
        "tajima_d_group1", "tajima_d_group2",
        "four_model_class_central25_dxy"]
d[keep].to_csv(OUT / "inputs/pap_psa_nuclear_loci.tsv", sep="\t", index=False)

# BED is zero-based, half-open. The source start/end are one-based inclusive.
with open(OUT / "inputs/pap_psa_nuclear_loci.bed", "w") as handle:
    for row in d.itertuples():
        handle.write(f"{row.chrom}\t{row.start-1}\t{row.end}\t{row.locus_id}\n")

samples = [f"PAP{i:03d}" for i in range(1, 79)] + [f"PSA{i:03d}" for i in range(1, 8)]
(OUT / "inputs/pap_psa_samples.txt").write_text("\n".join(samples) + "\n")
with open(OUT / "inputs/pap_psa_sample_bams.tsv", "w") as handle:
    handle.write("sample\tspecies\tbam\n")
    for sample in samples:
        handle.write(f"{sample}\t{sample[:3]}\torchid_paper2/04_bam/{sample}.sorted.bam\n")

counts = d["four_model_class_central25_dxy"].value_counts()
with open(OUT / "inputs/model_class_counts.tsv", "w") as handle:
    handle.write("model_class\tn_nuclear_loci\n")
    for model, count in counts.items():
        handle.write(f"{model}\t{count}\n")

