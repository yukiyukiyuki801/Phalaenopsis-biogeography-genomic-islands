#!/usr/bin/env python3
"""Infer a rooted ML tree for every PAP--PSA RAD locus (no locus filtering)."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "orchid_paper2/19_local_genealogy"
IQTREE = Path("/home/ps2327/miniconda3/bin/iqtree")
MANIFEST = BASE / "inputs/iqtree_all650_diploid_manifest.tsv"
OUT = BASE / "trees_all650_diploid_ml"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "logs").mkdir(exist_ok=True)
(OUT / "cleaned_alignments").mkdir(exist_ok=True)
manifest = pd.read_csv(MANIFEST, sep="\t")


def clean_alignment(source, destination):
    records, name, seq = [], None, []
    for line in source.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq).upper()))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq).upper()))
    # IQ-TREE treats a sequence made entirely of N/IUPAC ambiguity as missing.
    kept = [(n, s) for n, s in records if n == "REF_P_equestris" or any(x in "ACGT" for x in s)]
    with destination.open("w") as handle:
        for n, s in kept:
            handle.write(f">{n}\n{s}\n")
    return kept


def run(row):
    locus = row.locus_id
    prefix = OUT / locus
    log = OUT / f"logs/{locus}.driver.log"
    treefile = Path(f"{prefix}.treefile")
    report = Path(f"{prefix}.iqtree")
    # A treefile is written before every IQ-TREE finishing step. Require the
    # final report as the completion sentinel so interrupted fits are retried.
    if treefile.exists() and report.exists():
        return {"locus_id": locus, "model_class": row.model_class,
                "exit_code": 0, "status": "existing_complete",
                "retained_sequences": None, "retained_PAP": None, "retained_PSA": None,
                "treefile_exists": True}
    cleaned = OUT / f"cleaned_alignments/{locus}.fa"
    kept = clean_alignment(ROOT / row.alignment, cleaned)
    kept_names = [x[0] for x in kept]
    command = [str(IQTREE), "-s", str(cleaned), "-pre", str(prefix),
               "-m", "HKY+F+G4", "-nt", "1", "-seed", "20261002",
               "-o", "REF_P_equestris", "-redo", "-quiet"]
    try:
        with log.open("w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                                       timeout=3600)
        code, status = completed.returncode, "complete" if completed.returncode == 0 else "failed"
    except subprocess.TimeoutExpired:
        code, status = 124, "timeout_1h"
    return {"locus_id": locus, "model_class": row.model_class,
            "exit_code": code, "status": status,
            "retained_sequences": len(kept_names),
            "retained_PAP": sum(x.startswith("PAP") for x in kept_names),
            "retained_PSA": sum(x.startswith("PSA") for x in kept_names),
            "treefile_exists": treefile.exists()}


rows = []
with ThreadPoolExecutor(max_workers=30) as pool:
    futures = [pool.submit(run, row) for row in manifest.itertuples(index=False)]
    for future in as_completed(futures):
        rows.append(future.result())
pd.DataFrame(rows).sort_values("locus_id").to_csv(
    BASE / "results/iqtree_all650_diploid_ml_run_summary.tsv", sep="\t", index=False)
