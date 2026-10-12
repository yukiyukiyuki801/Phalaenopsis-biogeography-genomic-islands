#!/usr/bin/env python3
"""Prepare one representative protein per gene from NCBI protein FASTA files.

NCBI protein headers commonly include a bracketed gene identifier. When no gene
identifier is present, the protein accession is treated as its own gene. The
longest non-redundant amino-acid sequence is retained for each inferred gene.
"""

from __future__ import annotations

import argparse
import gzip
import re
from urllib.parse import unquote
from pathlib import Path


GENE_RE = re.compile(r"\[gene=([^\]]+)\]")


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def fasta_records(path: Path):
    header = None
    sequence = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
            yield header, "".join(sequence)


def attributes(field: str):
    parsed = {}
    for item in field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = unquote(value)
    return parsed


def protein_gene_map(gff_path: Path):
    transcript_to_gene = {}
    protein_to_transcript = {}
    with open_text(gff_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2]
            attrs = attributes(fields[8])
            if feature in {"mRNA", "transcript"}:
                transcript = attrs.get("ID")
                parent = attrs.get("Parent")
                if transcript and parent:
                    transcript_to_gene[transcript] = parent.split(",")[0]
            elif feature == "CDS":
                protein = attrs.get("protein_id")
                parent = attrs.get("Parent")
                if protein and parent:
                    protein_to_transcript[protein] = parent.split(",")[0]
    mapping = {}
    for protein, transcript in protein_to_transcript.items():
        mapping[protein] = transcript_to_gene.get(transcript, transcript)
    return mapping


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--gff", type=Path)
    parser.add_argument("--taxon", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    accession_to_gene = protein_gene_map(args.gff) if args.gff else {}
    best = {}
    for header, sequence in fasta_records(args.input):
        accession = header.split()[0]
        match = GENE_RE.search(header)
        gene = accession_to_gene.get(accession)
        if gene is None:
            gene = match.group(1) if match else accession
        clean = sequence.rstrip("*").replace("*", "X")
        if len(clean) < 30:
            continue
        if gene not in best or len(clean) > len(best[gene][1]):
            best[gene] = (accession, clean)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for gene in sorted(best):
            accession, sequence = best[gene]
            handle.write(f">{args.taxon}|{gene}|{accession}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    print(f"{args.taxon}\t{len(best)}\t{args.output}")


if __name__ == "__main__":
    main()
