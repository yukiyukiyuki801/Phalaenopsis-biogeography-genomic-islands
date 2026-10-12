#!/usr/bin/env python3
"""Nearest-mRNA annotation within 10 kb for PAP–PSA balancing-like RAD loci."""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
OUT = ROOT / '23_balancing_function'
for folder in ['inputs', 'tables']:
    (OUT / folder).mkdir(parents=True, exist_ok=True)
CLASS = ROOT / 'results/rad_locus_statistics/rad_locus_statistics.four_models_central25_dxy.tsv'
GFF = PROJECT / 'GEA/orchid.gff'
MAPPING = PROJECT / 'GEA/all_info_new_gene.txt'
FAI = ROOT / '01_reference/GCF_001263595.1_ASM126359v1_genomic.fna.fai'
DISTANCE = 10000



def write_table(path, rows, fields=None):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields or list(rows[0]), extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def checksum(path):
    sha = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            sha.update(block)
    return sha.hexdigest()


with CLASS.open() as handle:
    loci = [r for r in csv.DictReader(handle, delimiter='\t')
            if r['pair'] == 'PAP-PSA' and r['chrom'] != 'NC_017609.1']
assert len(loci) == 650
candidates = {r['locus_id'] for r in loci if r['four_model_class_central25_dxy'] == 'balancing_selection'}
assert len(candidates) == 14
write_table(OUT / 'inputs/old_pap_psa_all650_classification.tsv', loci)
write_table(OUT / 'inputs/old_pap_psa_balancing14.tsv', [r for r in loci if r['locus_id'] in candidates])
with FAI.open() as handle:
    lengths = {p[0]: int(p[1]) for line in handle if (p := line.split())}
assert all(r['chrom'] in lengths and 1 <= int(r['start']) <= int(r['end']) <= lengths[r['chrom']] for r in loci)

mapping = defaultdict(set)
symbols = defaultdict(set)
with MAPPING.open() as handle:
    for line in handle:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2:
            transcript = p[0].split('.')[0]
            for at in re.findall(r'AT[1-5MC]G\d{5}', p[1], flags=re.I):
                at = at.upper()
                mapping[transcript].add(at)
                if len(p) > 2:
                    symbols[at].add(p[2])
transcripts = defaultdict(list)
gff_regions = {}
assembly = ''
with GFF.open() as handle:
    for line in handle:
        if line.startswith('#!genome-build-accession'):
            assembly = line.strip()
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9:
            continue
        if p[2] == 'region':
            gff_regions[p[0]] = int(p[4])
        if p[2] != 'mRNA':
            continue
        attrs = dict(v.split('=', 1) for v in p[8].split(';') if '=' in v)
        transcript = attrs.get('transcript_id', attrs.get('Name', '')).split('.')[0]
        transcripts[p[0]].append(dict(transcript_id=transcript, feature_id=attrs.get('ID', ''),
                                      gene_id=attrs.get('Parent', ''), product=unquote(attrs.get('product', '')),
                                      gene_start=int(p[3]), gene_end=int(p[4]), strand=p[6]))
assert 'GCF_001263595.1' in assembly
assert all(gff_regions.get(r['chrom']) == lengths[r['chrom']] for r in loci), 'GFF/reference scaffold length mismatch'
print(f'Validated reference assembly and lengths of all {len({r["chrom"] for r in loci})} assayed scaffolds', flush=True)

annotations, locus_summary = [], []
for locus in loci:
    chrom, start, end = locus['chrom'], int(locus['start']), int(locus['end'])
    available = transcripts.get(chrom, [])
    distances = [max(t['gene_start']-end, start-t['gene_end'], 0) for t in available]
    nearest = min(distances) if distances else None
    eligible = nearest is not None and nearest <= DISTANCE
    selected = [t for t,d in zip(available,distances) if d == nearest] if eligible else []
    at_ids = set()
    for t in selected:
        ats = sorted(mapping.get(t['transcript_id'], set()))
        at_ids.update(ats)
        annotations.append(dict(locus_id=locus['locus_id'], model_class=locus['four_model_class_central25_dxy'],
                                chrom=chrom, locus_start=start, locus_end=end, distance_bp=nearest,
                                **t, arabidopsis_ids=';'.join(ats),
                                arabidopsis_symbols=';'.join(sorted({s for at in ats for s in symbols[at]}))))
    locus_summary.append(dict(locus_id=locus['locus_id'], model_class=locus['four_model_class_central25_dxy'],
                              chrom=chrom, start=start, end=end,
                              nearest_distance_bp=nearest if nearest is not None else '',
                              assigned_within_10kb=eligible, n_nearest_mrnas=len(selected),
                              n_orchid_genes=len({t['gene_id'] for t in selected}),
                              n_arabidopsis_ids=len(at_ids), arabidopsis_ids=';'.join(sorted(at_ids)),
                              status='mapped' if at_ids else 'nearest_mrna_no_AT_mapping' if selected else 'no_mrna_within_10kb'))
write_table(OUT / 'tables/all650_nearest_mrna_annotations.tsv', annotations)
write_table(OUT / 'tables/all650_locus_annotation_status.tsv', locus_summary)
candidate_annotations = [r for r in annotations if r['locus_id'] in candidates]
write_table(OUT / 'tables/balancing14_candidate_genes.tsv', candidate_annotations,
            fields=list(annotations[0]))
gene_groups = defaultdict(list)
for row in candidate_annotations:
    gene_groups[(row['locus_id'],row['gene_id'])].append(row)
gene_summary=[]
for (locus_id,gene_id),rows in sorted(gene_groups.items()):
    gene_summary.append(dict(locus_id=locus_id, orchid_gene=gene_id,
                             distance_bp=min(r['distance_bp'] for r in rows),
                             n_mrna_isoforms=len({r['transcript_id'] for r in rows}),
                             reference_products='; '.join(sorted({r['product'] for r in rows})),
                             arabidopsis_ids=';'.join(sorted({v for r in rows for v in r['arabidopsis_ids'].split(';') if v})),
                             arabidopsis_symbols=';'.join(sorted({v for r in rows for v in r['arabidopsis_symbols'].split(';') if v}))))
write_table(OUT / 'tables/balancing_candidate_gene_function_summary.tsv',gene_summary)
write_table(OUT / 'tables/balancing14_locus_annotation_status.tsv', [r for r in locus_summary if r['locus_id'] in candidates])
query = sorted({at for r in candidate_annotations for at in r['arabidopsis_ids'].split(';') if at})
background = sorted({at for r in annotations for at in r['arabidopsis_ids'].split(';') if at})
assert set(query) <= set(background)
for name, ids in [('candidate_arabidopsis_ids.txt',query),('all650_RAD_background_arabidopsis_ids.txt',background)]:
    (OUT / 'inputs' / name).write_text('\n'.join(ids)+'\n')
counts = dict(n_nuclear_loci=len(loci), n_balancing_loci=len(candidates),
              candidate_annotation_status=dict(Counter(r['status'] for r in locus_summary if r['locus_id'] in candidates)),
              n_candidate_orchid_genes=len({r['gene_id'] for r in candidate_annotations}),
              n_candidate_mrna_ids=len({r['transcript_id'] for r in candidate_annotations}),
              n_candidate_arabidopsis_ids=len(query), n_background_arabidopsis_ids=len(background),
              distance_cutoff_bp=DISTANCE, assembly=assembly)
(OUT / 'tables/annotation_summary.json').write_text(json.dumps(counts,indent=2)+'\n')
manifest = [dict(path=str(p.relative_to(PROJECT)), sha256=checksum(p)) for p in [CLASS,GFF,MAPPING,FAI]]
write_table(OUT / 'inputs/input_checksums.tsv',manifest)
print(json.dumps(counts,indent=2), flush=True)

