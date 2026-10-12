#!/usr/bin/env python3
"""Compute only the retained individual-level shared-SNP dosage distances."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
VCF=ROOT/'phalaenopsis_mac20260730/Phalaenopsis_eq_LD_filter.vcf'
OUT=ROOT/'orchid_paper2/results/neighbornet';OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for line in VCF.open():
    if line.startswith('##'):continue
    fields=line.rstrip().split('\t')
    if line.startswith('#CHROM'):samples=fields[9:];continue
    gi=fields[8].split(':').index('GT');row=[]
    for value in fields[9:]:
        gt=value.split(':')[gi].replace('|','/').split('/')
        row.append(sum(map(int,gt)) if len(gt)==2 and all(a in ('0','1') for a in gt) else np.nan)
    rows.append(row)
x=np.array(rows);n=len(samples);d=np.zeros((n,n))
assert x.shape==(6596,170)
for i in range(n):
    for j in range(i+1,n):
        called=np.isfinite(x[:,i])&np.isfinite(x[:,j])
        d[i,j]=d[j,i]=np.mean(np.abs(x[called,i]-x[called,j]))/2
pd.DataFrame(d,index=samples,columns=samples).rename_axis('sample').to_csv(OUT/'individual_prevosti_distance.tsv',sep='\t')
with (OUT/'splitstree_170_individuals_distances.nex').open('w') as out:
    out.write('#NEXUS\nBEGIN TAXA;\nDIMENSIONS NTAX=170;\nTAXLABELS\n')
    out.write('\n'.join("'"+s+"'" for s in samples)+'\n;\nEND;\nBEGIN DISTANCES;\nDIMENSIONS NTAX=170;\nFORMAT TRIANGLE=BOTH DIAGONAL LABELS;\nMATRIX\n')
    for s,row in zip(samples,d):out.write("'"+s+"' "+' '.join(f'{v:.10g}' for v in row)+'\n')
    out.write(';\nEND;\n')
print('Wrote 170-individual distance matrix and SplitsTree NEXUS.')
