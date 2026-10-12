#!/usr/bin/env python3
"""Pool completed BA3-SNP chains after burn-in and calculate summaries."""
from pathlib import Path
import re
import numpy as np
import pandas as pd

root = Path('./orchid_paper2')
run = root/'18_contemporary_migration/runs/all_ld_6596'
out = root/'18_contemporary_migration/results/pooled_five_chains'
out.mkdir(parents=True, exist_ok=True)
burn = 1_000_000
frames = []
for c in range(1, 6):
    f = run/f'chain_{c}/orchid_18pop_all_ld_pruned.trace.txt'
    df = pd.read_csv(f, sep='\t')
    df = df[df.State > burn].copy()
    df['chain'] = c
    frames.append(df)
    print(f'chain_{c}: {len(df)} post-burn-in samples')
all_df = pd.concat(frames, ignore_index=True)
# The individual traces are retained as the source record; avoid writing a
# large duplicate pooled trace during the summary step.

params = [x for x in all_df.columns if re.fullmatch(r'm\[\d+\]\[\d+\]', x)]
rows=[]
for p in params:
    vals = all_df[p].to_numpy(float)
    q = np.quantile(vals, [.025,.5,.975])
    by = [x[p].to_numpy(float) for x in frames]
    means = np.array([x.mean() for x in by]); vars_ = np.array([x.var(ddof=1) for x in by])
    n = len(by[0]); m = len(by)
    W = vars_.mean(); B = n*means.var(ddof=1)
    varhat = (n-1)/n*W + B/n
    rhat = np.sqrt(varhat/W) if W > 0 else np.nan
    # Fast conservative lag-1 ESS estimate; use chain-specific diagnostics
    # or batch-means ESS for publication-grade estimates if needed.
    x = vals - vals.mean(); den = np.dot(x,x)
    rho1 = np.dot(x[:-1], x[1:]) / den if den > 0 else 0.0
    rho1 = min(max(rho1, -0.99), 0.99)
    ess = len(vals)*(1-rho1)/(1+rho1) if rho1 >= 0 else len(vals)
    i,j = map(int, re.findall(r'\d+',p))
    rows.append((p,i,j,vals.mean(),q[1],q[0],q[2],ess,rhat,(vals>0.01).mean()))
cols=['parameter','recipient_index','source_index','mean','median','q025','q975','ESS','Rhat','Pr_gt_0.01']
summary=pd.DataFrame(rows, columns=cols)
names=pd.read_csv(root/'18_contemporary_migration/inputs/sample_to_ba3_population.tsv', sep='\t')['ba3_population'].drop_duplicates().tolist()
if len(names)==18:
    summary['recipient']=summary.recipient_index.map(dict(enumerate(names)))
    summary['source']=summary.source_index.map(dict(enumerate(names)))
    summary.to_csv(out/'migration_posterior_summary.tsv', sep='\t', index=False)
else:
    summary.to_csv(out/'migration_posterior_summary.tsv', sep='\t', index=False)
mat=summary.pivot(index='recipient_index', columns='source_index', values='median')
mat.to_csv(out/'migration_median_matrix.tsv', sep='\t')
print(f'pooled samples: {len(all_df)}; parameters: {len(params)}')
print(summary.sort_values('median', ascending=False).head(20).to_string(index=False))
