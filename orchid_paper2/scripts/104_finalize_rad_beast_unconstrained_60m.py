#!/usr/bin/env python3
"""Finalize completed 60M lineage-unconstrained RAD-BEAST chains."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'orchid_paper2/16_rad_beast'
OUT=BASE/'final_unconstrained_60m_summary'; OUT.mkdir(parents=True,exist_ok=True)
BURN=6_000_000

def split_rhat(a):
 n=min(map(len,a)); a=np.stack([x[-n:] for x in a]); h=n//2; z0=np.concatenate([a[:,:h],a[:,-h:]],axis=0)
 def basic(z):
  nn=z.shape[1]; W=z.var(1,ddof=1).mean(); B=nn*z.mean(1).var(ddof=1); return np.sqrt(((nn-1)/nn*W+B/nn)/W) if W else np.nan
 flat=z0.ravel(); z=norm.ppf((rankdata(flat)-.375)/(len(flat)+.25)).reshape(z0.shape)
 folded=np.abs(z0-np.median(flat)).ravel(); zf=norm.ppf((rankdata(folded)-.375)/(len(folded)+.25)).reshape(z0.shape)
 return max(basic(z),basic(zf))

def ess_bulk(a):
 n=min(map(len,a)); a=np.stack([x[-n:] for x in a]); m=len(a); flat=a.ravel(); z=norm.ppf((rankdata(flat)-.375)/(len(flat)+.25)).reshape(a.shape)
 W=z.var(1,ddof=1).mean(); B=n*z.mean(1).var(ddof=1); varp=(n-1)/n*W+B/n; ac=[]
 for row in z:
  y=row-row.mean(); f=np.fft.rfft(y,n=2*n); v=np.fft.irfft(f*f.conjugate(),n=2*n)[:n]/np.arange(n,0,-1); ac.append(v)
 rho=1-(W-np.asarray(ac).mean(0))/varp; rho[0]=1; total=0
 for t in range(1,n-1,2):
  pair=rho[t]+rho[t+1]
  if pair<0: break
  total+=pair
 return min(m*n,m*n/max(1,1+2*total))

def single_ess(x):
 x=np.array(x,dtype=float,copy=True); n=len(x); x-=x.mean(); f=np.fft.rfft(x,n=2*n); a=np.fft.irfft(f*f.conjugate(),n=2*n)[:n]/np.arange(n,0,-1); a/=a[0]; total=0
 for t in range(1,n-1,2):
  pair=a[t]+a[t+1]
  if pair<0: break
  total+=pair
 return min(n,n/max(1,1+2*total))

def run_dir(s,c): return BASE/f'runs/{s}_combined_rateprior_unconstrained_60m/chain{c}'
def read_log(s,c):
 f=next(run_dir(s,c).glob('*.log')); x=pd.read_csv(f,sep='\t',comment='#').sort_values('Sample'); return x[x.Sample>=BURN].drop_duplicates('Sample')
def header(f):
 z=[]
 for l in open(f):
  if l.startswith('tree STATE_'): break
  z.append(l)
 return z
def merge_trees(s):
 first=next(run_dir(s,1).glob('*.trees')); out=OUT/f'{s}_unconstrained_merged_postburn.trees'; n=0
 with open(out,'w') as h:
  h.writelines(header(first))
  for c in range(1,4):
   f=next(run_dir(s,c).glob('*.trees'))
   for l in open(f):
    if not l.startswith('tree STATE_'): continue
    st=int(re.search(r'STATE_(\d+)',l).group(1))
    if st<BURN: continue
    h.write(re.sub(r'tree STATE_\d+',f'tree STATE_{n}',l,1)); n+=1
  h.write('End;\n')
 return n

for s in ['zhang','givnish']:
 chains=[read_log(s,c) for c in range(1,4)]
 params=[p for p in chains[0] if p!='Sample' and all(pd.api.types.is_numeric_dtype(x[p]) for x in chains)]
 rows=[]
 for p in params:
  vals=[x[p].to_numpy(float) for x in chains]; pooled=np.concatenate(vals)
  rows.append({'parameter':p,'median':np.median(pooled),'q025':np.quantile(pooled,.025),'q975':np.quantile(pooled,.975),'ess_bulk':ess_bulk(vals),'rhat_rank_split':split_rhat(vals),'min_chain_ess':min(single_ess(v) for v in vals)})
 pd.DataFrame(rows).to_csv(OUT/f'{s}_unconstrained_scalar_diagnostics.tsv',sep='\t',index=False)
 for c,x in enumerate(chains,1): x.assign(chain=c).to_csv(OUT/f'{s}_unconstrained_chain{c}_postburn.log',sep='\t',index=False)
 (OUT/f'{s}_unconstrained_tree_count.txt').write_text(str(merge_trees(s))+'\n')
