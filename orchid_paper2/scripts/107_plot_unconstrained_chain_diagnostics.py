#!/usr/bin/env python3
"""Trace and marginal-density panels for final unconstrained chains."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'orchid_paper2/16_rad_beast/final_unconstrained_60m_summary';OUT=BASE/'figures';OUT.mkdir(exist_ok=True)
PARAM=[('posterior','Posterior'),('clockRate','Clock rate'),('ucldStdev','UCLN SD'),('treeHeight','PEQ–ingroup root age'),('mrca.age(ingroupTaxa)','29-tip ingroup crown age')];COL=['#1b9e77','#d95f02','#7570b3']
for s in ['zhang','givnish']:
 frames=[pd.read_csv(BASE/f'{s}_unconstrained_chain{c}_postburn.log',sep='\t') for c in range(1,4)]
 fig,axs=plt.subplots(len(PARAM),2,figsize=(13,12))
 for r,(p,label) in enumerate(PARAM):
  for c,x in enumerate(frames):
   axs[r,0].plot(x.Sample/1e6,x[p],lw=.45,alpha=.72,color=COL[c],label=f'Chain {c+1}')
   axs[r,1].hist(x[p],bins=55,density=True,histtype='step',lw=1.2,color=COL[c],label=f'Chain {c+1}')
  axs[r,0].set_ylabel(label);axs[r,0].grid(alpha=.2);axs[r,1].grid(alpha=.2)
 axs[-1,0].set_xlabel('State (millions; after 6M burn-in)');axs[-1,1].set_xlabel('Posterior value');axs[0,0].legend(ncol=3,fontsize=8);axs[0,1].legend(ncol=3,fontsize=8);fig.suptitle(f'Lineage-unconstrained RAD-BEAST diagnostics — {s.capitalize()} calibration',fontsize=14);fig.tight_layout(rect=(0,0,1,.98));fig.savefig(OUT/f'{s}_unconstrained_trace_density.pdf');fig.savefig(OUT/f'{s}_unconstrained_trace_density.png',dpi=300);plt.close(fig)
