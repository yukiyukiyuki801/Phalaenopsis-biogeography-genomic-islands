#!/usr/bin/env python3
import pandas as pd, numpy as np, json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, magma

ROOT=Path('.'); base=ROOT/'orchid_paper2/18_contemporary_migration'; res=base/'results/pooled_five_chains'; out=res/'figures'; out.mkdir(parents=True,exist_ok=True)
s=pd.read_csv(res/'migration_posterior_summary.tsv',sep='\t'); s=s[s.recipient!=s.source].copy()
meta=pd.read_csv(ROOT/'phalaenopsis_mac20260730/full_170samples_popmap.csv'); cls=pd.read_csv(base/'final_18_population_classification.tsv',sep='\t')
meta['species']=meta.id.str.extract(r'^(PAM|PAP|PSA)')
loc=meta.groupby('pop',as_index=False).agg(longitude=('longitude','first'),latitude=('latitude','first'),species=('species','first'))
loc=loc.merge(cls[['original_population','primary_ba3_population','action']],left_on='pop',right_on='original_population')
loc=loc[loc.action!='remove']; pop=loc.groupby('primary_ba3_population',as_index=False).agg(longitude=('longitude','mean'),latitude=('latitude','mean'),species=('species',lambda x:x.mode().iloc[0]))
pop=pop.set_index('primary_ba3_population')
# heatmap
p=s.pivot(index='recipient',columns='source',values='median'); p=p.reindex(index=sorted(p.index),columns=sorted(p.columns))
fig,ax=plt.subplots(figsize=(10,9)); im=ax.imshow(p.values,cmap='magma',aspect='auto'); fig.colorbar(im,ax=ax,label='Posterior median'); ax.set_xticks(range(len(p.columns)),p.columns,rotation=45,ha='right'); ax.set_yticks(range(len(p.index)),p.index); ax.set_title('BA3-SNP migration posterior\nRows are recipients; columns are sources (cyan: Rhat > 1.05)')
species_tick={'PAM':'#D95F02','PAP':'#1B9E77','PSA':'#7570B3'}
for lab in ax.get_xticklabels():
 lab.set_color(species_tick.get(pop.loc[lab.get_text(),'species'],'black'))
for lab in ax.get_yticklabels():
 lab.set_color(species_tick.get(pop.loc[lab.get_text(),'species'],'black'))
# Colored species strips make the species assignment visible even when labels are dense.
for i,name in enumerate(p.columns): ax.add_patch(plt.Rectangle((i-.5,-1.05),1,.12,color=species_tick.get(pop.loc[name,'species'],'black'),clip_on=False))
for i,name in enumerate(p.index): ax.add_patch(plt.Rectangle((-1.05,i-.5),.12,1,color=species_tick.get(pop.loc[name,'species'],'black'),clip_on=False))
for _,r in s[s.Rhat>1.05].iterrows():
 if r.source in p.columns and r.recipient in p.index: ax.plot(list(p.columns).index(r.source),list(p.index).index(r.recipient),'x',color='cyan',ms=7,mew=1.3)
fig.tight_layout(); fig.savefig(out/'ba3_pooled_migration_heatmap.png',dpi=400); fig.savefig(out/'ba3_pooled_migration_heatmap.pdf'); plt.close(fig)
# curved arrows routed away from centroid; top edges only to keep map readable
e=s[(s['median']>=.01)&(s['Pr_gt_0.01']>=.50)].sort_values('median',ascending=False).head(25).copy(); e['rank']=np.arange(len(e)); cx=pop.longitude.mean(); cy=pop.latitude.mean()
cols={'PAM':'#D95F02','PAP':'#1B9E77','PSA':'#7570B3'}; fig,ax=plt.subplots(figsize=(13,8.5)); ax.set_facecolor('#f7f7f7'); norm=Normalize(e['median'].min(),e['median'].max())
boundary_files=[base/'inputs/ne_10m_admin_0_countries.geojson']
for gf in boundary_files:
 if not gf.is_file():
  raise FileNotFoundError('Download the complete Natural Earth 1:10m country boundaries as documented; do not plot an incomplete background.')
 if gf.is_file():
  gj=json.loads(gf.read_text())
  for feat in gj.get('features',[]):
   geom=feat.get('geometry',{}); coords=geom.get('coordinates',[]); rings=coords if geom.get('type')=='MultiPolygon' else [coords]
   for poly in rings:
    for ring in poly:
     a=np.asarray(ring)
     if len(a)>1:
      ax.fill(a[:,0],a[:,1],color='#e8e8e8',zorder=0)
      ax.plot(a[:,0],a[:,1],color='#777',lw=.45,zorder=1)
for _,r in e.iterrows():
 if r.source not in pop.index or r.recipient not in pop.index: continue
 a,b=pop.loc[r.source],pop.loc[r.recipient]; dx=b.longitude-a.longitude; dy=b.latitude-a.latitude; d=np.hypot(dx,dy) or 1; px,py=-dy/d,dx/d; mx,my=(a.longitude+b.longitude)/2,(a.latitude+b.latitude)/2; sign=1 if (mx-cx)*px+(my-cy)*py>=0 else -1; rad=sign*(0.12+0.018*r['rank'])*(0.7+0.3*min(1,d/8)); con=cols.get(pop.loc[r.recipient,'species'],'#555'); patch=FancyArrowPatch((a.longitude,a.latitude),(b.longitude,b.latitude),connectionstyle=f'arc3,rad={rad:.3f}',arrowstyle='-|>',mutation_scale=10,linewidth=.5+3*(r['median']-e['median'].min())/(e['median'].max()-e['median'].min()+1e-9),color=con,alpha=.72,linestyle='--' if r.Rhat>1.05 else '-') ; ax.add_patch(patch)
for name,r in pop.iterrows(): ax.scatter(r.longitude,r.latitude,s=65,color=cols.get(r.species,'grey'),edgecolor='black',zorder=4,marker={'PAM':'o','PAP':'^','PSA':'s'}.get(r.species,'o')); ax.annotate(name,(r.longitude,r.latitude),xytext=(4,4),textcoords='offset points',fontsize=8,weight='bold')
ax.set_xlim(95,150); ax.set_ylim(-20,24); ax.set_aspect('equal'); ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude'); ax.set_title('Pooled BA3-SNP migration among 18 populations\nTop 25 directional edges ranked by migration strength (posterior median); arrows source → recipient'); ax.grid(color='white',lw=.7)
species_handles=[Line2D([0],[0],marker={'PAM':'o','PAP':'^','PSA':'s'}[sp],color='w',markerfacecolor=co,markeredgecolor='black',markersize=8,label=sp) for sp,co in cols.items()]
lo=float(e['median'].min()); hi=float(e['median'].max()); mid=(lo+hi)/2
width_handles=[Line2D([0],[0],color='#555',lw=0.8,label=f'{lo:.3f}'),Line2D([0],[0],color='#555',lw=1.8,label=f'{mid:.3f}'),Line2D([0],[0],color='#555',lw=3,label=f'{hi:.3f}'),Line2D([0],[0],color='#555',lw=1.5,ls='--',label='Rhat > 1.05')]
leg1=ax.legend(handles=species_handles,title='Population species',loc='upper left',frameon=True); ax.add_artist(leg1); ax.legend(handles=width_handles,title='Arrow encoding',loc='lower left',frameon=True)
fig.tight_layout(); fig.savefig(out/'ba3_pooled_migration_curved_map.png',dpi=400); fig.savefig(out/'ba3_pooled_migration_curved_map.pdf'); plt.close(fig)
# intervals
t=s.sort_values('median',ascending=False).head(20).iloc[::-1].copy(); fig,ax=plt.subplots(figsize=(10,7)); y=np.arange(len(t)); ax.hlines(y,t.q025,t.q975,color=['red' if x>1.05 else 'black' for x in t.Rhat]); ax.scatter(t['median'],y,c=['red' if x>1.05 else 'black' for x in t.Rhat]); ax.set_yticks(y,[f'{a} → {b}' for a,b in zip(t.source,t.recipient)]); ax.set_xlabel('Migration posterior (median and 95% interval)'); ax.set_title('Strongest pooled BA3-SNP migration edges'); fig.tight_layout(); fig.savefig(out/'ba3_pooled_migration_top_edges.png',dpi=400); fig.savefig(out/'ba3_pooled_migration_top_edges.pdf'); plt.close(fig)
e.to_csv(out/'ba3_curved_map_edges.tsv',sep='\t',index=False)
e[['rank','source','recipient','median','q025','q975','Pr_gt_0.01','ESS','Rhat']].rename(columns={'median':'migration_median','q025':'migration_q025','q975':'migration_q975','Pr_gt_0.01':'Pr_m_gt_0.01'}).to_csv(out/'ba3_final_ranked_migration_strength.tsv',sep='\t',index=False)
