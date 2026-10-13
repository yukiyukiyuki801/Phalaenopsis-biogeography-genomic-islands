#!/usr/bin/env python3
"""Final figures/tables for lineage-unconstrained 60M RAD-BEAST analyses."""
from pathlib import Path
from collections import defaultdict
import io,re,gc
import numpy as np, pandas as pd
from Bio import Phylo
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]; BASE=ROOT/'orchid_paper2/16_rad_beast'; IN=BASE/'final_unconstrained_60m_summary'; FIG=IN/'figures'; FIG.mkdir(exist_ok=True)
PAM=set('AMB PAL PNG PON QL SAB SK SUL SUR TUA WJ'.split()); PSA={'MDA'}; COL={'PAM':'#D95F02','PAP':'#1B9E77','PSA':'#7570B3','PEQ':'#444'}; BURN=6_000_000
def lineage(x): return 'PEQ' if x=='PEQ' else ('PAM' if x in PAM else ('PSA' if x in PSA else 'PAP'))
def trans(path):
 d={};on=False
 for l in open(path):
  if l.strip().lower()=='translate':on=True;continue
  if on:
   z=re.match(r'\s*(\d+)\s+([^,;]+)',l)
   if z:d[z.group(1)]=z.group(2).strip()
   if ';' in l:break
 return d
def parse(l,tr,comments=False):
 s=l.split('=',1)[1].strip();s=re.sub(r'^\[&R\]\s*','',s)
 if not comments:s=re.sub(r'\[&[^\]]*\]','',s)
 t=Phylo.read(io.StringIO(s),'newick')
 for q in t.get_terminals():q.name=tr.get(q.name,q.name)
 return t
def mcc(path):
 tr=trans(path);return parse(next(l for l in open(path) if l.startswith('tree ')),tr)
def support_map(path):
 tr=trans(path);t=parse(next(l for l in open(path) if l.startswith('tree ')),tr,True);d={}
 for c in t.get_nonterminals():
  z=re.search(r'(?:^|,)posterior=([0-9.eE+-]+)',c.comment or '');d[frozenset(q.name for q in c.get_terminals())]=float(z.group(1)) if z else 1
 return d
def layout(t,order=None):
 order=order or [q.name for q in t.get_terminals()];y={n:len(order)-i for i,n in enumerate(order)}
 def yf(c): return y[c.name] if c.is_terminal() else sum(yf(x) for x in c.clades)/len(c.clades)
 dep=t.depths();H=max(dep[q] for q in t.get_terminals());age={c:H-dep[c] for c in dep};return y,yf,age,H
def draw(ax,t,title,order=None,support=None):
 y,yf,age,H=layout(t,order)
 def rec(c):
  if not c.clades:return
  ys=[yf(x) for x in c.clades];ax.plot([age[c]]*2,[min(ys),max(ys)],c='#222',lw=1)
  for x in c.clades:ax.plot([age[c],age[x]],[yf(x)]*2,c='#222',lw=1);rec(x)
 rec(t.root)
 for n,v in y.items():ax.text(-H*.02,v,n,c=COL[lineage(n)],va='center',ha='left',fontsize=8,fontweight='bold')
 if support:
  for c in t.get_nonterminals():
   k=frozenset(q.name for q in c.get_terminals())
   if len(k)<len(y):ax.text(age[c]+H*.007,yf(c)-.16,f'{support.get(k,np.nan):.2f}',fontsize=5.8,c='#7A1F5C')
 ax.set_xlim(H*1.07,-H*.12);ax.set_ylim(0,len(y)+1);ax.set_yticks([]);ax.set_xlabel('Time before present (Ma)');ax.grid(axis='x',c='#ddd',lw=.5);ax.set_title(title)
 return y,yf,age,H
def save(fig,name):fig.tight_layout();fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');fig.savefig(FIG/f'{name}.png',dpi=350,bbox_inches='tight');plt.close(fig)
def sampled_and_frequencies(s):
 merged=IN/f'{s}_unconstrained_merged_postburn.trees';tr=trans(merged);nt=sum(l.startswith('tree STATE_') for l in open(merged));pick=set(np.linspace(0,nt-1,100,dtype=int));trees=[];counts=defaultdict(int);bychain=[defaultdict(int) for _ in range(3)];nper=nt//3;stats=np.zeros((3,4),int);pam=pap=paponly=both=0;i=-1
 for l in open(merged):
  if not l.startswith('tree STATE_'):continue
  i+=1
  t=parse(l,tr);names={q.name for q in t.get_terminals()};allnames=names; papset=names-PAM-{'PEQ'}
  pap_only=papset-PSA; pm={q.name for q in t.common_ancestor(PAM).get_terminals()}==PAM; pp={q.name for q in t.common_ancestor(papset).get_terminals()}==papset; po={q.name for q in t.common_ancestor(pap_only).get_terminals()}==pap_only
  pam+=pm;pap+=pp;paponly+=po;both+=pm and pp;ci=min(i//nper,2);stats[ci]+=np.array([pm,po,pp,pm and pp])
  if i in pick:trees.append(t)
  if i%10==0:
   for c in t.get_nonterminals():
    k=frozenset(q.name for q in c.get_terminals());counts[k]+=1;bychain[ci][k]+=1
  if i%500==0:gc.collect()
 n=nt; overall={'scenario':s,'chain':'pooled','trees':n,'pam_monophyly_frequency':pam/n,'pap_only_monophyly_frequency':paponly/n,'pap_psa_monophyly_frequency':pap/n,'reciprocal_monophyly_frequency':both/n}; chainrows=[]
 for c in range(3):chainrows.append({'scenario':s,'chain':c+1,'trees':nper,'pam_monophyly_frequency':stats[c,0]/nper,'pap_only_monophyly_frequency':stats[c,1]/nper,'pap_psa_monophyly_frequency':stats[c,2]/nper,'reciprocal_monophyly_frequency':stats[c,3]/nper})
 topo=[];maj=[{k for k,v in d.items() if v/max(d.values())>=.5} for d in bychain]
 for a in range(3):
  for b in range(a+1,3):topo.append({'scenario':s,'chain_a':a+1,'chain_b':b+1,'majority_clades_a':len(maj[a]),'majority_clades_b':len(maj[b]),'majority_clade_jaccard':len(maj[a]&maj[b])/len(maj[a]|maj[b])})
 return trees,counts,len(range(0,n,10)),overall,chainrows,topo
def tangle(left,right,name,title):
 lo=[q.name for q in left.get_terminals()];ro=[q.name for q in right.get_terminals()];yl,fl,al,hl=layout(left,lo);yr,fr,ar,hr=layout(right,ro);fig,ax=plt.subplots(figsize=(15,11))
 def half(t,age,yf,side):
  def rec(c):
   if not c.clades:return
   ys=[yf(x) for x in c.clades];xx=side*(1+age[c]);ax.plot([xx]*2,[min(ys),max(ys)],c='black',lw=.8)
   for q in c.clades:ax.plot([side*(1+age[q]),xx],[yf(q)]*2,c='black',lw=.8);rec(q)
  rec(t.root)
 half(left,al,fl,-1);half(right,ar,fr,1)
 for n in lo:ax.plot([-.95,.95],[yl[n],yr[n]],c=COL[lineage(n)],alpha=.45,lw=.7);ax.text(0,(yl[n]+yr[n])/2,n,ha='center',va='center',fontsize=7,c=COL[lineage(n)],bbox=dict(fc='white',ec='none',pad=.1))
 ax.axis('off');ax.set_ylim(0,len(lo)+1);ax.set_title(title,fontsize=14);save(fig,name)

M={};S={};CLOUD={};FREQ=[];CHAINFREQ=[];TOPO=[];COUNTS={};N={}
for s in ['zhang','givnish']:
 path=IN/f'{s}_unconstrained_final_mcc.tree';M[s]=mcc(path);S[s]=support_map(path);CLOUD[s],COUNTS[s],N[s],fr,cfr,topo=sampled_and_frequencies(s);FREQ.append(fr);CHAINFREQ.extend(cfr);TOPO.extend(topo)
 order=[q.name for q in M[s].get_terminals()];fig,ax=plt.subplots(figsize=(12,10));y,yf,age,H=draw(ax,M[s],f'Lineage-unconstrained RAD-seq MCC chronogram — {s.capitalize()}\nnode labels are posterior clade probabilities',order,S[s])
 ax.axvspan(0,.0117,color='#fff4bd',alpha=.25,zorder=-5);ax.axvspan(.0117,2.58,color='#d9ecff',alpha=.20,zorder=-5);ax.axvspan(2.58,5.333,color='#ffe0bf',alpha=.18,zorder=-5)
 diag=pd.read_csv(IN/f'{s}_unconstrained_scalar_diagnostics.tsv',sep='\t').set_index('parameter')
 for label,param,tips in [('PEQ–ingroup root','treeHeight',set(order)),('29-tip ingroup crown','mrca.age(ingroupTaxa)',set(order)-{'PEQ'})]:
  c=M[s].common_ancestor(list(tips));r=diag.loc[param];ax.plot([r.q025,r.q975],[yf(c),yf(c)],color='#2878a8',lw=6,alpha=.45,zorder=5);ax.plot(r['median'],yf(c),'o',ms=4,color='#124e70',zorder=6);ax.text(r.q975+.04,yf(c)+.14,f'{label}: {r["median"]:.2f} [{r.q025:.2f}–{r.q975:.2f}]',fontsize=7,color='#124e70')
 ax.set_xlim(diag.loc['treeHeight','q975']*1.10,-diag.loc['treeHeight','q975']*.12)
 save(fig,f'{s}_unconstrained_final_mcc_chronogram')
 fig,ax=plt.subplots(figsize=(12,10))
 for t in CLOUD[s]:
  _,yf2,a2,_=layout(t,order)
  def rec(c):
   if not c.clades:return
   ys=[yf2(x) for x in c.clades];ax.plot([a2[c]]*2,[min(ys),max(ys)],c='#5799c6',alpha=.025,lw=.4)
   for x in c.clades:ax.plot([a2[c],a2[x]],[yf2(x)]*2,c='#5799c6',alpha=.025,lw=.4);rec(x)
  rec(t.root)
 draw(ax,M[s],f'Lineage-unconstrained posterior tree cloud — {s.capitalize()}',order,None);save(fig,f'{s}_unconstrained_posterior_tree_cloud')
 with open(IN/f'{s}_unconstrained_100_posterior_full.newick','w') as a,open(IN/f'{s}_unconstrained_100_posterior_29tip.newick','w') as b:
  for t in CLOUD[s]:Phylo.write(t,a,'newick');u=Phylo.read(io.StringIO(t.format('newick')),'newick');u.prune('PEQ');Phylo.write(u,b,'newick')
 u=Phylo.read(io.StringIO(M[s].format('newick')),'newick');u.prune('PEQ');Phylo.write(u,IN/f'{s}_unconstrained_final_mcc_29tip.newick','newick')
 rows=[];dep=M[s].depths();HH=max(dep[q] for q in M[s].get_terminals())
 for i,c in enumerate(M[s].get_nonterminals()):
  k=frozenset(q.name for q in c.get_terminals());rows.append({'node':f'N{i+1}','n_descendants':len(k),'descendants':','.join(sorted(k)),'mcc_age':HH-dep[c],'posterior':S[s].get(k,1),'sampled_posterior_frequency':COUNTS[s].get(k,0)/N[s]})
 pd.DataFrame(rows).to_csv(IN/f'{s}_unconstrained_node_ages_and_support.tsv',sep='\t',index=False)
pd.DataFrame(FREQ).to_csv(IN/'unconstrained_lineage_monophyly_frequencies.tsv',sep='\t',index=False)
pd.DataFrame(CHAINFREQ).to_csv(IN/'unconstrained_lineage_monophyly_frequencies_by_chain.tsv',sep='\t',index=False)
pd.DataFrame(TOPO).to_csv(IN/'unconstrained_topology_agreement_by_chain.tsv',sep='\t',index=False)
tangle(M['zhang'],M['givnish'],'zhang_vs_givnish_unconstrained_tanglegram','Lineage-unconstrained MCC topologies: Zhang versus Givnish')
# Cross-calibration clade support.
ks=set(COUNTS['zhang'])|set(COUNTS['givnish']);pd.DataFrame([{'n_descendants':len(k),'descendants':','.join(sorted(k)),'zhang_frequency':COUNTS['zhang'].get(k,0)/N['zhang'],'givnish_frequency':COUNTS['givnish'].get(k,0)/N['givnish']} for k in sorted(ks,key=lambda x:(len(x),sorted(x))) ]).to_csv(IN/'zhang_vs_givnish_unconstrained_clade_frequencies.tsv',sep='\t',index=False)
