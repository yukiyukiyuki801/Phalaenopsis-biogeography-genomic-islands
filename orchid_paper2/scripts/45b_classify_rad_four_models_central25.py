#!/usr/bin/env python3
"""Four-model RAD classification with allopatric dXY restricted to Q37.5-Q62.5."""
import argparse, csv, math
from collections import Counter, defaultdict

def q(values, p):
    x=sorted(float(v) for v in values if v not in ("","NA","nan"))
    return x[max(0, math.ceil(p*len(x))-1)]

ap=argparse.ArgumentParser()
ap.add_argument("--input",required=True); ap.add_argument("--output",required=True)
ap.add_argument("--thresholds",required=True); ap.add_argument("--summary",required=True)
a=ap.parse_args()
with open(a.input) as f:
    reader=csv.DictReader(f,delimiter="\t"); fields=reader.fieldnames; rows=list(reader)
by=defaultdict(list)
for r in rows: by[r["pair"]].append(r)
thr={}
for pair,rr in by.items():
    for metric in ("hudson_fst","dxy"):
        vals=[r[metric] for r in rr]
        thr[pair,metric]={p:q(vals,p) for p in (.25,.375,.625,.75)}
    g1,g2=pair.split("-")
    for g,m in ((g1,"pi_group1"),(g2,"pi_group2")):
        vals=[r[m] for r in rr]
        thr[pair,"pi_"+g]={p:q(vals,p) for p in (.25,.75)}

def classify(fst,dxy,pi,t,p):
    if fst>=t[.75] and dxy>=p[.75] and pi<=p["pi25"]: return "divergence_with_gene_flow"
    if fst>=t[.75] and p[.375]<dxy<p[.625] and pi<=p["pi25"]: return "allopatric_selection"
    if fst>=t[.75] and dxy<=p[.25] and pi<=p["pi25"]: return "recurrent_selection"
    if fst<=t[.25] and dxy>=p[.75] and pi>=p["pi75"]: return "balancing_selection"
    return "unclassified"

for r in rows:
    pair=r["pair"]; g1,g2=pair.split("-"); fst=float(r["hudson_fst"]); dxy=float(r["dxy"])
    ft=thr[pair,"hudson_fst"]; dt=thr[pair,"dxy"]
    models=[]
    for g,m in ((g1,"pi_group1"),(g2,"pi_group2")):
        pt=thr[pair,"pi_"+g]
        params=dict(dt); params.update(pi25=pt[.25],pi75=pt[.75])
        models.append((g,classify(fst,dxy,float(r[m]),ft,params)))
    assigned=[m for _,m in models if m!="unclassified"]
    r["four_model_class_central25_dxy"]=assigned[0] if assigned else "unclassified"
    r["model_species1_central25_dxy"]=models[0][1]
    r["model_species2_central25_dxy"]=models[1][1]
    r["supporting_species_central25_dxy"]=",".join(g for g,m in models if m!="unclassified") or "none"
    r["n_supporting_species_central25_dxy"]=sum(m!="unclassified" for _,m in models)

extra=["four_model_class_central25_dxy","model_species1_central25_dxy",
       "model_species2_central25_dxy","supporting_species_central25_dxy",
       "n_supporting_species_central25_dxy"]
with open(a.output,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields+extra,delimiter="\t");w.writeheader();w.writerows(rows)
with open(a.thresholds,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t");w.writerow(["pair","metric","q25","q37.5","q62.5","q75"])
    for pair in sorted(by):
        for metric in ("hudson_fst","dxy"):
            z=thr[pair,metric];w.writerow([pair,metric,z[.25],z[.375],z[.625],z[.75]])
with open(a.summary,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t");w.writerow(["pair","model","n","fraction_all_loci","fraction_classified_loci"])
    order=("divergence_with_gene_flow","allopatric_selection","recurrent_selection","balancing_selection","unclassified")
    for pair,rr in sorted(by.items()):
        c=Counter(r["four_model_class_central25_dxy"] for r in rr); classified=len(rr)-c["unclassified"]
        for model in order:
            w.writerow([pair,model,c[model],c[model]/len(rr),
                        c[model]/classified if model!="unclassified" and classified else "NA"])
