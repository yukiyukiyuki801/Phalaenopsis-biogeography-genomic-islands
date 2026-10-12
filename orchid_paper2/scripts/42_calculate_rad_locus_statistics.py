#!/usr/bin/env python3
"""Calculate pairwise population-genetic summaries in callable RAD clusters.

Primary estimators:
  pi: finite-sample-corrected pairwise diversity per callable site
  dxy: between-population pairwise difference per callable site
  FST: Hudson ratio of sums, (sum(dxy)-sum(mean pi))/sum(dxy)
  Tajima's D: standard formula using minimum called chromosomes in the locus
"""
from __future__ import annotations
import argparse, csv, gzip, math
from collections import defaultdict
from pathlib import Path

GROUP_N = {"PAM": 88, "PAP": 78, "PSA": 7}
PAIRS = (("PAM", "PAP"), ("PAM", "PSA"), ("PAP", "PSA"))


def load_bed(path, gap=500):
    by_chrom, bounds = defaultdict(list), {}
    cluster = 0
    last_chrom, last_end = None, None
    with open(path) as f:
        for line in f:
            chrom, s, e = line.rstrip().split("\t")[:3]
            s, e = int(s), int(e)
            if chrom != last_chrom or last_end is None or s - last_end > gap:
                cluster += 1
            by_chrom[chrom].append((s + 1, e, cluster))
            if cluster not in bounds:
                bounds[cluster] = [chrom, s + 1, e]
            else:
                bounds[cluster][1] = min(bounds[cluster][1], s + 1)
                bounds[cluster][2] = max(bounds[cluster][2], e)
            last_chrom, last_end = chrom, e
    return by_chrom, bounds


def parse_info(s):
    return dict(x.split("=", 1) for x in s.split(";") if "=" in x)


def scalar(info, key):
    x = info.get(key, ".").split(",")[0]
    return float(x) if x not in ("", ".") else 0.0


def pi_site(ac, an):
    if an < 2:
        return math.nan
    p = ac / an
    return 2.0 * p * (1.0 - p) * an / (an - 1.0)


def tajima_d(theta_pi, segregating, n):
    n = int(n)
    if n < 2 or segregating < 2:
        return math.nan
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n*n + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2.0) / (a1*n) + a2 / (a1*a1)
    e1 = c1 / a1
    e2 = c2 / (a1*a1 + a2)
    den = math.sqrt(e1*segregating + e2*segregating*(segregating-1.0))
    return (theta_pi - segregating/a1) / den if den > 0 else math.nan


def fmt(x):
    return "NA" if not math.isfinite(x) else f"{x:.10g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--callable-bed", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--call-rate", type=float, default=.70)
    ap.add_argument("--min-bases", type=int, default=100)
    ap.add_argument("--min-snps", type=int, default=3)
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    intervals, bounds = load_bed(args.callable_bed)
    cursor = defaultdict(int)
    # pair/cluster -> bases, informative, sums plus population-specific fields
    dat = defaultdict(lambda: {"bases":0,"snps":0,"dxy_sum":0.0,
        "pi1_sum":0.0,"pi2_sum":0.0,"s1":0,"s2":0,"min_an1":10**9,
        "min_an2":10**9,"positions":[]})
    with gzip.open(args.vcf, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            z = line.rstrip().split("\t", 9)
            chrom, pos, ref, alt = z[0], int(z[1]), z[3], z[4]
            iv = intervals[chrom]; j = cursor[chrom]
            while j < len(iv) and pos > iv[j][1]: j += 1
            cursor[chrom] = j
            if j == len(iv) or pos < iv[j][0]:
                raise RuntimeError(f"Position outside callable BED: {chrom}:{pos}")
            cluster = iv[j][2]; info = parse_info(z[7])
            vals = {g:(scalar(info,f"AN_{g}"),scalar(info,f"AC_{g}")) for g in GROUP_N}
            for g1,g2 in PAIRS:
                an1,ac1=vals[g1]; an2,ac2=vals[g2]
                if an1 < 2*GROUP_N[g1]*args.call_rate or an2 < 2*GROUP_N[g2]*args.call_rate: continue
                p1,p2=ac1/an1,ac2/an2
                pi1,pi2=pi_site(ac1,an1),pi_site(ac2,an2)
                x=dat[(f"{g1}-{g2}",cluster)]
                x["bases"]+=1; x["dxy_sum"]+=p1*(1-p2)+p2*(1-p1)
                x["pi1_sum"]+=pi1; x["pi2_sum"]+=pi2
                x["min_an1"]=min(x["min_an1"],int(an1)); x["min_an2"]=min(x["min_an2"],int(an2))
                informative = alt != "." and 0 < ac1+ac2 < an1+an2
                if informative: x["snps"]+=1
                if alt != "." and 0 < ac1 < an1: x["s1"]+=1
                if alt != "." and 0 < ac2 < an2: x["s2"]+=1
                if alt != "." and informative: x["positions"].append((chrom,pos,ref,alt))

    fields=("pair locus_id chrom start end callable_bases informative_snps "
            "pi_group1 pi_group2 pi_pair_mean dxy hudson_fst segregating_group1 "
            "segregating_group2 min_an_group1 min_an_group2 tajima_d_group1 "
            "tajima_d_group2 tajima_d_pair_mean").split()
    rows=[]
    for (pair,cluster),x in dat.items():
        if x["bases"] < args.min_bases or x["snps"] < args.min_snps: continue
        g1,g2=pair.split("-"); b=x["bases"]
        pi1=x["pi1_sum"]/b; pi2=x["pi2_sum"]/b; dxy=x["dxy_sum"]/b
        fst=(x["dxy_sum"]-(x["pi1_sum"]+x["pi2_sum"])/2)/x["dxy_sum"] if x["dxy_sum"]>0 else math.nan
        td1=tajima_d(x["pi1_sum"],x["s1"],x["min_an1"])
        td2=tajima_d(x["pi2_sum"],x["s2"],x["min_an2"])
        tdmean=sum(v for v in (td1,td2) if math.isfinite(v))/sum(math.isfinite(v) for v in (td1,td2)) if any(math.isfinite(v) for v in (td1,td2)) else math.nan
        chrom,start,end=bounds[cluster]
        row=dict(pair=pair,locus_id=f"RAD{cluster:05d}",chrom=chrom,start=start,end=end,
            callable_bases=b,informative_snps=x["snps"],pi_group1=pi1,
            pi_group2=pi2,pi_pair_mean=(pi1+pi2)/2,dxy=dxy,hudson_fst=fst,
            segregating_group1=x["s1"],segregating_group2=x["s2"],
            min_an_group1=x["min_an1"],min_an_group2=x["min_an2"],
            tajima_d_group1=td1,tajima_d_group2=td2,tajima_d_pair_mean=tdmean,
            positions=x["positions"])
        rows.append(row)
    rows.sort(key=lambda r:(r["pair"],r["chrom"],r["start"]))

    metrics=("hudson_fst","dxy","pi_pair_mean","tajima_d_pair_mean")
    thresholds={}
    for pair in ("PAM-PAP","PAM-PSA","PAP-PSA"):
        rr=[r for r in rows if r["pair"]==pair]
        for m in metrics:
            v=sorted(r[m] for r in rr if math.isfinite(r[m]))
            thresholds[(pair,m)]=(v[max(0,math.ceil(.05*len(v))-1)],
                                  v[math.ceil(.95*len(v))-1])
    with open(out/"rad_locus_statistics.tsv","w",newline="") as f:
        extra=[f"top5_{m}" for m in metrics]+[f"bottom5_{m}" for m in metrics]+[f"directional_outlier_{m}" for m in metrics]
        w=csv.DictWriter(f,fieldnames=fields+extra,delimiter="\t")
        w.writeheader()
        for r in rows:
            q={k:(fmt(v) if isinstance(v,float) else v) for k,v in r.items() if k in fields}
            q.update({f"top5_{m}":int(math.isfinite(r[m]) and r[m]>=thresholds[(r['pair'],m)][1]) for m in metrics})
            q.update({f"bottom5_{m}":int(math.isfinite(r[m]) and r[m]<=thresholds[(r['pair'],m)][0]) for m in metrics})
            q.update({f"directional_outlier_{m}":int(math.isfinite(r[m]) and (r[m]>=thresholds[(r['pair'],m)][1] if m in ('hudson_fst','dxy') else r[m]<=thresholds[(r['pair'],m)][0])) for m in metrics})
            w.writerow(q)
    with open(out/"top5_thresholds.tsv","w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(("pair","metric","q05","q95","directional_tail"))
        for k,v in sorted(thresholds.items()): w.writerow((*k,fmt(v[0]),fmt(v[1]),"upper" if k[1] in ('hudson_fst','dxy') else "lower"))
    with gzip.open(out/"rad_snp_plot_table.tsv.gz","wt",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(("pair","locus_id","chrom","pos","ref","alt",*metrics,*[f"directional_outlier_{m}" for m in metrics]))
        for r in rows:
            for chrom,pos,ref,alt in r["positions"]:
                w.writerow((r["pair"],r["locus_id"],chrom,pos,ref,alt,*[fmt(r[m]) for m in metrics],*[int(math.isfinite(r[m]) and (r[m]>=thresholds[(r['pair'],m)][1] if m in ('hudson_fst','dxy') else r[m]<=thresholds[(r['pair'],m)][0])) for m in metrics]))
    with open(out/"analysis_parameters.tsv","w") as f:
        f.write(f"call_rate_per_population\t{args.call_rate}\nmin_callable_bases\t{args.min_bases}\nmin_informative_snps\t{args.min_snps}\nFST_estimator\tHudson ratio of sums\n")

if __name__ == "__main__": main()
