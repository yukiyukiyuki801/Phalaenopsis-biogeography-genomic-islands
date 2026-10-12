#!/usr/bin/env python3
"""Build a 29-population + P. equestris complete-site RAD supermatrix."""
from __future__ import annotations
import argparse, csv, gzip
from collections import defaultdict
from pathlib import Path

IUPAC = {frozenset("A"):"A",frozenset("C"):"C",frozenset("G"):"G",frozenset("T"):"T",
         frozenset("AG"):"R",frozenset("CT"):"Y",frozenset("GC"):"S",frozenset("AT"):"W",
         frozenset("GT"):"K",frozenset("AC"):"M"}

def load_clusters(path, gap=500):
    intervals, bounds = defaultdict(list), {}
    cluster=0; last_chrom=None; last_end=None
    with open(path) as fh:
        for line in fh:
            chrom,s,e=line.rstrip().split("\t")[:3]; s=int(s); e=int(e)
            if chrom != last_chrom or last_end is None or s-last_end > gap: cluster += 1
            intervals[chrom].append((s+1,e,cluster))
            bounds.setdefault(cluster,[chrom,s+1,e]); bounds[cluster][2]=max(bounds[cluster][2],e)
            last_chrom,last_end=chrom,e
    return intervals,bounds

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--vcf',required=True); ap.add_argument('--bed',required=True)
    ap.add_argument('--popmap',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--site-pop-call',type=float,default=.70)
    ap.add_argument('--min-locus-sites',type=int,default=100)
    ap.add_argument('--max-tip-missing',type=float,default=.50)
    a=ap.parse_args(); out=Path(a.outdir); (out/'loci').mkdir(parents=True,exist_ok=True)
    smap={}; lineage={}
    with open(a.popmap) as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            smap[r['sample']]=r['population']; lineage[r['population']]=r['major_lineage']
    pops=sorted(set(smap.values()),key=lambda p:(lineage[p],p)); pop_samples=defaultdict(list)
    for s,p in smap.items(): pop_samples[p].append(s)
    intervals,bounds=load_clusters(a.bed); cursor=defaultdict(int)
    seq=defaultdict(lambda:{p:[] for p in pops+['PEQ']}); coords=defaultdict(list)
    with gzip.open(a.vcf,'rt') as fh:
        for line in fh:
            if line.startswith('##'): continue
            if line.startswith('#CHROM'):
                header=line.rstrip().split('\t'); col={s:i for i,s in enumerate(header)}
                missing=set(smap)-set(col)
                if missing: raise SystemExit('Samples absent from VCF: '+','.join(sorted(missing)))
                continue
            z=line.rstrip().split('\t'); chrom=z[0]; pos=int(z[1]); ref=z[3].upper(); alt=z[4].upper()
            if len(ref)!=1 or ref not in 'ACGT' or (alt!='.' and (len(alt)!=1 or alt not in 'ACGT')): continue
            iv=intervals[chrom]; j=cursor[chrom]
            while j<len(iv) and pos>iv[j][1]: j+=1
            cursor[chrom]=j
            if j==len(iv) or pos<iv[j][0]: continue
            cid=iv[j][2]; alleles=[ref] if alt=='.' else [ref,alt]
            for p in pops:
                observed=[]; called=0
                for sample in pop_samples[p]:
                    gt=z[col[sample]].split(':',1)[0].replace('|','/')
                    if gt in ('.','./.'): continue
                    aa=gt.split('/')
                    try: bases=[alleles[int(x)] for x in aa if x!='.']
                    except (ValueError,IndexError): bases=[]
                    if len(bases)==2: called+=1; observed.extend(bases)
                if called/len(pop_samples[p]) < a.site_pop_call or not observed:
                    base='N'
                else:
                    base=IUPAC.get(frozenset(observed),'N')
                seq[cid][p].append(base)
            seq[cid]['PEQ'].append(ref)
            coords[cid].append(pos)

    kept=[]; concat={p:[] for p in pops+['PEQ']}; offset=0
    rejected_short=0; rejected_missing=0; missing_fail_by_pop=defaultdict(int)
    for cid in sorted(seq):
        n=len(coords[cid])
        miss={p:seq[cid][p].count('N')/n for p in pops}
        if n<a.min_locus_sites:
            rejected_short += 1; continue
        failing=[p for p,v in miss.items() if v>a.max_tip_missing]
        if failing:
            rejected_missing += 1
            for p in failing: missing_fail_by_pop[p] += 1
            continue
        chrom,start,end=bounds[cid]; locus=f'RAD{cid:05d}'
        with open(out/'loci'/f'{locus}.fasta','w') as f:
            for p in pops+['PEQ']: f.write(f'>{p}\n{"".join(seq[cid][p])}\n')
        for p in concat: concat[p].extend(seq[cid][p])
        kept.append([locus,chrom,start,end,n,offset+1,offset+n,max(miss.values()),sum(miss.values())/len(miss)])
        offset+=n
    if not kept: raise SystemExit('No loci passed filters')
    with open(out/'rad29_peq_concat.fasta','w') as f:
        for p in pops+['PEQ']: f.write(f'>{p}\n{"".join(concat[p])}\n')
    with open(out/'locus_partitions.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['locus','chrom','genomic_start','genomic_end','callable_sites','alignment_start','alignment_end','max_tip_missing','mean_tip_missing']); w.writerows(kept)
    with open(out/'taxa.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['taxon','lineage','n_individuals'])
        for p in pops: w.writerow([p,lineage[p],len(pop_samples[p])])
        w.writerow(['PEQ','OUTGROUP',1])
    with open(out/'matrix_summary.tsv','w') as f:
        f.write(f'populations\t{len(pops)}\noutgroup\tPEQ\nindividuals_collapsed\t{len(smap)}\ncallable_bed_clusters\t{len(bounds)}\nclusters_with_filtered_vcf_sites\t{len(seq)}\nrejected_fewer_than_min_sites\t{rejected_short}\nrejected_tip_missingness\t{rejected_missing}\nretained_loci\t{len(kept)}\nalignment_sites\t{offset}\nsite_population_call_threshold\t{a.site_pop_call}\nmin_locus_callable_sites\t{a.min_locus_sites}\nmax_tip_locus_missing\t{a.max_tip_missing}\n')
    with open(out/'locus_filter_attrition.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['stage','n_loci'])
        w.writerows([['callable_bed_clusters',len(bounds)],['clusters_with_filtered_vcf_sites',len(seq)],['rejected_fewer_than_min_sites',rejected_short],['rejected_tip_missingness',rejected_missing],['retained',len(kept)]])
    with open(out/'missingness_rejection_by_population.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['population','lineage','loci_with_tip_missing_gt_threshold'])
        for p in sorted(pops,key=lambda p:(-missing_fail_by_pop[p],p)): w.writerow([p,lineage[p],missing_fail_by_pop[p]])

if __name__=='__main__': main()
