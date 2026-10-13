# Analysis workflows retained in the manuscript

All paths below are relative to `orchid_paper2/` unless stated otherwise.
These are the executed analysis choices; batch resources are examples.

## 1. Read processing and invariant-inclusive calls

Inputs are 173 demultiplexed paired-end RAD datasets from PRJNA1320963 and the
versioned *P. equestris* reference. Scripts 00–05 download, trim, map, and call
reference-confidence genotypes. Trimmomatic 0.39 uses SLIDINGWINDOW:4:15 and
MINLEN:36; these reconstruct the published trimming framework rather than
claiming an identical lost command. BWA 0.7.17 and SAMtools 1.21 generate sorted,
indexed BAMs. GATK 4.6.2.0 performs HaplotypeCaller, GenomicsDBImport, and
GenotypeGVCFs with invariant calls on empirically callable intervals.

Scripts 34–36 divide joint calling into disjoint shards and merge them. The
curated script 38 retains only the filter actually used for the reported VCF
landscapes: biallelic SNP/invariant records and total INFO/DP ≥30. Pair-specific
filters are applied downstream. The canonical single-base target has 512,526
positions (43,262 variable; 469,264 invariant), distinct from 518,090 VCF records.
Full all-sites VCFs are generated from public reads rather than shipped in Git.

## 2. Individual relationships and contemporary migration — Figure 1, S1

Use the published 170-sample, 6,596-SNP VCF at the repository root, not an
additional neutral/MAF subset. `scripts/prepare_network_distances.py` generates
pairwise mean absolute dosage difference/2 at jointly called SNPs and a NEXUS
distance input. SplitsTree4 4.19.2 infers NeighborNet with its native EqualAngle
layout. The saved individual distance/NEXUS files are included.

`71_prepare_ba3_snps_input.py` converts this VCF to the final 18 units:

```bash
python orchid_paper2/scripts/71_prepare_ba3_snps_input.py \
  --vcf phalaenopsis_mac20260730/Phalaenopsis_eq_LD_filter.vcf \
  --mapping orchid_paper2/18_contemporary_migration/final_18_population_classification.tsv \
  --metadata phalaenopsis_mac20260730/full_170samples_popmap.csv \
  --outdir orchid_paper2/18_contemporary_migration/inputs
```

This packaged converter exports only all LD-pruned SNPs. Scripts 72/73 tune and
run BA3-SNPs 1.1. Five runs use 10 million iterations, one million burn-in,
sampling every 1,000. The source-to-recipient parameter is m per generation.
The pooled summary includes all 306 directions, median m, credible intervals,
Pr(m >0.01), approximate lag-one pooled ESS, and ordinary variance-based R-hat.
The retained summarizer is `summarize_ba3_chains.py`. This ESS approximation
does not establish convergence or replace full-autocorrelation diagnostics.
The map selects the strongest **25 directions**, not
25%, with median m ≥0.01 and Pr(m >0.01) ≥0.50. Dashed arrows flag R-hat >1.05.
Species is encoded by recipient color and width by numeric median m. Tables
S2–S3 document assignments, including the retained two-sample LAM_DO exception.

The migration map requires the complete Natural Earth 1:10m country dataset,
including every polygon of Indonesia. The packaged renderer fails explicitly
when it is absent and uses equal longitude/latitude axis scaling:

```bash
curl -fL https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_0_countries.geojson \
  -o orchid_paper2/18_contemporary_migration/inputs/ne_10m_admin_0_countries.geojson
```

## 3. Single-copy orthologue calibration — Figure S2

Scripts 70–75 retrieve ten public taxa, assemble PAM RNA, predict PAM/PAP
proteins, and select primary proteins from annotated genomes (Table S4).
Script 81 uses OrthoFinder 2.5.4's DIAMOND/MCL/DendroBLAST workflow and STAG/STRIDE
species topology. There is no subsequent IQ-TREE or FastTree topology rerun.
Scripts 84/85 align and trim strict single-copy genes with MAFFT/trimAl. Dating
uses 353 genes and 169,668 amino-acid positions, not all multicopy orthogroups.

MCMCTree (PAML 4.10.7) uses WAG+F, five gamma categories, an independent-rates
clock, approximate likelihood, and time units of 100 Myr. Two chains per
scenario each use 100,000 burn-in states, sampling every 20 states for 100,000
retained samples. The two **separate** scenarios each contain three soft
secondary bounds, not a combined six-calibration run. Tables S5–S8 record
bounds, convergence, ages, and the lognormal priors transferred to RAD dating.

## 4. Unconstrained population RAD dating — Figures 2, S3–S7

Script 86 builds 29 population-consensus sequences from 94 selected samples
and adds the true reference outgroup. It retains invariant sites, ≥70% site
callability per population, ≥100 retained bases per locus, and ≤50% missingness
per population. The final alignment contains 274 loci and 79,225 bases.

Run the supplied `16_rad_beast/xml/*combined_rateprior_unconstrained_60m.xml`
directly in BEAST 2.7.4, with three independent seeds per calibration. Settings
are HKY+Γ4, empirical frequencies, IUPAC likelihoods, an uncorrelated lognormal
clock, a Yule prior, two secondary age priors, and the broad external nuclear
rate prior. Only ingroup monophyly relative to the outgroup is imposed.
The eight plastid loci were not separately partitioned; the alignment is not
strictly nuclear. Each chain has 60 million states, logging every 5,000 and
discarding six million. Scripts 104–110 summarize diagnostics, MCC trees with
common-ancestor heights, and 100 posterior trees separately per scenario.
Low ESS and Zhang chain disagreement remain in the diagnostic table (S9).

## 5. Historical ranges and stochastic mapping — Figures 2–4a, S8–S11

The ten sampled regions use codes A/B/C/E/F/G/H/I/J/K, with no sampled D tip.
The mapping and geometry files are retained under `19_biogeobears/bioregions`
and supply only the reported geographic inputs. Script 96 builds vector polygons
and tip assignments, not raster distributions or regionalization. Codes, colors,
adjacency, and tip geography are explicit inputs.

Script 103 refreshes MCC and 100 posterior-tree inputs from the reported
unconstrained summaries while retaining the supplied ten-area, graph-connected
range definition. It does not use an earlier constrained tree.

For each calibration, set `BGB_BASE` to
`$ORCHID_PROJECT_ROOT/orchid_paper2/19_biogeobears/unconstrained_29tip/<calibration>/max_range_3`.
Script 104 runs one model on one tree:

```bash
export ORCHID_PAPER2_ROOT="$ORCHID_PROJECT_ROOT"
export BGB_BASE="$ORCHID_PROJECT_ROOT/orchid_paper2/19_biogeobears/unconstrained_29tip/givnish/max_range_3"
Rscript orchid_paper2/scripts/104_run_biogeobears_unified29.R mcc DIVALIKEJ 3
```

Repeat for DEC, DECJ, DIVALIKE, DIVALIKEJ, BAYAREALIKE, and BAYAREALIKEJ on the
MCC and posterior_001–100 trees. There are 606 fits per calibration. Maximum
range size is **3**, with 10 single, 13 connected double, and 21 connected
triple ranges plus null. There is **no time stratification**. Script 106 ranks
AICc; script 107 generates 100 BSM histories under the selected MCC model.
The modal transition map (script 141) independently selects each node's modal
complete range. Its time bins are visualization bins, not model strata.
The PSA loop denotes retention, not a BSM sympatry assignment.

## 6. Cycle-specific climate association — Figures 4b, S12, S14

Script 111 reads de Boer et al. (2014) sea level and isotope values, identifies
111 sea-level cycles, and standardizes isotope values **within cycles**, not
over the pooled record. Standardized δ18O ≥0.75 is glacial, ≤0.25 interglacial,
and otherwise intermediate. Each of 100 posterior trees contributes its 28
internal ages. There are 49,999 common-fraction circular shifts. The main
figure reports Holm correction of the two exploratory focal glacial tests;
error bars are SD among posterior trees/null shifts. The sea-level comparison
and all category results are retained numerically. Raw versus standardized
curves are generated by script 114. These are not formal marine-isotope-stage
boundaries or 2,800 independent divergence events.

## 7. Population-lineage branching models — Figures S13, S15; Table 2

Scripts 132–138 prepare rounded ultrametric MCC inputs, BAMM controls, and
RPANDA fits. BAMM uses 24 independent runs: two calibrations × expected-shift
priors 0.5/1/2 × four runs, each with four coupled chains, 10 million generations,
logging every 2,000, and 25% burn-in. Priors are tree-specific. The coupled
chains are not independent replicates.

RPANDA fits exactly nine models per calibration: N1 (constant branching and
loss), E1–E4, and T1–T4. E1/T1 varies branching with zero loss; E2/T2 varies
branching with constant loss; E3/T3 has constant branching and variable loss;
E4/T4 varies both. E models use standardized InfTemp and T models time. Six
starts per fit, crown conditioning, sampling fraction one, spline df 30, and
dt 0.001 Myr are used. AICc comparison retains constant/time alternatives.
Rates concern sampled population lineages, not directly species extinction.

## 8. RAD differentiation signatures — Figures 5, S16–S20

Script 42 calculates finite-sample-corrected π, allele-frequency dXY, Hudson
FST, and approximate short-locus Tajima's D. Loci require ≥100 shared callable
positions, ≥70% called individuals per lineage, and ≥3 informative SNPs.
Scripts 45b/46b/47b classify and plot the original VCF data. Upper/lower
quartiles define high/low, and intermediate dXY is Q37.5–Q62.5. π is evaluated
separately in each species using OR logic. Signatures are divergence-like
(high FST/high dXY/low π), allopatric-like (high FST/intermediate dXY/low π),
recurrent-like (high FST/low dXY/low π), and balancing-like (low FST/high dXY/high π).

Scripts 117–125 estimate folded SFS-based indices with ANGSD 0.941. Minimum
base/mapping qualities are 20/30, proper pairs are required, and per-group
callability is ≥70%; upper depth uses the group-specific 99.5th percentile.
Script 127 compares the three-statistic quartile rule and the Shang et al.
(2023) FST/dXY-only 5%-tail rule. Class ratios are candidate patterns, not a
validated partition of selection mechanisms. Dataset counts are in Table S1.

## 9. Shared outlier loci — Figures S22–S23

Scripts 51/52 use 576 shared nuclear loci and top/bottom 5% sets. FST/dXY sets
are lineage-pair comparisons; π and Tajima's D sets are single lineages. Both
tails of D are retained. Two-way intersections include triple-shared loci.
There are 49,999 within-stratum permutations preserving set sizes in 4 ×4
callable-length/information strata. Significance uses one-sided empirical p
and BH q across 20 tests, not simply crossing the upper 97.5% null line.

## 10. PAP–PSA local genealogies and annotation — Figure 6, S21; Table S10

Scripts 100/102/106/107 prepare invariant-inclusive sequences and infer IQ-TREE
3.0.1 HKY+F+Γ4 trees for all 650 nuclear loci, rooted with the true reference.
All 650 trees were completed; 647 loci yield normalized metrics. Script
`108_summarize_all650_genealogies.py` retains only the complete-data analysis.
The original BH family of 12 contrasts is preserved; only relative mean
depth and nearest-neighbor mixing are plotted in the main figure. Distances
use 500 balanced seven-versus-seven resamples. Relative depth is ingroup mean
pairwise distance/reference divergence, not an absolute age or a tree of
height one. Five read-phased, tip-aligned trees provide qualitative illustrations
only; no separate read-phased sensitivity analysis is included.

`143_annotate_balancing_candidate_genes.py` is the annotation-only portion of
the original workflow. It matches transcript overlaps/nearest mRNAs within
10 kb in the versioned reference, collapses isoforms to genes, and joins
Arabidopsis homolog assignments. It does not contact g:Profiler or conduct GO
enrichment. The functional table contains nine orchid genes and eight homologs.
