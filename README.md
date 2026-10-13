# Phalaenopsis biogeography and genomic islands

Code and processed inputs for analyses of the *Phalaenopsis amabilis* complex
retained in the manuscript. PAM = *P. amabilis*, PAP = *P. aphrodite*, and
PSA = *P. sanderiana*. Previously published sampling and laboratory procedures
are described by Sun et al. (2026).

## Scope

The repository includes individual NeighborNet, contemporary migration,
single-copy orthologue dating, unconstrained population RAD dating, historical
biogeography and BSM, climate associations, population-lineage branching models,
RAD differentiation signatures, outlier recurrence, local genealogies, and
candidate-gene annotation. Only the Givnish and Zhang unconstrained age-plus-rate
prior scenarios are included for RAD dating. Herodotools inference and raster
preparation, ENM, SNAP(P)ER, constrained RAD
trees, fixed-rate experiments, and GO enrichment are not included.

## Organization

The `orchid_paper2/` layout is retained so scripts can find their inputs using
relative paths. Scripts are in `orchid_paper2/scripts/`, processed data in the
corresponding analysis directories, and separately numbered supplementary
tables in `tables/`. `docs/analyses.md` identifies execution order, inputs,
settings, and the source of each manuscript display. `code_provenance.tsv`
records copied/adapted sources, and `checksums.sha256` records release contents.

## Getting started

```bash
git clone https://github.com/yukiyukiyuki801/Phalaenopsis-biogeography-genomic-islands.git
cd Phalaenopsis-biogeography-genomic-islands
export ORCHID_PROJECT_ROOT="$PWD"
python scripts/check_package.py
```

Run commands from the repository root unless an analysis explicitly states
otherwise. Create required output directories before running batch examples.
Scheduler partitions, module names, CPU counts, and memory requests in `.sbatch`
files describe the original computing environment and must be adapted locally.
The repository does not submit jobs automatically.

Use the final supplied BEAST XMLs directly; generating them does not require an
earlier constrained analysis. Long MCMC fits and sequencing-data reconstruction
were not rerun as part of code packaging. Validation checks file integrity,
Python/shell syntax, key scientific counts, and selected preprocessing outputs;
it is not a claim of end-to-end execution on a clean machine. See
`VALIDATION.md` for the checks actually completed.

## Data access

- Raw demultiplexed RAD reads: NCBI BioProject **PRJNA1320963**; the run manifest
  and sample identifiers are included. FASTQ, BAM, gVCF, and full all-sites VCF
  files are not committed; use the read-processing and calling workflow.
- Mapping/reference outgroup: *P. equestris* assembly **GCF_001263595.1**.
  Reference checksum/index metadata are included; reference preparation is
  provided in the workflow.
- Published SNP input: the 170-individual, 6,596-SNP LD-pruned VCF and population
  metadata are included for NeighborNet and BA3-SNPs. No additional MAF filter
  or neutral-only restriction is used.
- Public orthologue inputs: accessions and input types for ten taxa are listed
  in `tables/Table_S4.tsv`. The 353-gene concatenated protein alignment,
  calibration trees, and control settings are included.
- RAD dating: 29 population consensuses plus the real reference outgroup,
  274 loci, 79,225 positions, two final XMLs, MCC trees, 100 posterior trees per
  calibration, and diagnostics are included.
- GFF annotation: `GEA/orchid.gff.gz` contains the reference annotation used for
  coordinate matching. Decompress with `gzip -dk GEA/orchid.gff.gz` before
  rerunning candidate annotation. The homolog-assignment table is included.

Input sets differ across analyses. In particular, the original VCF landscape
classification includes plastid loci, whereas genealogy and outlier recurrence
use nuclear-only sets. `tables/Table_S1.tsv` records these distinctions.

## Software

Python dependencies: NumPy, pandas, SciPy, matplotlib, seaborn, Biopython,
PyMuPDF, requests, and pysam where required. Geographic plotting additionally
uses GDAL/OGR. R dependencies include data.table, ggplot2, ape, phangorn,
BioGeoBEARS, BAMMtools, RPANDA, and supporting dependencies. Command-line tools
include BWA, SAMtools, Trimmomatic, GATK, BCFtools, SRA Toolkit, SplitsTree,
BA3-SNPs, OrthoFinder/DIAMOND, RNA-SPAdes, TransDecoder, MAFFT, trimAl, PAML,
BEAST2/TreeAnnotator, ANGSD, IQ-TREE, WhatsHap, and BAMM. Analysis documents
identify versions recovered from the completed runs.

No software binaries, environments, or third-party package source trees are
redistributed. No new license has been assigned by this export; contact the
repository owner regarding reuse of original project code. Third-party data
and tools remain subject to their original terms.
