#!/usr/bin/env bash
#SBATCH --job-name=angsd_safmaf
#SBATCH --partition=day,scavenge
#SBATCH --array=1-3
#SBATCH --cpus-per-task=30
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/safmaf_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/safmaf_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load angsd/0.940-GCC-12.2.0
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
WORK="$ROOT/21_angsd_rad_landscapes"
REF="$ROOT/01_reference/GCF_001263595.1_ASM126359v1_genomic.fna"
case "$SLURM_ARRAY_TASK_ID" in
  1) GROUP=PAM ;;
  2) GROUP=PAP ;;
  3) GROUP=PSA ;;
esac
MININD=$(awk -v g="$GROUP" 'BEGIN{FS="\t"}$1==g{print $2}' "$WORK/inputs/depth_filters.tsv")
MAXDEPTH=$(awk -v g="$GROUP" 'BEGIN{FS="\t"}$1==g{print $6}' "$WORK/inputs/depth_filters.tsv")
[[ "$MININD" =~ ^[0-9]+$ && "$MAXDEPTH" =~ ^[0-9]+$ ]]

# REF is supplied as a temporary orientation for the SAF, which is folded in
# the next stage.  It is not interpreted as a true ancestral sequence.
angsd \
  -bam "$WORK/inputs/${GROUP}.bamlist" \
  -ref "$REF" -anc "$REF" \
  -sites "$WORK/inputs/callable_sites_alleles.txt" \
  -rf "$WORK/inputs/target_rad_regions.txt" \
  -GL 1 -doCounts 1 -doSaf 1 -doMajorMinor 3 -doMaf 1 \
  -minMapQ 30 -minQ 20 -only_proper_pairs 1 -remove_bads 1 \
  -minInd "$MININD" -setMaxDepth "$MAXDEPTH" \
  -nThreads "$SLURM_CPUS_PER_TASK" \
  -out "$WORK/saf/${GROUP}"

test -s "$WORK/saf/${GROUP}.saf.idx"
test -s "$WORK/saf/${GROUP}.mafs.gz"
ln -sfn "../saf/${GROUP}.mafs.gz" "$WORK/maf/${GROUP}.mafs.gz"

