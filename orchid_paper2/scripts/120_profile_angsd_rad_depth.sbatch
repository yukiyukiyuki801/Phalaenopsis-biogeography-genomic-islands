#!/usr/bin/env bash
#SBATCH --job-name=angsd_depth
#SBATCH --partition=day,scavenge
#SBATCH --array=1-3
#SBATCH --cpus-per-task=30
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/depth_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/depth_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load angsd/0.940-GCC-12.2.0
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
WORK="$ROOT/21_angsd_rad_landscapes"
REF="$ROOT/01_reference/GCF_001263595.1_ASM126359v1_genomic.fna"
case "$SLURM_ARRAY_TASK_ID" in
  1) GROUP=PAM; MININD=62 ;;
  2) GROUP=PAP; MININD=55 ;;
  3) GROUP=PSA; MININD=5 ;;
esac

angsd \
  -bam "$WORK/inputs/${GROUP}.bamlist" \
  -ref "$REF" -sites "$WORK/inputs/callable_sites.txt" \
  -rf "$WORK/inputs/target_rad_regions.txt" \
  -doCounts 1 -doDepth 1 -maxDepth 100000 \
  -minMapQ 30 -minQ 20 -only_proper_pairs 1 -remove_bads 1 \
  -minInd "$MININD" -nThreads "$SLURM_CPUS_PER_TASK" \
  -out "$WORK/depth/${GROUP}"

test -s "$WORK/depth/${GROUP}.depthGlobal"

