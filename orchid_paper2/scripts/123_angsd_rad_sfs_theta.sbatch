#!/usr/bin/env bash
#SBATCH --job-name=angsd_theta
#SBATCH --partition=day,scavenge
#SBATCH --array=1-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/theta_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/theta_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load angsd/0.940-GCC-12.2.0
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
WORK="$ROOT/21_angsd_rad_landscapes"
case "$SLURM_ARRAY_TASK_ID" in
  1) GROUP=PAM ;;
  2) GROUP=PAP ;;
  3) GROUP=PSA ;;
esac

realSFS "$WORK/saf/${GROUP}.saf.idx" -fold 1 -cores "$SLURM_CPUS_PER_TASK" \
  > "$WORK/sfs/${GROUP}.folded.sfs"
realSFS saf2theta "$WORK/saf/${GROUP}.saf.idx" \
  -sfs "$WORK/sfs/${GROUP}.folded.sfs" -fold 1 \
  -outname "$WORK/theta/${GROUP}"
thetaStat print "$WORK/theta/${GROUP}.thetas.idx" | gzip -c \
  > "$WORK/theta/${GROUP}.thetas.tsv.gz"

test -s "$WORK/theta/${GROUP}.thetas.tsv.gz"

