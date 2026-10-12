#!/usr/bin/env bash
#SBATCH --job-name=angsd_fst
#SBATCH --partition=day,scavenge
#SBATCH --array=1-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/fst_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/fst_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load angsd/0.940-GCC-12.2.0
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
WORK="$ROOT/21_angsd_rad_landscapes"
case "$SLURM_ARRAY_TASK_ID" in
  1) G1=PAM; G2=PAP ;;
  2) G1=PAM; G2=PSA ;;
  3) G1=PAP; G2=PSA ;;
esac
PAIR=${G1}-${G2}

realSFS "$WORK/saf/${G1}.saf.idx" "$WORK/saf/${G2}.saf.idx" \
  -fold 1 -cores "$SLURM_CPUS_PER_TASK" > "$WORK/sfs/${PAIR}.folded.2dsfs"
realSFS fst index "$WORK/saf/${G1}.saf.idx" "$WORK/saf/${G2}.saf.idx" \
  -sfs "$WORK/sfs/${PAIR}.folded.2dsfs" -fold 1 \
  -fstout "$WORK/fst/${PAIR}"
realSFS fst print "$WORK/fst/${PAIR}.fst.idx" | gzip -c \
  > "$WORK/fst/${PAIR}.per_site.tsv.gz"
realSFS fst stats "$WORK/fst/${PAIR}.fst.idx" \
  > "$WORK/fst/${PAIR}.global.tsv"

test -s "$WORK/fst/${PAIR}.per_site.tsv.gz"

