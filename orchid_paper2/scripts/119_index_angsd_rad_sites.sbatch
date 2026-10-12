#!/usr/bin/env bash
#SBATCH --job-name=angsd_index
#SBATCH --partition=day,scavenge
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/index_%j.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/21_angsd_rad_landscapes/logs/index_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load angsd/0.940-GCC-12.2.0
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
INPUT="$ROOT/21_angsd_rad_landscapes/inputs"

angsd sites index "$INPUT/callable_sites.txt"
angsd sites index "$INPUT/callable_sites_alleles.txt"
test -s "$INPUT/callable_sites.txt.idx"
test -s "$INPUT/callable_sites_alleles.txt.idx"

