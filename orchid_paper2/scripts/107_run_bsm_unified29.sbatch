#!/usr/bin/env bash
#SBATCH --job-name=bgb_u29_bsm
#SBATCH --partition=ycga,week
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/bgb_u29_bsm_%j.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/bgb_u29_bsm_%j.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load R-bundle-CRAN/2024.06-foss-2022b
# Set R_LIBS_USER for your installation if needed.
export ORCHID_PAPER2_ROOT=${ORCHID_PROJECT_ROOT}
Rscript orchid_paper2/scripts/107_run_bsm_unified29.R
