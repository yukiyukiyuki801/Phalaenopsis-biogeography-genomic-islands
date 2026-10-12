#!/usr/bin/env bash
#SBATCH --job-name=bgb_u29
#SBATCH --partition=day,scavenge
#SBATCH --array=1-606
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/bgb_u29_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/bgb_u29_%A_%a.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load R-bundle-CRAN/2024.06-foss-2022b
# Set R_LIBS_USER for your installation if needed.
export ORCHID_PAPER2_ROOT=${ORCHID_PROJECT_ROOT}
models=(DEC DECJ DIVALIKE DIVALIKEJ BAYAREALIKE BAYAREALIKEJ)
idx=$((SLURM_ARRAY_TASK_ID-1))
tree_num=$((idx/6))
model=${models[$((idx%6))]}
if (( tree_num == 0 )); then tree_id=mcc; else tree_id=$tree_num; fi
Rscript orchid_paper2/scripts/104_run_biogeobears_unified29.R "$tree_id" "$model"
