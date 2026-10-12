#!/usr/bin/env bash
#SBATCH --job-name=rad_uncon_final_mcc
#SBATCH --partition=day,ycga,scavenge
#SBATCH --array=1-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/rad_uncon_final_mcc_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/rad_uncon_final_mcc_%A_%a.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load Beast/2.7.4-GCC-12.2.0
export _JAVA_OPTIONS=-Xmx80g
base=${ORCHID_PROJECT_ROOT}/orchid_paper2/16_rad_beast/final_unconstrained_60m_summary
if [[ $SLURM_ARRAY_TASK_ID == 1 ]]; then s=zhang; else s=givnish; fi
treeannotator -lowMem -heights median -burnin 0 -limit 0.0 "$base/${s}_unconstrained_merged_postburn.trees" "$base/${s}_unconstrained_final_mcc.tree"
