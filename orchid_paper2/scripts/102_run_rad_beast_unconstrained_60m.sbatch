#!/usr/bin/env bash
#SBATCH --job-name=rad_uncon60
#SBATCH --partition=ycga,week
#SBATCH --array=1-6
#SBATCH --cpus-per-task=30
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/rad_uncon60_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/rad_uncon60_%A_%a.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
module load Beast/2.7.4-GCC-12.2.0
if (( SLURM_ARRAY_TASK_ID <= 3 )); then
    scenario=zhang
    chain=$SLURM_ARRAY_TASK_ID
else
    scenario=givnish
    chain=$((SLURM_ARRAY_TASK_ID-3))
fi
run="$ROOT/16_rad_beast/runs/${scenario}_combined_rateprior_unconstrained_60m/chain${chain}"
mkdir -p "$run"
cd "$run"
seed=$((202608900 + SLURM_ARRAY_TASK_ID))
state="rad29_peq_${scenario}_combined_rateprior_unconstrained_60m.xml.state"
xml="$ROOT/16_rad_beast/xml/rad29_peq_${scenario}_combined_rateprior_unconstrained_60m.xml"
if [[ -s "$state" ]]; then
    beast -threads "$SLURM_CPUS_PER_TASK" -seed "$seed" -resume -statefile "$state" "$xml"
else
    beast -threads "$SLURM_CPUS_PER_TASK" -seed "$seed" -statefile "$state" "$xml"
fi
