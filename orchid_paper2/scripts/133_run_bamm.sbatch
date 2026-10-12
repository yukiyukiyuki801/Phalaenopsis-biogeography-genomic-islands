#!/bin/bash
#SBATCH --job-name=orchid_bamm
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=1-24
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/bamm_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/bamm_%A_%a.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load GCCcore/12.2.0
cd ${ORCHID_PROJECT_ROOT}
run_dir=$(awk -F '\t' -v idx="$SLURM_ARRAY_TASK_ID" 'NR==idx+1 {print $4}' orchid_paper2/22_diversification/inputs/bamm_run_registry.tsv)
test -n "$run_dir"
cd "$run_dir"
grep -Eq '^numberOfGenerations[[:space:]]*=[[:space:]]*10000000[[:space:]]*$' control.txt
${ORCHID_PROJECT_ROOT}/orchid_paper2/22_diversification/software/bamm/bamm -c control.txt
