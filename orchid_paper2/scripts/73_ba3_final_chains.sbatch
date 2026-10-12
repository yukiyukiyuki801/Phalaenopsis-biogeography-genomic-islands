#!/usr/bin/env bash
#SBATCH --job-name=ba3_final
#SBATCH --partition=ycga,day,scavenge
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-5
#SBATCH --output=orchid_paper2/logs/ba3_final_%A_%a.out
#SBATCH --error=orchid_paper2/logs/ba3_final_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"

project_dir=${ORCHID_PROJECT_ROOT}
base="$project_dir/orchid_paper2/18_contemporary_migration/runs/all_ld_6596"
input_source="$project_dir/orchid_paper2/18_contemporary_migration/inputs/orchid_18pop_all_ld_pruned.immanc"
params="$base/autotune/orchid_18pop_all_ld_pruned.immanc.finalParams.txt"
binary="$project_dir/orchid_paper2/software/bin/BA3-SNPS"
chain=${SLURM_ARRAY_TASK_ID}
run_dir="$base/chain_${chain}"

test -s "$params"
read -r delta_m delta_a delta_f < <(awk '!/^#/ && NF>=3 {print $1,$2,$3; exit}' "$params")
test -n "$delta_m" && test -n "$delta_a" && test -n "$delta_f"

mkdir -p "$run_dir"
cp "$input_source" "$run_dir/orchid_18pop_all_ld_pruned.immanc"
cd "$run_dir"
seed=$((75000 + chain))

"$binary" \
  -F orchid_18pop_all_ld_pruned.immanc \
  -l 6596 \
  -i 10000000 \
  -b 1000000 \
  -n 1000 \
  -s "$seed" \
  -m "$delta_m" \
  -a "$delta_a" \
  -f "$delta_f" \
  -t -u -g \
  -o "chain_${chain}.out.txt" \
  > "chain_${chain}.stdout.txt" \
  2> "chain_${chain}.stderr.txt"

test -s "chain_${chain}.out.txt"
test -s orchid_18pop_all_ld_pruned.trace.txt
