#!/usr/bin/env bash
#SBATCH --job-name=ba3_tune
#SBATCH --partition=ycga,day,scavenge
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=orchid_paper2/logs/ba3_autotune_%j.out
#SBATCH --error=orchid_paper2/logs/ba3_autotune_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"

project_dir=${ORCHID_PROJECT_ROOT}
run_dir="$project_dir/orchid_paper2/18_contemporary_migration/runs/all_ld_6596/autotune"
input_source="$project_dir/orchid_paper2/18_contemporary_migration/inputs/orchid_18pop_all_ld_pruned.immanc"
binary_dir="$project_dir/orchid_paper2/software/bin"
autotune="$project_dir/orchid_paper2/software/BA3-SNPS-autotune/BA3-SNPS-autotune.py"

mkdir -p "$run_dir" "$project_dir/orchid_paper2/logs"
cp "$input_source" "$run_dir/orchid_18pop_all_ld_pruned.immanc"
cd "$run_dir"
export PATH="$binary_dir:$PATH"

python "$autotune" \
  -i orchid_18pop_all_ld_pruned.immanc \
  -l 6596 \
  -g 100000 \
  -b 10000 \
  -r 10 \
  -s 74001 \
  -o autotune.out.txt \
  > autotune.log 2>&1

test -s orchid_18pop_all_ld_pruned.immanc.finalParams.txt
