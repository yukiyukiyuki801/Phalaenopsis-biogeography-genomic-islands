#!/usr/bin/env bash
#SBATCH --job-name=rad_concat29
#SBATCH --partition=day,scavenge
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/rad_concat29_%j.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/logs/rad_concat29_%j.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
python "$ROOT/scripts/86_build_rad_population_concat.py" \
  --vcf "$ROOT/05_all_sites_vcf/filtered/orchid_173.paper1_like.biallelic_snps_invariant.vcf.gz" \
  --bed "$ROOT/05_all_sites_vcf/callable/cohort_callable.bed" \
  --popmap "$ROOT/14_snapper/inputs/whole29/sample_to_population.tsv" \
  --outdir "$ROOT/16_rad_beast/inputs"
