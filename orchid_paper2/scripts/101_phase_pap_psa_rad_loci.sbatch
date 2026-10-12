#!/usr/bin/env bash
#SBATCH --job-name=pap_psa_phase
#SBATCH --partition=ycga,day,scavenge
#SBATCH --array=1-650
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=orchid_paper2/19_local_genealogy/logs/phase_%A_%a.out
#SBATCH --error=orchid_paper2/19_local_genealogy/logs/phase_%A_%a.err
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"

project=${ORCHID_PROJECT_ROOT}
cd "$project"
module load BCFtools/1.21-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0

base=orchid_paper2/19_local_genealogy
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$base/inputs/pap_psa_nuclear_loci.tsv")
locus=$(cut -f1 <<< "$line")
chrom=$(cut -f2 <<< "$line")
start=$(cut -f3 <<< "$line")
end=$(cut -f4 <<< "$line")
env=orchid_paper2/envs/local_genealogy-pip/bin
outdir="$base/phased_loci"
mkdir -p "$outdir"

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/pap_psa_phase.${locus}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

bcftools view -r "${chrom}:${start}-${end}" -Oz \
  -o "$tmpdir/${locus}.vcf.gz" "$base/inputs/pap_psa_rad_variants.vcf.gz"
bcftools index -f -t "$tmpdir/${locus}.vcf.gz"

mapfile -t bams < <(tail -n +2 "$base/inputs/pap_psa_sample_bams.tsv" | cut -f3)
"$env/whatshap" phase --no-reference --tag PS \
  -o "$outdir/${locus}.phased.vcf.gz" "$tmpdir/${locus}.vcf.gz" "${bams[@]}"
bcftools index -f -t "$outdir/${locus}.phased.vcf.gz"
