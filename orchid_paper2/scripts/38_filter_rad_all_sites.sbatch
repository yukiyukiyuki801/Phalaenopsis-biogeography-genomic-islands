#!/usr/bin/env bash
#SBATCH --job-name=rad_filter
#SBATCH --partition=ycga,day,scavenge,week
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=orchid_paper2/logs/rad_filter_%j.out
#SBATCH --error=orchid_paper2/logs/rad_filter_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load BCFtools/1.21-GCC-13.3.0

project=${ORCHID_PROJECT_ROOT}
cd "$project"
root=orchid_paper2/05_all_sites_vcf
input="$root/orchid_173.callable_all_sites.vcf.gz"
outdir="$root/filtered"
groups="$outdir/sample_groups.tsv"
mkdir -p "$outdir"

test -s "$input" && test -s "$input.tbi"

# The three released taxa are unambiguously encoded in the sample prefix.
bcftools query -l "$input" | awk '
  /^PAM/{print $0 "\tPAM"; next}
  /^PAP/{print $0 "\tPAP"; next}
  /^PSA/{print $0 "\tPSA"; next}
  {print "Unexpected sample prefix: " $0 > "/dev/stderr"; exit 2}
' > "$groups"

run_profile() {
  profile=$1
  mask=$2
  quality=$3
  out="$outdir/orchid_173.${profile}.biallelic_snps_invariant.vcf.gz"

  if [[ "$mask" == none ]]; then
    bcftools view -Ou -e 'N_ALT>1 || (N_ALT=1 && TYPE!="snp")' "$input"
  else
    bcftools view -Ou -e 'N_ALT>1 || (N_ALT=1 && TYPE!="snp")' "$input" \
      | bcftools +setGT -Ou -- -t q -n . -i "$mask"
  fi \
    | bcftools +fill-tags -Ou -- -S "$groups" \
        -t AC,AN,AF,NS,F_MISSING \
    | bcftools view -Oz --threads "$SLURM_CPUS_PER_TASK" \
        -e "$quality" -o "$out"

  bcftools index --threads "$SLURM_CPUS_PER_TASK" -t "$out"
  bcftools stats -s - "$out" > "$outdir/${profile}.bcftools.stats.txt"
  bcftools view -Oz --threads "$SLURM_CPUS_PER_TASK" \
    -i 'F_MISSING<=0.70' "$out" \
    -o "$outdir/orchid_173.${profile}.global_call30.vcf.gz"
  bcftools index --threads "$SLURM_CPUS_PER_TASK" -t \
    "$outdir/orchid_173.${profile}.global_call30.vcf.gz"
}

# Paper 1 sensitivity: retain existing GATK calls and reproduce its permissive
# 70%-missing-site ceiling. Pair-specific filters are applied downstream.
run_profile paper1_like none 'INFO/DP<30'

