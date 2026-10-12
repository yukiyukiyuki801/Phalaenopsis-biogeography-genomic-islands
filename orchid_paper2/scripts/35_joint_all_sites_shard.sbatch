#!/usr/bin/env bash
#SBATCH --job-name=joint_shard
#SBATCH --partition=ycga,day,scavenge,week
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --array=1-24
#SBATCH --output=orchid_paper2/logs/joint_shard_%A_%a.out
#SBATCH --error=orchid_paper2/logs/joint_shard_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17

project=${ORCHID_PROJECT_ROOT}
cd "$project"
root=orchid_paper2/05_all_sites_vcf
ref=orchid_paper2/01_reference/GCF_001263595.1_ASM126359v1_genomic.fna
sample_map="$root/gvcf.sample_map.tsv"
tag=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
intervals="$root/shards/intervals/shard_${tag}.bed"
workspace="$root/shards/genomicsdb/shard_${tag}"
out="$root/shards/vcf/shard_${tag}.vcf.gz"
final="$root/orchid_173.callable_all_sites.vcf.gz"

# A successful original monolithic call makes the fallback shards unnecessary.
if [[ -s "$final" && -s "$final.tbi" ]]; then
  echo "Validated monolithic output already exists; shard not required"
  exit 0
fi

test -s "$intervals"
test -s "$sample_map"
mkdir -p "$root/shards/genomicsdb" "$root/shards/vcf"
if [[ -e "$workspace" ]]; then
  echo "Refusing to overwrite existing workspace: $workspace" >&2
  exit 2
fi

gatk --java-options "-Xmx130g" GenomicsDBImport \
  --genomicsdb-workspace-path "$workspace" \
  --sample-name-map "$sample_map" \
  --intervals "$intervals" \
  --merge-input-intervals true \
  --reader-threads "$SLURM_CPUS_PER_TASK" \
  --genomicsdb-shared-posixfs-optimizations true

gatk --java-options "-Xmx130g" GenotypeGVCFs \
  -R "$ref" \
  -V "gendb://$workspace" \
  -L "$intervals" \
  --include-non-variant-sites true \
  -O "$out"

test -s "$out"
test -s "$out.tbi"
