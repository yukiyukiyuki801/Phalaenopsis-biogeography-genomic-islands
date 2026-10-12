#!/usr/bin/env bash
#SBATCH --job-name=joint_merge
#SBATCH --partition=ycga,day,scavenge,week
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=orchid_paper2/logs/joint_merge_%j.out
#SBATCH --error=orchid_paper2/logs/joint_merge_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load BCFtools/1.21-GCC-13.3.0

project=${ORCHID_PROJECT_ROOT}
cd "$project"
root=orchid_paper2/05_all_sites_vcf
final="$root/orchid_173.callable_all_sites.vcf.gz"
list="$root/shards/vcf/shard_vcfs.list"

# Preserve a successfully completed monolithic result.
if [[ -s "$final" && -s "$final.tbi" ]]; then
  echo "Validated monolithic output already exists; no shard merge required"
  exit 0
fi

: > "$list"
for i in $(seq -w 1 24); do
  f="$root/shards/vcf/shard_0${i}.vcf.gz"
  test -s "$f" && test -s "$f.tbi"
  printf '%s\n' "$f" >> "$list"
done

tmp="$root/orchid_173.callable_all_sites.sharded.tmp.vcf.gz"
bcftools concat --threads "$SLURM_CPUS_PER_TASK" -f "$list" -Oz -o "$tmp"
bcftools index --threads "$SLURM_CPUS_PER_TASK" -t "$tmp"

samples=$(bcftools query -l "$tmp" | wc -l)
[[ "$samples" -eq 173 ]] || { echo "Expected 173 samples, found $samples" >&2; exit 3; }

duplicates=$(bcftools query -f '%CHROM\t%POS\n' "$tmp" | uniq -d | wc -l)
[[ "$duplicates" -eq 0 ]] || { echo "Found $duplicates adjacent duplicate coordinates" >&2; exit 4; }

mv "$tmp" "$final"
mv "$tmp.tbi" "$final.tbi"
bcftools index -n "$final"
echo "Merged 24 shards into $final with $samples samples"
