#!/usr/bin/env bash
#SBATCH --job-name=orchid_ref
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=orchid_paper2/logs/reference_%j.out
#SBATCH --error=orchid_paper2/logs/reference_%j.err

set -euo pipefail
module load BWA/0.7.17-GCCcore-12.2.0
module load SAMtools/1.21-GCC-13.3.0
module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17

root=orchid_paper2/01_reference
gz="$root/GCF_001263595.1_ASM126359v1_genomic.fna.gz"
ref="$root/GCF_001263595.1_ASM126359v1_genomic.fna"

gzip -t "$gz"
if [[ ! -s "$ref" ]]; then
  gzip -dc "$gz" > "$ref"
fi

bwa index "$ref"
samtools faidx "$ref"
gatk CreateSequenceDictionary -R "$ref" -O "${ref%.fna}.dict"
sha256sum "$gz" "$ref" > "$root/reference.sha256"

