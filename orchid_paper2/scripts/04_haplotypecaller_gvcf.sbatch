#!/usr/bin/env bash
#SBATCH --job-name=orchid_gvcf
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-173%12
#SBATCH --output=orchid_paper2/logs/gvcf_%A_%a.out
#SBATCH --error=orchid_paper2/logs/gvcf_%A_%a.err

set -euo pipefail
module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17

manifest=${MANIFEST:-orchid_paper2/00_metadata/run_manifest.tsv}
line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$manifest")
sample=$(awk '{print $2}' <<< "$line")
ref=orchid_paper2/01_reference/GCF_001263595.1_ASM126359v1_genomic.fna
bam=orchid_paper2/04_bam/${sample}.sorted.bam
outdir=orchid_paper2/05_all_sites_vcf/gvcf
mkdir -p "$outdir"

gatk --java-options "-Xmx20g" HaplotypeCaller \
  -R "$ref" -I "$bam" -O "$outdir/${sample}.g.vcf.gz" \
  -ERC GVCF --native-pair-hmm-threads "$SLURM_CPUS_PER_TASK"

# HaplotypeCaller creates the Tabix index with the bgzipped gVCF. A second
# IndexFeatureFile pass is redundant and can exhaust its small default heap on
# reference assemblies containing many contigs.
test -s "$outdir/${sample}.g.vcf.gz.tbi"
