#!/usr/bin/env bash
#SBATCH --job-name=orchid_map
#SBATCH --partition=day
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-173%12
#SBATCH --output=orchid_paper2/logs/map_%A_%a.out
#SBATCH --error=orchid_paper2/logs/map_%A_%a.err

set -euo pipefail
module load BWA/0.7.17-GCCcore-12.2.0
module load SAMtools/1.21-GCC-12.2.0
module load Trimmomatic/0.39-Java-11

manifest=${MANIFEST:-orchid_paper2/00_metadata/run_manifest.tsv}
line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$manifest")
run=$(awk '{print $1}' <<< "$line")
sample=$(awk '{print $2}' <<< "$line")
library=$(awk '{print $3}' <<< "$line")

raw=orchid_paper2/02_sra/fastq
clean=orchid_paper2/03_filtered_reads
bamdir=orchid_paper2/04_bam
tmp=orchid_paper2/tmp/map_${sample}
ref=orchid_paper2/01_reference/GCF_001263595.1_ASM126359v1_genomic.fna
mkdir -p "$clean" "$bamdir/qc" "$tmp"

r1="$raw/${run}_1.fastq.gz"
r2="$raw/${run}_2.fastq.gz"
p1="$clean/${sample}_R1.fastq.gz"
u1="$clean/${sample}_R1.unpaired.fastq.gz"
p2="$clean/${sample}_R2.fastq.gz"
u2="$clean/${sample}_R2.unpaired.fastq.gz"

# SRA FASTQs are already demultiplexed by sample; do not run process_radtags.
# Paper 1 only records trimming bases below Phred 15. These explicit settings
# reconstruct that intent; they are not claimed to be the lost original command.
java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" PE \
  -threads "$SLURM_CPUS_PER_TASK" -phred33 \
  "$r1" "$r2" "$p1" "$u1" "$p2" "$u2" \
  SLIDINGWINDOW:4:15 MINLEN:36

rg="@RG\\tID:${library}.${sample}\\tSM:${sample}\\tLB:${library}\\tPL:ILLUMINA"
bwa mem -t "$SLURM_CPUS_PER_TASK" -R "$rg" "$ref" "$p1" "$p2" \
  | samtools sort -@ 4 -m 4G -T "$tmp/${sample}" \
      -o "$bamdir/${sample}.sorted.bam" -

samtools index -@ 4 "$bamdir/${sample}.sorted.bam"
samtools quickcheck -v "$bamdir/${sample}.sorted.bam"
samtools flagstat -@ 4 "$bamdir/${sample}.sorted.bam" \
  > "$bamdir/qc/${sample}.flagstat.txt"
samtools stats -@ 4 "$bamdir/${sample}.sorted.bam" \
  > "$bamdir/qc/${sample}.stats.txt"
