#!/usr/bin/env bash
#SBATCH --job-name=orchid_sra
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-173%12
#SBATCH --output=orchid_paper2/logs/sra_%A_%a.out
#SBATCH --error=orchid_paper2/logs/sra_%A_%a.err

set -euo pipefail
module load SRA-Toolkit/3.2.0-gompi-2024a

manifest=${MANIFEST:-orchid_paper2/00_metadata/run_manifest.tsv}
line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$manifest")
run=$(awk '{print $1}' <<< "$line")
sample=$(awk '{print $2}' <<< "$line")

if [[ -z "$run" || -z "$sample" ]]; then
  echo "Missing run/sample for task $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

archive=orchid_paper2/02_sra/archive
fastq=orchid_paper2/02_sra/fastq
tmp=orchid_paper2/tmp/fasterq_${run}
mkdir -p "$archive" "$fastq" "$tmp"

prefetch --type sra --max-size u -O "$archive" "$run"
if [[ ! -s "$archive/$run/$run.sra" ]]; then
  echo "Full-quality SRA archive unavailable after prefetch: $run" >&2
  echo "A Lite/no-quality object is not accepted for Phred-based trimming." >&2
  exit 4
fi
fasterq-dump --split-files --threads "$SLURM_CPUS_PER_TASK" \
  --temp "$tmp" --outdir "$fastq" "$archive/$run/$run.sra"

if [[ -s "$fastq/${run}_1.fastq" && -s "$fastq/${run}_2.fastq" ]]; then
  gzip -f "$fastq/${run}_1.fastq" "$fastq/${run}_2.fastq"
else
  echo "Expected paired FASTQs were not produced for $run" >&2
  exit 3
fi

sha256sum "$fastq/${run}_1.fastq.gz" "$fastq/${run}_2.fastq.gz" \
  > "$fastq/${run}.sha256"
