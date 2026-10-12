#!/usr/bin/env bash
#SBATCH --job-name=pam_fastp
#SBATCH --partition=day
#SBATCH --time=23:59:00
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=orchid_paper2/15_ortholog_phylogenomics/logs/pam_fastp_%A_%a.out
#SBATCH --error=orchid_paper2/15_ortholog_phylogenomics/logs/pam_fastp_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load fastp/0.23.2-GCCcore-10.2.0

PROJECT_ROOT="${1:-${ORCHID_PROJECT_ROOT}}"
WORK_ROOT="${PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics"
INPUT_ROOT="${WORK_ROOT}/downloads/pam_baiyuzan_fastq"
OUTPUT_ROOT="${WORK_ROOT}/transcriptomes/pam/trimmed"
mkdir -p "${OUTPUT_ROOT}" "${WORK_ROOT}/logs"

runs=(SRR7125083 SRR7125084 SRR7125085)
run="${runs[${SLURM_ARRAY_TASK_ID}]}"

fastp \
    --in1 "${INPUT_ROOT}/${run}_1.fastq.gz" \
    --in2 "${INPUT_ROOT}/${run}_2.fastq.gz" \
    --out1 "${OUTPUT_ROOT}/${run}_1.clean.fastq.gz" \
    --out2 "${OUTPUT_ROOT}/${run}_2.clean.fastq.gz" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --qualified_quality_phred 20 --unqualified_percent_limit 30 \
    --length_required 50 \
    --thread "${SLURM_CPUS_PER_TASK}" \
    --html "${OUTPUT_ROOT}/${run}.fastp.html" \
    --json "${OUTPUT_ROOT}/${run}.fastp.json"

gzip -t "${OUTPUT_ROOT}/${run}_1.clean.fastq.gz"
gzip -t "${OUTPUT_ROOT}/${run}_2.clean.fastq.gz"

