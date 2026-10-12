#!/usr/bin/env bash
#SBATCH --job-name=pam_rnaspades
#SBATCH --partition=day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=240G
#SBATCH --output=orchid_paper2/15_ortholog_phylogenomics/logs/pam_rnaspades_%j.out
#SBATCH --error=orchid_paper2/15_ortholog_phylogenomics/logs/pam_rnaspades_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load SPAdes/3.15.5-GCC-12.2.0

PROJECT_ROOT="${1:-${ORCHID_PROJECT_ROOT}}"
WORK_ROOT="${PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics"
TRIM_ROOT="${WORK_ROOT}/transcriptomes/pam/trimmed"
OUTPUT_ROOT="${WORK_ROOT}/transcriptomes/pam/rnaspades_v3"

rnaspades.py \
    --pe-1 1 "${TRIM_ROOT}/SRR7125083_1.clean.fastq.gz" \
    --pe-2 1 "${TRIM_ROOT}/SRR7125083_2.clean.fastq.gz" \
    --pe-1 2 "${TRIM_ROOT}/SRR7125084_1.clean.fastq.gz" \
    --pe-2 2 "${TRIM_ROOT}/SRR7125084_2.clean.fastq.gz" \
    --pe-1 3 "${TRIM_ROOT}/SRR7125085_1.clean.fastq.gz" \
    --pe-2 3 "${TRIM_ROOT}/SRR7125085_2.clean.fastq.gz" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 230 \
    -o "${OUTPUT_ROOT}"

test -s "${OUTPUT_ROOT}/transcripts.fasta"
