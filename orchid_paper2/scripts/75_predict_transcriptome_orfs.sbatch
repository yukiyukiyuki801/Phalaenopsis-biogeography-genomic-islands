#!/usr/bin/env bash
#SBATCH --job-name=orchid_orfs
#SBATCH --partition=day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=orchid_paper2/15_ortholog_phylogenomics/logs/orfs_%x_%j.out
#SBATCH --error=orchid_paper2/15_ortholog_phylogenomics/logs/orfs_%x_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 TAXON INPUT_TRANSCRIPTS_FASTA[.gz] [PROJECT_ROOT]" >&2
    exit 2
fi

TAXON="$1"
INPUT_FASTA="$2"
PROJECT_ROOT="${3:-${ORCHID_PROJECT_ROOT}}"
WORK_ROOT="${PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics"
TD_ROOT="${WORK_ROOT}/software/TransDecoder-TransDecoder-v5.7.1"
OUTPUT_ROOT="${WORK_ROOT}/transcriptomes/${TAXON}/transdecoder"
TRANSCRIPTS="${OUTPUT_ROOT}/${TAXON}.transcripts.fasta"

mkdir -p "${OUTPUT_ROOT}" "${WORK_ROOT}/logs"

if [[ "${INPUT_FASTA}" == *.gz ]]; then
    gzip -dc "${INPUT_FASTA}" > "${TRANSCRIPTS}"
else
    cp "${INPUT_FASTA}" "${TRANSCRIPTS}"
fi

cd "${OUTPUT_ROOT}"
"${TD_ROOT}/TransDecoder.LongOrfs" -t "${TRANSCRIPTS}" -m 100
"${TD_ROOT}/TransDecoder.Predict" -t "${TRANSCRIPTS}" --single_best_only

test -s "${TRANSCRIPTS}.transdecoder.pep"
test -s "${TRANSCRIPTS}.transdecoder.cds"
cp "${TRANSCRIPTS}.transdecoder.pep" "${WORK_ROOT}/proteomes/${TAXON}.faa"
cp "${TRANSCRIPTS}.transdecoder.cds" "${WORK_ROOT}/proteomes/${TAXON}.cds.fna"

