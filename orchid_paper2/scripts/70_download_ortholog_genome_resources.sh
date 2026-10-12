#!/usr/bin/env bash
#SBATCH --job-name=orchid_ortholog_dl
#SBATCH --partition=day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=orchid_paper2/15_ortholog_phylogenomics/logs/download_genomes_%j.out
#SBATCH --error=orchid_paper2/15_ortholog_phylogenomics/logs/download_genomes_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"

PROJECT_ROOT="${1:-${ORCHID_PROJECT_ROOT}}"
WORK_ROOT="${PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics"
DATASETS_BIN="${DATASETS_BIN:-/tmp/ncbi_datasets.zip}"
ACCESSIONS="${WORK_ROOT}/metadata/ncbi_annotated_genome_accessions.txt"
ARCHIVE="${WORK_ROOT}/downloads/ncbi_annotated_genomes.zip"
UNPACKED="${WORK_ROOT}/downloads/ncbi_annotated_genomes"

mkdir -p "${WORK_ROOT}/downloads" "${WORK_ROOT}/logs" "${UNPACKED}"

if [[ ! -x "${DATASETS_BIN}" ]]; then
    echo "ERROR: NCBI datasets executable not found: ${DATASETS_BIN}" >&2
    exit 1
fi

"${DATASETS_BIN}" download genome accession \
    --inputfile "${ACCESSIONS}" \
    --include protein,cds,gff3,seq-report \
    --filename "${ARCHIVE}"

unzip -o "${ARCHIVE}" -d "${UNPACKED}"
find "${UNPACKED}" -type f -print0 | sort -z | xargs -0 sha256sum \
    > "${WORK_ROOT}/metadata/ncbi_annotated_genomes.sha256"

