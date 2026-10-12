#!/usr/bin/env bash
#SBATCH --job-name=orthofinder10_fast
#SBATCH --partition=day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=180G
#SBATCH --output=orchid_paper2/15_ortholog_phylogenomics/logs/orthofinder10_fast_%j.out
#SBATCH --error=orchid_paper2/15_ortholog_phylogenomics/logs/orthofinder10_fast_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module purge
module load OrthoFinder/2.5.4-foss-2020b
module load DIAMOND/2.0.15-GCCcore-10.2.0
module load MCL/14.137-GCCcore-10.2.0

PROJECT_ROOT="${1:-${ORCHID_PROJECT_ROOT}}"
WORK_ROOT="${PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics"
INPUT_DIR="${WORK_ROOT}/orthofinder_input_10taxa"

for executable in orthofinder diamond mcl; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "ERROR: required executable is unavailable: ${executable}" >&2
        exit 1
    fi
    echo "DEPENDENCY_OK ${executable}=$(command -v "${executable}")"
done

expected=(PAM PAP PEQ APO DCAT DNOB VPLA PGUA ASPA RICE)
for taxon in "${expected[@]}"; do
    file="${INPUT_DIR}/${taxon}.faa"
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: missing or empty OrthoFinder input: ${file}" >&2
        exit 1
    fi
done

# DendroBLAST is OrthoFinder's fast default gene-tree workflow. Unlike the
# prior -M msa/-T iqtree run, it does not infer a full ML tree for every one
# of the >23,000 orthogroups. IQ-TREE will be reserved for selected strict
# single-copy orthologues downstream.
orthofinder -f "${INPUT_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -a "${SLURM_CPUS_PER_TASK}" \
    -S diamond \
    -M dendroblast \
    -n 10taxa_dendroblast
