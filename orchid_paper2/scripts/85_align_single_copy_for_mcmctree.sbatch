#!/usr/bin/env bash
#SBATCH --job-name=sco_align
#SBATCH --partition=day,scavenge
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=90G
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/logs/sco_align_%j.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/logs/sco_align_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load R-bundle-CRAN/2024.06-foss-2022b
ROOT="${ORCHID_PROJECT_ROOT}/orchid_paper2"
OF="${ROOT}/15_ortholog_phylogenomics/orthofinder_input_10taxa/OrthoFinder/Results_10taxa_dendroblast"
ENV="${ROOT}/15_ortholog_phylogenomics/envs/fasttree"
WORK="${ROOT}/15_ortholog_phylogenomics/single_copy_dating"
SRC="${OF}/Single_Copy_Orthologue_Sequences"
ALN="${WORK}/alignments"
TRIM="${WORK}/trimmed"
MATRIX="${WORK}/matrix"
mkdir -p "${ALN}" "${TRIM}" "${MATRIX}"

test -x "${ENV}/bin/mafft"
test -x "${ENV}/bin/trimal"
test "$(find "${SRC}" -maxdepth 1 -type f -name '*.fa' | wc -l)" -eq 353

export SRC ALN TRIM ENV
find "${SRC}" -maxdepth 1 -type f -name '*.fa' -print0 | \
  xargs -0 -n 1 -P "${SLURM_CPUS_PER_TASK}" bash -c '
    f="$1"; b=$(basename "$f" .fa)
    "${ENV}/bin/mafft" --auto --thread 1 --quiet "$f" > "${ALN}/${b}.aln.fa"
    "${ENV}/bin/trimal" -in "${ALN}/${b}.aln.fa" -out "${TRIM}/${b}.trim.fa" -automated1
  ' _

Rscript "${ROOT}/scripts/84_concatenate_single_copy_proteins.R" "${TRIM}" "${MATRIX}"
echo "Completed strict single-copy protein matrix"
cat "${MATRIX}/strict_single_copy_concat_summary.tsv"
