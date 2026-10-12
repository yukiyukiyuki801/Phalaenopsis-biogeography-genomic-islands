#!/usr/bin/env bash
#SBATCH --job-name=mcmctree_hessian
#SBATCH --partition=day,scavenge
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/logs/mcmctree_hessian_%j.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/logs/mcmctree_hessian_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load PAML/4.10.7-GCCcore-12.2.0
ROOT="${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/single_copy_dating"
cd "${ROOT}"
mcmctree config/mcmctree_hessian.ctl
test -s out.BV
cp out.BV in.BV
echo "MCMCTree approximate-likelihood Hessian completed: ${ROOT}/in.BV"
