#!/usr/bin/env bash
#SBATCH --job-name=mcmctree_orchid
#SBATCH --partition=day,scavenge
#SBATCH --array=1-4
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/logs/mcmctree_orchid_%A_%a.out
#SBATCH --error=${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/logs/mcmctree_orchid_%A_%a.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
module load PAML/4.10.7-GCCcore-12.2.0
ROOT="${ORCHID_PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics/single_copy_dating"

case "${SLURM_ARRAY_TASK_ID}" in
  1) scenario=zhang_multipoint; chain=1; seed=20260811; tree=zhang_multipoint.trees ;;
  2) scenario=zhang_multipoint; chain=2; seed=20260812; tree=zhang_multipoint.trees ;;
  3) scenario=givnish_multipoint; chain=1; seed=20260821; tree=givnish_multipoint.trees ;;
  4) scenario=givnish_multipoint; chain=2; seed=20260822; tree=givnish_multipoint.trees ;;
esac

RUN="${ROOT}/dating/${scenario}/chain${chain}"
mkdir -p "${RUN}"
cp "${ROOT}/in.BV" "${RUN}/in.BV"
sed -e "s/SEED_VALUE/${seed}/" -e "s/TREE_FILE/${tree}/" \
  "${ROOT}/config/mcmctree_chain.template.ctl" > "${RUN}/mcmctree.ctl"
cd "${RUN}"
mcmctree mcmctree.ctl
test -s mcmc.txt
test -s out.txt
echo "Completed ${scenario} chain ${chain}"
