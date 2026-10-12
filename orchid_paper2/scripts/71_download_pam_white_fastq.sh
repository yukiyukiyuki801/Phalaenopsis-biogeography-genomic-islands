#!/usr/bin/env bash
#SBATCH --job-name=pam_white_fastq
#SBATCH --partition=day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --output=orchid_paper2/15_ortholog_phylogenomics/logs/pam_fastq_%j.out
#SBATCH --error=orchid_paper2/15_ortholog_phylogenomics/logs/pam_fastq_%j.err

set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"

PROJECT_ROOT="${1:-${ORCHID_PROJECT_ROOT}}"
WORK_ROOT="${PROJECT_ROOT}/orchid_paper2/15_ortholog_phylogenomics"
FASTQ_ROOT="${WORK_ROOT}/downloads/pam_baiyuzan_fastq"
mkdir -p "${FASTQ_ROOT}" "${WORK_ROOT}/logs"

# ENA mirrors provide the submitted paired FASTQ files directly. These paths
# were resolved from the ENA filereport API on 2026-08-05.
for run in SRR7125083 SRR7125084 SRR7125085; do
    suffix="00${run: -1}"
    for mate in 1 2; do
        output="${FASTQ_ROOT}/${run}_${mate}.fastq.gz"
        url="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR712/${suffix}/${run}/${run}_${mate}.fastq.gz"
        curl --fail --location --retry 8 --retry-delay 20 --continue-at - \
            --output "${output}" "${url}"
        gzip -t "${output}"
    done
done

(
    cd "${FASTQ_ROOT}"
    printf '%s  %s\n' aed7f4113e9266f0d090876f6547cca6 SRR7125083_1.fastq.gz
    printf '%s  %s\n' aff3261083c3082b411a724508b5c529 SRR7125083_2.fastq.gz
    printf '%s  %s\n' 9d64779a6ea3a7256d7f4f7f402d851d SRR7125084_1.fastq.gz
    printf '%s  %s\n' 3dcec593966d109aff411244942c6cac SRR7125084_2.fastq.gz
    printf '%s  %s\n' 5eaa05107b4f2aa03db15655d383a652 SRR7125085_1.fastq.gz
    printf '%s  %s\n' eb4e62a97a64ebaff06483e193d53116 SRR7125085_2.fastq.gz
) > "${WORK_ROOT}/metadata/pam_baiyuzan_fastq.expected.md5"

(
    cd "${FASTQ_ROOT}"
    md5sum -c "${WORK_ROOT}/metadata/pam_baiyuzan_fastq.expected.md5"
)

sha256sum "${FASTQ_ROOT}"/*.fastq.gz > \
    "${WORK_ROOT}/metadata/pam_baiyuzan_fastq.sha256"
