#!/usr/bin/env bash
set -euo pipefail

project=${1:-PRJNA1320963}
out=${2:-orchid_paper2/00_metadata/sra_runinfo.csv}
mkdir -p "$(dirname "$out")"

url="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${project}"
curl -L --fail --retry 5 --retry-delay 10 "$url" -o "$out"

if [[ ! -s "$out" ]] || [[ $(wc -l < "$out") -lt 2 ]]; then
  echo "No public SRA runs found for ${project}" >&2
  exit 2
fi

echo "Saved SRA run table to $out"

