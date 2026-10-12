#!/usr/bin/env bash
set -euo pipefail
: "${ORCHID_PROJECT_ROOT:?Set ORCHID_PROJECT_ROOT to the repository root}"
ROOT=${ORCHID_PROJECT_ROOT}/orchid_paper2
OUT="$ROOT/17_herodotools/inputs/gadm41"
mkdir -p "$OUT"
for code in PHL TWN IDN MYS BRN TLS PNG AUS; do
  file="$OUT/gadm41_${code}.gpkg"
  [[ -s "$file" ]] || curl -fL --retry 3 \
    "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_${code}.gpkg" \
    -o "$file"
done
sha256sum "$OUT"/*.gpkg > "$OUT/sha256sums.txt"
