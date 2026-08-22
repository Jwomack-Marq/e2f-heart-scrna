#!/usr/bin/env bash
# Runs the whole 2026-08-21 deliverable: DE -> GO/GSEA/figures -> Excel.
#
# Everything happens inside e2f-enrich:latest because this machine has neither R
# nor a Python scientific stack. The bundle is mounted READ-ONLY: unlike the
# shiny_app/build_*.R scripts, nothing here rewrites app_data.rds.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO/deliverables/2026-08-21"
IMG="e2f-enrich:latest"
DIR="$REPO/analysis/2026-08-21_email"

mkdir -p "$OUT/csv" "$OUT/plots/part1" "$OUT/plots/part2"

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "== building $IMG (first run only; the Bioconductor layer takes a while) =="
  docker build -f "$DIR/Dockerfile" -t "$IMG" "$DIR"
fi

run() {
  docker run --rm \
    -v "$REPO/shiny_app/app_data.rds:/in/app_data.rds:ro" \
    -v "$REPO:/repo:ro" \
    -v "$OUT:/out" \
    -u "$(id -u):$(id -g)" -e HOME=/tmp \
    "$IMG" Rscript "/repo/analysis/2026-08-21_email/$1"
}

for step in 01_de.R 02_enrich.R 05_mito_sensitivity.R 03_excel.R 04_cover.R 99_verify.R; do
  echo; echo "===================== $step ====================="
  run "$step"
done

echo; echo "== outputs =="
ls -lh "$OUT"/*.xlsx
echo "figures: $(find "$OUT/plots" -type f | wc -l)"
echo "csv:     $(find "$OUT/csv"   -type f | wc -l)"
