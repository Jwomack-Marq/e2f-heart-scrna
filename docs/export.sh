#!/usr/bin/env bash
# Regenerate every figure and generated table the methods book embeds.
#
# Runs ONCE, offline -- it is not part of rendering the book. Its outputs
# (docs/assets/*.svg|png, docs/_generated/*.md) are committed, so `docs/render.sh`
# works on a machine with no R and no data bundle.
#
# The bundle is mounted READ-ONLY. Unlike shiny_app/build_*.R, nothing here
# rewrites app_data.rds -- the book documents the analysis, it does not redo it.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMG="e2f-export:latest"
BUNDLE="$REPO/shiny_app/app_data.rds"

if [[ ! -f "$BUNDLE" ]]; then
  echo "error: $BUNDLE not found." >&2
  echo "       The bundle is a git-ignored build artifact; see README.md 'Data pipeline'." >&2
  exit 1
fi

if ! docker image inspect e2f-enrich:latest >/dev/null 2>&1; then
  echo "== building e2f-enrich:latest (first run only; the Bioconductor layer takes a while) =="
  docker build -f "$REPO/analysis/2026-08-21_email/Dockerfile" -t e2f-enrich:latest \
    "$REPO/analysis/2026-08-21_email"
fi
if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "== building $IMG (adds SeuratObject, for reading the original run's objects) =="
  docker build -f "$REPO/docs/Dockerfile.export" -t "$IMG" "$REPO/docs"
fi

mkdir -p "$REPO/docs/assets" "$REPO/docs/_generated" "$REPO/slides/assets"

# --user so the SVGs are not root-owned; HOME=/tmp because that leaves $HOME
# unwritable (same reason analysis/2026-08-21_email/run.sh does it).
docker run --rm \
  -v "$BUNDLE:/in/app_data.rds:ro" \
  -v "$REPO:/repo:ro" \
  -v "$REPO/original_Han_analysis:/orig:ro" \
  -v "$REPO/docs/assets:/out/assets" \
  -v "$REPO/docs/_generated:/out/generated" \
  -v "$REPO/slides/assets:/out/slides" \
  -e SLIDE_PNG_DIR=/out/slides \
  -u "$(id -u):$(id -g)" -e HOME=/tmp \
  "$IMG" Rscript /repo/docs/export_assets.R "$@"

echo
echo "== outputs =="
echo "figures: $(find "$REPO/docs/assets" -type f | wc -l)  ($(du -sh "$REPO/docs/assets" | cut -f1))"
echo "tables:  $(find "$REPO/docs/_generated" -type f | wc -l)"
