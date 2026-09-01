#!/usr/bin/env bash
# Render the meeting deck to PowerPoint (and a PDF preview if LibreOffice is around).
#
# Needs only the Quarto CLI: the deck is `engine: markdown`, so nothing executes at
# render time. Figures come from slides/assets/, which docs/export.sh fills with raster
# twins of the methods-book figures (pandoc's pptx writer cannot be trusted with SVG).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMG="e2f-docs:latest"
SRC="${1:-qc-and-processing-review.qmd}"

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "== building $IMG =="
  docker build -f "$REPO/docs/Dockerfile" -t "$IMG" "$REPO/docs"
fi

# fail early on a figure the deck references but nobody exported
missing=0
while read -r a; do
  [[ -f "$REPO/slides/$a" ]] || { echo "missing asset: $a" >&2; missing=1; }
done < <(grep -o 'assets/[A-Za-z0-9._-]*' "$REPO/slides/$SRC" | sort -u)
if [[ $missing -eq 1 ]]; then
  echo "run docs/export.sh to regenerate slides/assets/" >&2
  exit 1
fi

docker run --rm -v "$REPO/slides:/slides" -u "$(id -u):$(id -g)" -e HOME=/tmp \
  "$IMG" quarto render "/slides/$SRC" --to pptx

OUT="$REPO/slides/${SRC%.qmd}.pptx"
echo "Wrote $OUT  ($(du -h "$OUT" | cut -f1))"

# Optional PDF, purely for a quick read-through without PowerPoint.
if command -v libreoffice >/dev/null 2>&1; then
  libreoffice --headless --convert-to pdf --outdir "$REPO/slides" "$OUT" >/dev/null 2>&1 || true
  [[ -f "${OUT%.pptx}.pdf" ]] && echo "Wrote ${OUT%.pptx}.pdf (preview only)"
fi
