#!/usr/bin/env bash
# Render the methods book to docs/_site/.
#
# Needs only Docker: the image carries the Quarto CLI and nothing else, because
# no chapter executes code (see docs/Dockerfile). Figures must already be in
# docs/assets/ -- regenerate them with docs/export.sh if the bundle changed.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMG="e2f-docs:latest"

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "== building $IMG (first run only) =="
  docker build -f "$REPO/docs/Dockerfile" -t "$IMG" "$REPO/docs"
fi

echo "== checking chapter/asset coverage =="
python3 "$REPO/tools/check_docs_coverage.py"

echo
echo "== quarto render =="
docker run --rm \
  -v "$REPO/docs:/docs" \
  -u "$(id -u):$(id -g)" -e HOME=/tmp \
  "$IMG" quarto render --output-dir _site

echo
echo "Wrote $REPO/docs/_site/index.html"
echo "Open it directly, or serve the folder:  python3 -m http.server -d docs/_site 8080"
