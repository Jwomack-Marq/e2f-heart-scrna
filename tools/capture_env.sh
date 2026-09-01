#!/usr/bin/env bash
# Record the package versions inside each Docker image that produced results.
#
# WHY THIS EXISTS AND WHAT IT IS NOT. There is no renv.lock here and an renv.lock
# would be the wrong tool: no R is installed on the host, the Docker images ARE the
# environment, and each of them installs from a floating repo (install2.r /
# BiocManager with no version constraints). So `docker build` today and `docker build`
# in six months produce different package versions from the same Dockerfile.
#
# This does not fix that. It RECORDS it: for each built image, the R version, the
# platform, the image id, and every installed package version. That is enough to say
# afterwards which versions a given result was computed with -- which is the thing we
# could not say at all before -- and enough to notice when a rebuild has moved.
#
# It reads the images that exist locally. It cannot recover versions for an image that
# was never built or has been pruned, and it does NOT rebuild anything.
#
#   bash tools/capture_env.sh            # every known image
#   bash tools/capture_env.sh e2f-enrich # just one
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO/env"
mkdir -p "$OUT"

# image -> what it is used for. Order is roughly the dependency order of the builds.
declare -A ROLE=(
  [e2f-enrich:latest]="enrichment (clusterProfiler/fgsea/msigdbr); base for the two below"
  [e2f-seurat-full:latest]="Seurat pipeline: SCTransform, Harmony, DESeq2 -- our_analysis/"
  [e2f-export:latest]="docs/export.sh -- regenerates the methods-book figures"
  [lab-server-e2f-heart-scrna-dev:latest]="the deployed Shiny app"
)
WANT=("${@:-}")

for img in $(printf '%s\n' "${!ROLE[@]}" | sort); do
  if [[ -n "${WANT[0]:-}" ]] && [[ "$img" != *"${WANT[0]}"* ]]; then continue; fi
  if ! docker image inspect "$img" >/dev/null 2>&1; then
    echo "skip (not built): $img" >&2; continue
  fi
  f="$OUT/$(echo "$img" | tr ':/' '__').tsv"
  echo "== $img -> env/$(basename "$f")"
  {
    echo "# image	$img"
    echo "# role	${ROLE[$img]}"
    echo "# image_id	$(docker image inspect --format '{{.Id}}' "$img")"
    echo "# image_created	$(docker image inspect --format '{{.Created}}' "$img")"
    echo "# captured	$(date -u +%Y-%m-%dT%H:%M:%SZ) by tools/capture_env.sh"
  } > "$f"
  # installed.packages() returns one row PER LIBRARY PATH, so a package installed in
  # two libs appears twice with two versions and the file cannot tell you which one
  # library() actually loads. These images all have that (Matrix ships with R and is
  # then upgraded into the site library). Report the version that WINS -- the first
  # match walking .libPaths() in order, which is what library() resolves to -- and
  # record the shadowed copy in a third column rather than dropping it silently.
  docker run --rm --network none "$img" R -q --no-echo -e '
    cat(sprintf("# r_version\t%s\n", R.version.string))
    cat(sprintf("# platform\t%s\n", R.version$platform))
    cat(sprintf("# libpaths\t%s\n", paste(.libPaths(), collapse = " : ")))
    ip <- as.data.frame(installed.packages()[, c("Package","Version","LibPath")],
                        stringsAsFactors = FALSE)
    ip$rank <- match(ip$LibPath, .libPaths())
    ip <- ip[order(ip$Package, ip$rank), ]
    eff <- ip[!duplicated(ip$Package), ]
    sh  <- ip[duplicated(ip$Package), ]
    eff$shadowed <- vapply(eff$Package, function(p) {
      v <- sh$Version[sh$Package == p]
      if (length(v)) paste(v, collapse = ",") else ""
    }, character(1))
    eff <- eff[order(tolower(eff$Package)), ]
    cat(sprintf("# n_packages\t%d\n", nrow(eff)))
    cat(sprintf("# n_shadowed\t%d\n", sum(nzchar(eff$shadowed))))
    cat("package\tversion\tshadowed_versions\n")
    cat(sprintf("%s\t%s\t%s\n", eff$Package, eff$Version, eff$shadowed), sep = "")' >> "$f"
  echo "   $(grep -c $'^[^#]' "$f" | awk '{print $1-1}') packages, R $(grep '^# r_version' "$f" | cut -f2)"
done
echo
echo "Wrote $OUT/. Commit these -- they are the record of what produced the results."
