# Environment record

What package versions the results in this repo were actually computed with.

Regenerate with `bash tools/capture_env.sh` (reads the built Docker images; rebuilds
nothing). One `.tsv` per image: a `#`-prefixed header with the R version, platform,
image id and capture time, then `package version shadowed_versions`.

## What this is, and what it is not

It is a **record**, not a lockfile. There is no `renv.lock` and one would be the wrong
tool: no R is installed on the host, the Docker images *are* the environment, and every
Dockerfile installs from a floating repo — `install2.r` and `BiocManager::install()`
with no version constraints anywhere. So the same Dockerfile built today and built in
six months gives different package versions.

This does not fix that. It makes the drift **visible and attributable**: afterwards you
can say which versions produced a given figure or table, and you can diff two captures
to see what a rebuild moved. Pinning the repos to a dated Posit Package Manager snapshot
is the actual fix, and it is deliberately not done here — it changes what a rebuild
installs, so it needs to happen at a moment when the images can be rebuilt and the
results re-verified, not as a side effect of writing them down.

`our_analysis/Dockerfile.seurat` already says this in its own header, and it is the
reason `processing/seurat.cm.subclustered.rds` cannot be compared across rebuilds.

## Images

| Image | R | Packages | Role |
|---|---|---|---|
| `e2f-enrich` | 4.5.1 | 180 | enrichment (clusterProfiler / fgsea / msigdbr); base for the two below |
| `e2f-seurat-full` | 4.5.1 | 271 | Seurat pipeline — SCTransform, Harmony, DESeq2 (`our_analysis/`) |
| `e2f-export` | 4.5.1 | 191 | `docs/export.sh`, the methods-book figures |
| `lab-server-e2f-heart-scrna-dev` | 4.5.1 | 110 | the deployed Shiny app |

`e2f-docs` carries only the Quarto CLI — no R, deliberately, so a stray executable chunk
in a chapter fails loudly instead of running.

## Two things the first capture surfaced

**The deployed app is not on the same versions as anything tested locally.** It has
bslib 0.12.0, plotly 4.12.1 and Matrix 1.7-6, against 0.11.0 / 4.12.0 / 1.7-5 in the
three analysis images. bslib is the one to watch: its fill layout is exactly what was
crushing the enrichment and variant-explorer plots, so a layout verified against 0.11.0
locally is not strictly a verification of what users load. Checking a plot on dev after
deploying stays necessary — it is not belt-and-braces.

**`Matrix` is installed twice in every image**, once as the copy shipped with R and once
upgraded into the site library. `installed.packages()` returns a row per library path,
so a naive capture lists both versions with no way to tell which `library()` gets. The
`version` column here is the one that wins on `.libPaths()` order; the losing copy is
kept in `shadowed_versions` rather than dropped, because a package resolving to an older
copy than you think is a real failure mode and hiding it would defeat the point.
