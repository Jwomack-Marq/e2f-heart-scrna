# E2F7/8 heart scRNA-seq — interactive cell browser

A server-side [Shiny](https://shiny.posit.co/) app for exploring single-cell
RNA-seq of E2F7/8 knockout (KO) vs wild-type (WT) developing mouse heart at P0
and P7. Colour the UMAP by any gene or metadata, compare KO vs WT, inspect
differential expression by cell type and by cardiomyocyte subcluster, view
subcluster identity / cell-cycle state, and browse pathway/enrichment results.

> **Descriptive pilot — n = 1 per condition, sex-confounded, KO not transcript-confirmed.**
> All KO-vs-WT differences are hypothesis-generating only. See the **About / caveats**
> tab in the app.

## Layout

```
shiny_app/
  app.R           # the whole app (UI + server)
  app_data.rds    # data bundle (built by build_app_data.R) — NOT in git, see below
  rsconnect/      # shinyapps.io deployment record
```

At runtime the app needs only `app_data.rds` (no Seurat / source objects). The bundle
is a **build artifact**: gzip-compressed it is ~108 MB, over GitHub's 100 MB file
limit, so it is **git-ignored** and shipped to shinyapps.io via `rsconnect` (working
tree), not committed. Regenerate it with `build_app_data.R` from the analysis pipeline
(and `build_expr_full.R` / the enrichment scripts) — gzip, not xz, for faster
cold-start deserialization.

## Run locally

```r
install.packages(c("shiny", "bslib", "ggplot2", "Matrix", "plotly", "DT",
                   "svglite", "shinycssloaders"))
# presto (descriptive Wilcoxon for the interactive "Subset & DEGs" tab):
remotes::install_github("immunogenomics/presto")

shiny::runApp("shiny_app")   # needs shiny_app/app_data.rds present
```

## Deploy (shinyapps.io)

```r
rsconnect::deployApp("shiny_app")   # account: jwomackmu, app: e2f-heart-scrna
```

## Data pipeline

`app_data.rds` is produced by the upstream analysis pipeline (not in this repo;
see `our_analysis/06_outputs/app/build_app_data.R` in the project workspace).
It bundles a downsampled cell × metadata table, a curated expression panel plus a
broader matrix for on-the-fly DE, precomputed cell-type / subcluster DE tables,
marker heatmaps, enrichment results, and per-gene info.

## History

This project previously also shipped as a [shinylive](https://posit-dev.github.io/r-shinylive/)
static (WebAssembly/webR) site served via GitHub Pages. That export is archived
at the git tag **`shinylive-static-archive`**:

```bash
git checkout shinylive-static-archive   # recover app.json + shinylive/ runtime
```
