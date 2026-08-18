# E2F7/8 heart scRNA-seq — interactive cell browser

A server-side [Shiny](https://shiny.posit.co/) app for exploring single-cell
RNA-seq of E2F7/8 knockout (KO) vs wild-type (WT) developing mouse heart at P0
and P7. Colour the UMAP by any gene or metadata, compare KO vs WT, inspect
differential expression by cell type and by cardiomyocyte subcluster, view
subcluster identity / cell-cycle state, and browse pathway/enrichment results.

> **Descriptive pilot — n = 1 per condition, sex-confounded, KO not transcript-confirmed.**
> All KO-vs-WT differences are hypothesis-generating only. See the **About / caveats**
> tab in the app.

The repo also holds a **mechanistic model** of the cardiomyocyte cell-cycle fate
decision, built from two FUCCI papers — see [model/](model/) and its
[RESULTS.md](model/RESULTS.md) / [TODO.md](model/TODO.md). The app is the
descriptive half; `model/` is the quantitative half.

## Layout

```
shiny_app/
  app.R           # the whole app (UI + server)
  app_data.rds    # data bundle (built by build_app_data.R) — NOT in git, see below
  rsconnect/      # shinyapps.io deployment record
model/
  cmcycle/        # cell-cycle fate model + re-analysis of the FUCCI papers (Python, stdlib only)
  figures/        # generated SVGs
  tests/          # 32 tests
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

### Per-subcluster enrichment (`shiny_app/build_subcluster_enrichment.R`)

The cardiomyocyte deep-dive's per-subcluster GO/GSEA tabs and summary sheet read
`app$enrich$sub`, added by `build_subcluster_enrichment.R`. Unlike `build_app_data.R`
this runs **locally against the existing bundle** (no Seurat objects needed): it
reads `app_data.rds`, computes per-res-0.2-subcluster KO-vs-WT GO/GSEA (from the
in-bundle `sub_DE` tables) and cluster-identity GO (from markers detected on the
broad `deg_expr` matrix), and rewrites `app_data.rds` in place (backing up to
`app_data.pre_subenrich.bak.rds` first). Needs `clusterProfiler`, `org.Mm.eg.db`,
`fgsea`, `msigdbr`, `presto`. Re-run it whenever `build_app_data.R` regenerates the
bundle, then redeploy:

```r
# from the repo root, in a real R session (not a sandbox — the 103 MB rds is large)
source("shiny_app/build_subcluster_enrichment.R")   # or: Rscript shiny_app/build_subcluster_enrichment.R
```

### Extra bundle-only analyses (module scores, communication, annotation check)

Three further builders run **locally against the existing bundle** (like
`build_subcluster_enrichment.R` — no Seurat objects needed). Each reads `app_data.rds`,
computes new results, backs up to its own `.bak.rds`, and rewrites `app_data.rds` in
place. All are **descriptive / hypothesis-generating only** (n = 1, sex-confounded).

- **`build_signature_scores.R`** — per-cell module scores (base R + Matrix, no new
  packages): cardiomyocyte proliferation / cytokinesis / cell-cycle-exit / a
  polyploidization proxy, plus CM maturation and a glycolysis→FAO metabolic-switch
  score. Written as new `app$meta` columns (`sig_*`) that surface automatically in the
  UMAP "colour by" dropdown and drive the **Cell-cycle exit & ploidy** and
  **Maturation & metabolism** tabs. Run once with `--probe` first to print gene coverage
  in the curated vs broad matrices and exit without writing.
- **`build_communication.R`** — curated ligand→receptor scoring (VEGF/angiogenesis
  focus) → `app$commun`, drives the **Cell–cell signalling** tab. A descriptive
  NATMI-style score, not a permutation-tested CellChat/LIANA run.
- **`build_refmap.R`** — reference-marker annotation concordance → `app$refmap`, drives
  the **Annotation check** tab (Help menu). A marker-signature check, not Seurat anchor
  transfer.
- **`build_fourgroup.R`** — the four-group (WT-P0 / WT-P7 / KO-P0 / KO-P7) analysis within
  CM subclusters → `app$fourgroup`, driving the **Four-group (WT/KO × P0/P7)** and
  **Maturation ∩ P7 KO** tabs. Per-subcluster cell counts, cell-cycle phase composition,
  maturation-score summaries, descriptive Wilcoxon DE for four contrasts × two phase strata,
  a gene-level maturation axis, and its intersection with the P7 KO response. Also adds a
  `cm_subcluster` column to `app$meta`, which makes CM subcluster a filter in the
  **Subset & DEGs** tab and a colour/split option on the UMAP.
  Must run **after** `build_signature_scores.R` — it consumes the `sig_*` columns.

  Two design points worth knowing:
  - **P0-vs-P7 contrasts are phase-matched by default.** P7 was FACS cycling-enriched
    4.5–5.2× and P0 essentially unenriched, so a raw P0-vs-P7 contrast reads out the sort
    as much as development. Every contrast exists in a `G1` stratum (the app's default) and
    an `all` stratum (raw, labelled sort-confounded in the UI).
  - **The maturation axis is cycle-free.** `sig_maturation`'s immature program contains
    `Mki67` / `Top2a` / `Ccnd1`, so using it to argue "less mature ⇒ more cycling-competent"
    would be partly circular. `build_signature_scores.R` now also emits
    `sig_maturation_nocc` (those three dropped) and the intersection defaults to it.

  Arms too thin to support a contrast are flagged rather than silently reported — note that
  an arm can clear the 10-cell floor overall and still be a handful of cells once restricted
  to G1 (CM2's KO-P0 arm is 31 cells, ~12 of them G1).

Re-run order after `build_app_data.R` regenerates the bundle (each is idempotent and
safe to skip; the app guards absent slots with a "run the builder" message):

```r
# real R session, from the repo root
Rscript shiny_app/build_signature_scores.R --probe   # coverage gate (writes nothing)
Rscript shiny_app/build_signature_scores.R           # then compute + save
Rscript shiny_app/build_communication.R
Rscript shiny_app/build_refmap.R
source("shiny_app/build_subcluster_enrichment.R")    # (existing) per-subcluster enrichment
Rscript shiny_app/build_fourgroup.R --probe          # group sizes + size estimate, writes nothing
Rscript shiny_app/build_fourgroup.R                  # then compute + save
# then: rsconnect::deployApp("shiny_app")
```

`build_fourgroup.R` runs DE on the broad `app$deg_expr` matrix by default; pass
`--matrix=curated` to use the full-cell curated panel instead, or `--seurat=<path>` to
compute over all genes × all CM cells from the upstream Seurat object (the only route to a
genome-wide DEG table — `deg_expr` is downsampled to ~8k cells, which thins the smaller
subcluster arms considerably).

## History

This project previously also shipped as a [shinylive](https://posit-dev.github.io/r-shinylive/)
static (WebAssembly/webR) site served via GitHub Pages. That export is archived
at the git tag **`shinylive-static-archive`**:

```bash
git checkout shinylive-static-archive   # recover app.json + shinylive/ runtime
```
