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

  It ships **two DE grids** over the same contrasts, because the two available matrices
  trade against each other and neither wins outright. `app$fourgroup$de` uses the broad
  matrix (~24k genes, but a downsampled ~5.8k CM cells, so CM2's KO-P0 arm falls to 9
  cells and drops out and CM4/CM9 lose their G1 strata); `app$fourgroup$de2` uses the
  curated panel (2,181 genes but all 21,598 CM cells, so every contrast runs). Build the
  second with `Rscript shiny_app/build_fourgroup.R --de2` **after** the main run — it
  recomputes only the DE grid and leaves the rest of the slot alone. The app offers both
  as a "DE matrix" choice and, when a contrast is missing from one, says whether the
  other has it. Across the four priority subclusters this takes coverage of the twelve
  requested contrasts from 7/12 to 12/12.

  The app's **Gene-set Venn** tab crosses any two or three of these sets — DE contrasts,
  the maturation/metabolic/cycling axes, the intersection quadrants, or the curated panels —
  and reports every pairwise overlap against its chance expectation, because a Venn on its
  own hides the threshold that built each set, the direction of change, and the null. Its
  cycling circle defaults to the **curated canonical** genes rather than the data-driven
  axis: the axis calls 531 genes cycling-associated but only 46 are canonical, the rest
  being largely housekeeping (`Ran`, `Nap1l1`, `Calm1`, `Ppia`) because cycling cells are
  globally more transcriptionally active. The universe used for the hypergeometric is the
  tested gene space (`app$deg_genes` or the curated panel), **not** the gated DE tables —
  those exclude expressed-but-unchanging genes, which is exactly what a universe needs.

  Arms too thin to support a contrast are flagged rather than silently reported — note that
  an arm can clear the 10-cell floor overall and still be a handful of cells once restricted
  to G1 (CM2's KO-P0 arm is 31 cells, ~12 of them G1).

  It also emits **`app$fourgroup$geneaxes`**, the gene map behind the *Gene map* sub-tab of
  **Maturation & metabolism**. Every gene gets two coordinates, each an AUC from a tertile
  split of the cells computed **within each timepoint and averaged** (so neither axis becomes
  a restatement of P0-vs-P7, which the cycling sort confounds):

  - `mat_auc` — how strongly the gene marks mature vs immature cardiomyocytes
  - `met_auc` — how strongly it marks oxidative (FAO) vs glycolytic metabolism

  plus `quadrant`, `distance` from the centre (the ranking of how strongly a gene defines
  the joint program), and `in_score_set`. Three things are easy to get wrong here:

  - **The axes split at each one's median, not at 0.5.** `wilcoxauc`'s AUC carries a small
    global offset because the two tertile groups differ in detection rate (here mat +0.009,
    met −0.014). Most genes sit within ~0.02 of the median, so splitting at a hard 0.5 put
    65% of them in a single corner — an artifact, not biology. The centre is stored on the
    table as `attr(geneaxes, "centre")`.
  - **The axes use a looser row gate than the DE tables** (expression level only, no
    significance requirement). The DE gate exists to keep 77 tables small; applying it to the
    axes dropped Gapdh, Aldoa, Pgk1, Eno1, Hk1 and Cpt1a off the map entirely. A gene with no
    association belongs at the origin, not missing.
  - **`in_score_set` flags circularity.** Genes inside the sets that define a score sit at the
    extremes of their own axis by construction; the app hides them by default. `Cox6a2` is in
    both `mat_mature` and `faox`, so it is doubly circular. With set genes hidden the
    mature↔oxidative / immature↔glycolytic diagonal still holds 57% of genes (50% would mean
    the two axes are independent), so the coupling is not an artifact of the inputs.

- **`build_fourgroup_enrichment.R`** — GO (BP/MF/CC) and Hallmark+KEGG GSEA for the
  **four-group contrasts** → `app$enrich$fourgroup`, driving the **Enrichment** panel of the
  Four-group tab. Covers all four contrasts × both strata × every subcluster (77 tables,
  154 direction-tests), enriching each direction **separately**.

  This is not the same question as `build_subcluster_enrichment.R`, and that is the point:
  `app$enrich$sub` is KO-vs-WT **pooled across P0 and P7**, which cannot answer "what changes
  in the KO *at P7*". The DE for the timepoint-specific contrasts has existed in
  `app$fourgroup$de` since `build_fourgroup.R`; GO had simply never been run on it.

  One thing to know before reading its output: **the universe is recomputed from the
  expression matrix, not taken from the DE table.** `app$fourgroup$de` is row-gated (expressed
  in ≥ 5 % of one arm AND (padj < 0.05 OR |log2FC| ≥ 0.5)), so most expressed-but-unchanging
  genes — exactly the background a hypergeometric test needs — are already missing from it.
  Using it as the universe would inflate every fold enrichment. The builder instead takes
  genes detected in ≥ 5 % of the cells of one arm of that specific contrast, cluster and
  stratum: ~11,000 genes rather than the ~3,000 the gated table would have supplied.

  Direction labels come from `FG$built$contrasts$pos`/`$neg`, so a `WT: P0 vs P7` panel reads
  "up at P7" / "up at P0" rather than the KO-vs-WT wording. `gsea_barplot_gg()` takes
  `up_lab`/`down_lab` for the same reason.

  `--probe` reports input and universe sizes and exits; `--ont=BP` restricts the ontologies;
  `--contrast=<key>` restricts to one contrast; `--grid=de2` uses the curated-panel DE grid.

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
Rscript shiny_app/build_fourgroup.R --de2            # second DE grid on the curated panel
Rscript shiny_app/build_fourgroup_enrichment.R --probe   # input/universe sizes, writes nothing
Rscript shiny_app/build_fourgroup_enrichment.R           # GO + GSEA for the four-group contrasts
# then: rsconnect::deployApp("shiny_app")
```

`build_fourgroup_enrichment.R` is the slow one (~460 `enrichGO` calls over 77 DE tables,
a couple of hours). Everything else finishes in minutes.

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
