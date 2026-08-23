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

## Methods notebooks (`docs/`)

A Quarto book documenting every analysis behind the app — one chapter per analysis,
each stating the question, the method with its real parameters, the result, the
reasoning behind the method choice, and what the result cannot support. Written for
reviewing the analysis rather than operating the app.

```bash
docs/render.sh                       # -> docs/_site/index.html  (Quarto in Docker, no R needed)
docs/export.sh                       # regenerate figures + generated tables from the bundle
python3 tools/check_docs_coverage.py  # asserts all 18 app tabs are documented exactly once
```

Rendering executes nothing: every chapter is `engine: markdown`, so the render image
carries the Quarto CLI and no R at all. Figures are exported once, offline, by
`docs/export_assets.R`, which loads the app's **own** plotting functions out of `app.R`
and calls them against the bundle (mounted read-only) — so a figure in the book is the
figure in the app, and every number in a table is read from the data rather than typed.
See [docs/README.md](docs/README.md).

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

### Where the four contrasts show up in the app

They are reachable from **two** places, deliberately.

**Cardiomyocytes → Cardiomyocyte deep-dive → DE (per subcluster)** is where you land when
you are reading one subcluster. Its *Comparison* dropdown carries the pooled KO-vs-WT
tables (`app$tables$sub_DE`, the original behaviour) plus the four timepoint-specific
contrasts from `app$fourgroup`, and the volcano, table, gene card, the top-genes ×
subclusters heatmap and the *Subcluster enrichment* tab all follow it. Note the
*Split map by* control next to it facets the **map only** — the pooled `sub_DE` tables had
their timepoint dimension collapsed at build time, so no display-time split could ever
recover a per-timepoint contrast. That is why the comparison is its own dropdown.

**Cardiomyocytes → Four-group (WT/KO × P0/P7)** remains the place for everything that is
about the four groups rather than about one subcluster: group sizes and the underpowered-arm
audit, the G1/maturation panels, the contrast and enrichment XLSX workbooks, and the
Coverage audit behind an empty GO result.

Both read the same `app$fourgroup$de` / `$de2` grids, so the same cluster × contrast ×
stratum × matrix must give the same table on either tab. `tools/test_downloads.R` asserts
that equality; if it ever fails, one of the two tabs is answering a different question than
its label claims.

Per-subcluster enrichment follows the same split. `app$enrich$sub` (pooled KO-vs-WT) and
`app$enrich$fourgroup` (per contrast) are offered under one *Comparison* dropdown on the
deep-dive's enrichment tab; both enrich the two directions **separately**, so that tab has
a *GO — up* and a *GO — down* panel rather than one showing only the up direction.

### WT programs ∩ KO clusters

A tab for one specific request: cross the WT P0→P7 change in a curated gene category
against the P7 KO-vs-WT response in a named group of CM subclusters, four ways.

- **Step 1** — WT P0→P7 (default all cardiomyocytes), split by curated category
  (`CM maturation`, canonical cell cycle = S ∪ G2/M ∪ E2F targets) and direction.
- **Step 2** — P7 KO-vs-WT in the *maturation* clusters (CM1/2/3/7/8) and the *cycling*
  clusters (CM2/4/5), per cluster and as a gene × cluster pivot.
- **The four comparisons** — each WT category list crossed against the **opposite**
  category's cluster group, drawn with the Gene-set Venn tab's `vn_plot()` and scored with
  its hypergeometric, plus a per-gene table carrying both sides' evidence on one row.

The KO side is defined by **clusters, not by a gene category**. Filtering both sides by
category would make comparisons 1 and 2 intersect the maturation set with the cell-cycle
set, and those two share only `Mki67` and `Top2a` — the Venns would read empty by
construction rather than by biology. `CM2` is in both cluster groups, so the two KO unions
are not independent; the tab says so.

Three things are easy to get wrong here, and the tab reports each rather than leaving it
to be reconstructed:

- **The effect-size measure decides the answer.** `presto`'s log2FC is a difference of mean
  log-normalised expression, so it scales with expression level. From WT P0 to P7 `Mcm3`
  quadruples its detection rate (6.6 % → 27.2 %) for a log2FC of 0.16, while `Myh7` moves
  2.5. A symmetric `|log2FC|` cut therefore classifies maturation genes and can **never**
  classify a cell-cycle one — it reports "no cell-cycle gene changes", which is false. The
  tab defaults to **AUC ≥ 0.60**, rank-based and scale-free, the same resolution
  `build_fourgroup.R` reached for its maturation axis (`MAT_AUC`). Step 1 always prints
  what *both* measures would have given, so an empty list is never mistaken for a result.
- **`CM4` has no G1 stratum** — for any contrast. Choosing the phase-matched stratum drops
  it from the cycling group, so the tab names the dropped clusters instead of quietly
  computing a two-cluster answer.
- **The KO contrast is strongly one-directional.** At AUC ≥ 0.60, CM1 has ~1,700 genes down
  in KO and ~50 up. Seven independently clustered populations do not agree that hard by
  biology; it is the same one-directional signature the 2026-08-21 deliverable traced to a
  library read-fraction difference between the two samples. The tab flags it when the median
  up/down ratio exceeds 5×, and hides mt- genes by default for the same reason.

Set sizes are lopsided by design — one curated category (a few genes) against a whole
cluster group (hundreds) — so the fold enrichment and hypergeometric p on the **Overlap
statistics** tab, not the picture, are what carry the result.

## Docker: the three images, and the lab-server dev deploy

**There is no R on this host.** Everything — running the app, parsing `app.R`, running the
tests — goes through one of three images. Knowing which does what saves a lot of time.

| image | has | use it for |
|---|---|---|
| `lab-server-e2f-heart-scrna-dev` | shiny, bslib, ggplot2, Matrix, plotly, DT, svglite, shinycssloaders, **openxlsx**, presto | the dev deploy, and **the only image that can run `tools/test_downloads.R`** |
| `lab-server-e2f-heart-scrna` | same **minus openxlsx** | production. `app.R` calls `library(openxlsx)` at line 20, so the test suite cannot run here |
| `e2f-enrich` | clusterProfiler, org.Mm.eg.db, fgsea, msigdbr, presto | the offline `analysis/` pipeline, plus quick parse and coverage checks (no shiny needed) |

Mount the repo **read-only** and use `--rm` for checks, so a throwaway container can never
touch the working tree:

```bash
# syntax check — the fastest way to catch a broken edit, ~5 s
docker run --rm -v "$PWD:/repo:ro" -w /repo e2f-enrich:latest \
  Rscript -e 'parse("shiny_app/app.R")'

# static download coverage: every rendered table and plot must have a download
docker run --rm -v "$PWD:/repo:ro" -w /repo e2f-enrich:latest \
  Rscript tools/check_download_coverage.R

# runtime suite via shiny::testServer — needs openxlsx, so the dev image. ~3 min
# (it readRDS's the 119 MB bundle once)
docker run --rm -v "$PWD:/repo:ro" -w /repo lab-server-e2f-heart-scrna-dev:latest \
  Rscript tools/test_downloads.R

# Figure Studio handoff: environment stripping, one figspec per figure, TTL. ~3 min
docker run --rm -v "$PWD:/repo:ro" -w /repo lab-server-e2f-heart-scrna-dev:latest \
  Rscript tools/test_studio_handoff.R
```

Ad-hoc poking at the bundle follows the same shape — write the script to a scratch dir,
mount it, do **one** `readRDS` per script:

```bash
docker run --rm -v "$PWD:/repo:ro" -v /tmp/scratch:/sp:ro \
  lab-server-e2f-heart-scrna-dev:latest Rscript /sp/probe.R
```

### The lab-server stack

The app is served by a separate compose project at **`/home/justin/Projects/lab-server`**,
which is not part of this repo. Two services run the same app from two payload directories:

| service | payload | route | restart |
|---|---|---|---|
| `e2f-heart-scrna` | `apps/e2f-heart-scrna/` | `/e2f-heart-scrna/` | `unless-stopped` |
| `e2f-heart-scrna-dev` | `apps/e2f-heart-scrna-dev/` | `/e2f-heart-scrna-dev/` | **`"no"`** — up only while someone is testing |

Both are `FROM rocker/shiny:4.5.1` with `COPY . /srv/shiny-server/`, so **the app is baked
into the image — there are no bind mounts and nothing updates without a rebuild.** nginx
strips the path prefix before proxying to `3838`, which is why shiny-server serves at `/`.
Separate payload dirs are the point: a dev experiment cannot disturb production or its
bundle.

Everything behind the front door sits behind the Han Lab sign-in (the `labgate` cookie), so
`curl http://localhost/e2f-heart-scrna-dev/` returns the login page, not the app. To check
the app itself, go **inside** the container (below).

### Deploying to dev

```bash
# 1. copy the payload across. Copy ALL of shiny_app/*.R, not a hand-listed subset:
#    app.R source()s its helpers at startup and reads build_signature_scores.R for the
#    curated gene lists, so a missing .R file is a container that will not boot. This has
#    already happened once. app_data.rds is left alone — the payload has its own copy.
cp shiny_app/*.R /home/justin/Projects/lab-server/apps/e2f-heart-scrna-dev/

# 2. rebuild and restart (recreates the container; ~30 s, R packages are cached)
cd /home/justin/Projects/lab-server
docker compose up -d --build e2f-heart-scrna-dev

# 3. stop it when you are done testing — restart:"no" means it will not come back by itself
docker compose stop e2f-heart-scrna-dev
```

**Verify the deploy took, rather than trusting the build log.** A successful build says
nothing about whether the file you meant to ship is the one inside:

```bash
# the payload actually in the running container must match the working tree
sha256sum shiny_app/app.R
docker exec lab-server-e2f-heart-scrna-dev-1 sha256sum /srv/shiny-server/app.R

# render the real UI from inside the container, past the sign-in, and grep it
docker exec lab-server-e2f-heart-scrna-dev-1 Rscript -e '
  h <- paste(readLines(url("http://localhost:3838/", open="rb"), warn=FALSE), collapse="\n")
  cat(nchar(h), "bytes\n"); cat(grepl("WT programs", h), "\n")'

# and confirm the app started clean
docker logs lab-server-e2f-heart-scrna-dev-1 2>&1 | grep -iE 'error|exception'
```

A healthy load is ~450 KB of HTML and takes ~11 s (the 119 MB `app_data.rds` is read at
startup). The container name is `lab-server-e2f-heart-scrna-dev-1` — compose appends `-1`.

**The dev bundle is its own copy.** `apps/e2f-heart-scrna-dev/app_data.rds` is a separate
119 MB file from `shiny_app/app_data.rds`; copying `app.R` does not update it. Re-copy the
bundle too after any builder run, or the app will show its "run the builder" guards.
Production is further behind still — its payload is a 23 MB bundle and a 105 KB `app.R`
from 2026-07-21.

### Running the app locally instead

Unrelated to the lab-server stack, the repo ships its own single-service compose file that
mounts the bundle at runtime rather than baking it in. See [DOCKER.md](DOCKER.md):

```bash
docker compose up --build      # http://localhost:3838
```

It needs `shiny_app/app_data.rds` to exist **before** `up` — Docker will otherwise create an
empty directory at the mount point and the app will fail to read it.

## Module scores: how maturation and metabolism are defined

Every `sig_*` score is built from a **hand-picked list of canonical marker genes, hardcoded
in [`build_signature_scores.R`](shiny_app/build_signature_scores.R)**. Nothing here is
data-driven and nothing comes from a database — MSigDB and GO appear in this repo only in
the *enrichment* builders, and never touch a module score.

**On provenance, so nobody over-claims it downstream:** the specific gene choices have no
recorded source. They arrived in one commit (`be1f1bd`, 2026-07-21, "more data anaylsis")
with no rationale and no citation, and a comment on the app's own copy
([app.R:55](shiny_app/app.R#L55)) attributes it to "the pipeline's `_common.R`", a file not
in this repo. The defensible description is *canonical markers curated by the analyst* —
not "derived from *(reference)*". If these lists are going into a figure legend, that
sentence is the one to write.

### Reading it in the app

**Every `sig_*` score documents itself.** Under the violin on **Maturation & metabolism**,
and under the score distributions on **Cell-cycle exit & ploidy**, a collapsed block gives
*what the score is for* (the question it exists to answer), *how it is built* (which two
sets, which matrix, how many genes were actually found), and *the genes themselves* —
naming any set gene absent from the scoring matrix, since those contributed nothing.
**Help → QC & normalization** carries the same lists as one downloadable reference table.

The app does not keep its own copy of those lists. `app_data.rds` never stored them, so
`SCORE_SETS` parses the `SETS <- list(...)` literal straight out of
`build_signature_scores.R` at startup — the same trick
`analysis/2026-08-21_email/01_de.R` uses — and a test asserts the parsed lists equal that
literal. A second hardcoded copy in `app.R` would drift the first time either was edited,
and the app would then confidently document genes it had not scored.



Every figure on **Maturation & metabolism**, and the four Venns on **WT programs ∩ KO
clusters**, carries a collapsed **"How this plot was made"** block underneath it —
what a point is, how the coordinate was computed, what the annotations mean, and the
functions that build it. The blocks are live: the gene map's quotes the selected panel's
actual axis centre, and the violin's names the gene sets behind whichever score is chosen.

They are written as free functions (`gm_method_note()`, `mat_violin_method_note()`,
`mat_scatter_method_note()`, `xc_venn_method_note()`) rather than inline `renderUI` bodies,
so the tests can call them with a panel or a score and check they still follow it —
`testServer` snapshots `output$` values, so a note that quietly stopped tracking its
dropdown would look correct forever.

### The lists

| set | n | genes |
|---|---|---|
| `mat_mature` | 10 | Myh6, Tnni3, Pln, Atp2a2, Ckm, Myl2, Cox6a2, Ckmt2, Actn2, Csrp3 |
| `mat_immature` | 9 | Myh7, Tnni1, Nppa, Nppb, **Ccnd1, Mki67, Top2a**, Myl7, Actc1 |
| `mat_immature_nocc` | 6 | as above minus the three cell-cycle genes |
| `glycolysis` | 11 | Slc2a1, Hk1, Hk2, Pfkm, Pfkl, Pkm, Ldha, Gapdh, Eno1, Aldoa, Pgk1 |
| `faox` | 14 | Cpt1a, Cpt1b, Cpt2, Acadm, Acadvl, Acadl, Hadha, Hadhb, Ppargc1a, Cox6a2, Ndufa4, Sdha, Acaa2, Etfa |
| `prolif` | 17 | Mki67, Top2a, Ccnb1, Ccnb2, Cdk1, Cdc20, Aurka, Aurkb, Bub1, Birc5, Cenpa, Cenpe, Cenpf, Ube2c, Cks2, Nusap1, Tpx2 |
| `cytokinesis` | 14 | Anln, Ect2, Racgap1, Kif23, Cit, Aurkb, Kif20b, Prc1, Cep55, Mklp1, Sept7, Sept9, Cdca8, Incenp |
| `ccexit` | 8 | Cdkn1a, Cdkn1c, Cdkn2a, Cdkn1b, Meis1, Rb1, Btg2, Gadd45a |

The maturation logic is the fetal→adult isoform switches (Myh7→Myh6, Tnni1→Tnni3,
Myl7→Myl2) plus calcium handling (Pln, Atp2a2), energetics (Ckm, Ckmt2, Cox6a2), sarcomere
(Actn2, Csrp3) and the fetal stress markers (Nppa, Nppb). Glycolysis is the pathway walked
end to end; FAO is carnitine shuttle → acyl-CoA dehydrogenases → trifunctional protein →
thiolase.

**There is no separate OXPHOS set.** Ppargc1a, Cox6a2, Ndufa4, Sdha and Etfa sit inside
`faox`, so "metabolic maturation" and "oxidative phosphorylation" are not separable on
these scores.

### How a score is computed

An `AddModuleScore` equivalent in base R + Matrix, no Seurat
([build_signature_scores.R:102](shiny_app/build_signature_scores.R#L102)):

1. Drop `app$confound` (Xist, Tsix, Y-genes, ROSA26) from every set. None of these sets
   contain any, so the scores are confounder-free by construction.
2. Bin all genes into **24** quantile bins by mean expression.
3. For each set gene, sample **100** control genes from its own bin (`seed = 1`).
4. **score = mean(set genes) − mean(pooled controls)**, per cell.

Composites are differences, never ratios: `sig_maturation` = mature − immature,
`sig_metabolic` = faox − glycolysis (positive = more oxidative). No z-scoring, no min-max,
no rescaling — the values in `app$meta` are raw. `pick_mat()` chooses per set whichever
matrix covers more of its genes, forcing both halves of a difference onto the same one, so
scores computed on the broad matrix are `NA` for cells outside its 8,026-cell downsample.
`app$score_meta` records genes-used / genes-in-set / cells-scored per score and is rendered
under **Help → QC & normalization**.

One deviation from Seurat worth knowing: the control pool is **not de-duplicated**, so genes
in dense expression bins carry extra weight in the background. It shifts the baseline, not
the ranking of cells.

### Three traps

**1. There are three different maturation lists in this codebase and they disagree.**

| where | n | difference |
|---|---|---|
| `SETS$mat_mature ∪ mat_immature` — drives the `sig_*` scores | 19 | the full list |
| `FG$built$score_set_genes$maturation` — the gene map | 16 | drops Ccnd1, Mki67, Top2a |
| `GENE_SETS[["CM maturation"]]` ([app.R:76](shiny_app/app.R#L76)) — the app dropdowns, and the **WT programs ∩ KO clusters** tab | 14 | drops Ckmt2, Actn2, Csrp3, Myl7, Actc1 |

This is an inconsistency, not a design choice. Reconcile it before these counts go into a
figure.

**2. `Cox6a2` is in both `mat_mature` and `faox`**, so the maturation and metabolic axes
share an input gene — `build_fourgroup.R` flags it as doubly circular.

**3. `mat_immature` contains Mki67, Top2a and Ccnd1**, so `sig_maturation` is partly a
cell-cycle score. Using it to argue "less mature ⇒ more cycling-competent" is circular;
that is why `sig_maturation_nocc` exists and why everything downstream defaults to it.

### The one place gene selection *is* data-driven

The **Gene map** axes ([build_fourgroup.R:353](shiny_app/build_fourgroup.R#L353)) invert the
logic: rather than scoring cells from genes, they rank *genes* against the per-cell score.
CM cells are split on tertiles of `sig_maturation_nocc` **within each timepoint**,
`presto::wilcoxauc` runs top-third vs bottom-third, and the AUCs are averaged across P0 and
P7 — within-then-average is what stops the axis becoming a P0-vs-P7 axis, which the FACS
sort confounds. Classification is AUC ≥ 0.60 / ≤ 0.40 at padj < 0.05, and yields only **48**
maturation-classified genes out of 11,047. That thin yield is why the curated lists are
still what the tabs use.

### Effect-size measure: why the app offers AUC

`presto`'s `logFC` is a difference of **mean log-normalised expression**, so it scales with
how highly expressed a gene is. Between WT P0 and P7, `Mcm3` quadruples its detection rate
(6.6 % → 27.2 %) for a log2FC of 0.16, while `Myh7` moves 2.5. A symmetric `|log2FC|` cut
therefore classifies maturation genes and can **never** classify a sparse cell-cycle one —
it returns "no cell-cycle gene changes", which is false. Measured on this bundle
(AllCM, all cells, padj < 0.05):

| cut | maturation up / down | cell cycle up / down |
|---|---|---|
| AUC ≥ 0.50 | 4 / 7 | **32** / 1 |
| AUC ≥ 0.55 | 4 / 5 | **14** / 0 |
| AUC ≥ 0.60 *(default)* | 2 / 5 | **3** / 0 |
| AUC ≥ 0.65 | 1 / 4 | 0 / 0 |
| `|log2FC|` ≥ 0.25 | 1 / 4 | **0** / 0 |
| `|log2FC|` ≥ 1.0 | 0 / 2 | 0 / 0 |

Both overlap tabs — **Gene-set Venn** and **WT programs ∩ KO clusters** — expose the choice
as a radio plus a slider, because the answer moves across the plausible range and the cut is
something to sweep rather than set once. `build_fourgroup.R` reached the same resolution
independently for its maturation axis (`MAT_AUC = 0.60`).

Separately, every volcano has a **"Colour genes at |log2FC| ≥"** slider (default 1). That is
a *display* threshold — it decides which points are coloured up/down versus grey, and never
filters the plot or the table. It defaulted to a hardcoded 1, which on the P7 KO-vs-WT
contrasts colours almost nothing, because at `|log2FC| ≥ 1` those contrasts have **zero**
genes down in KO.

An unrelated set that is easy to confuse: `model/cmcycle/baniol.py` carries its own
FAO/glycolysis pair for the Baniol FUCCI re-analysis, and the genes differ (Hmgcs2, Fabp3,
Pdk4, Ucp2, Cd36, Ech1, Decr1 …). Different dataset, different index — do not cross-cite.

## History

This project previously also shipped as a [shinylive](https://posit-dev.github.io/r-shinylive/)
static (WebAssembly/webR) site served via GitHub Pages. That export is archived
at the git tag **`shinylive-static-archive`**:

```bash
git checkout shinylive-static-archive   # recover app.json + shinylive/ runtime
```
