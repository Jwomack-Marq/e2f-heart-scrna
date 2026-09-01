# Normalization & preprocessing — procedure and current-vs-previous comparison

Reference summary of how the **current** pipeline (`rerun/`) normalizes and preprocesses
the E2F7/8 KO mouse-heart scRNA-seq data, with a side-by-side comparison to the
**previous** versions. Focuses on the pieces that cause the most confusion: **lanes,
tags, scrublet, and the SCTransform choices.** Companion docs: [README.md](README.md)
(pipeline + run order), [METHODS_COMPARISON.md](METHODS_COMPARISON.md) (full method diff),
[REPLICATES.md](REPLICATES.md) (n=1 / sex / ROSA26 caveats).

---

## 0. The data, in one line

**8 libraries = 2 timepoints (P0, P7) × 2 genotypes (WT, KO) × 2 sequencing lanes (lane1, lane6).**
Understanding "lanes" is the key to everything below.

## 1. Pipeline order (current)

```
PIPseeker v3.3.0  →  Scrublet (per lane)  →  QC filter + merge lanes  →  SCTransform
(cell calling)       (doublet detection)     (per group)                (glmGamPoi)
   →  per-condition merge  →  4-group Harmony integration  →  clustering + UMAP
   →  marker-based cell-type annotation  →  CellCycleScoring
```
Scripts: `scripts/pipseeker_scripts.txt`, `scripts/scrublets.py`,
`scripts/scRNA_lane_merge.Rmd`, `scripts/scRNA_mergeP0.Rmd` / `scRNA_mergeP7.Rmd`,
`rerun/analysis/combined.R`, `rerun/analysis/annotate.R`, `rerun/analysis/cell_cycle.R`.

## 2. Lanes  — the most important concept

- **What they are here:** `lane1` and `lane6` are the **same library (same animal, same
  PIPseq reaction) sequenced on two flow-cell lanes** to add depth — *not* two animals.
  Proof: both lanes of a condition point to the **identical Novogene library ID** (e.g.
  `P0WT_CKDL250003755`) and barcode overlap between lanes is **97–100%** (see
  `rerun/analysis/tables/lane_barcode_overlap.csv`). Lab-confirmed (June 2026).
- **Consequence: n = 1 biological replicate per condition.** All KO-vs-WT differential
  expression is **descriptive / hypothesis-generating only** — no method produces valid
  p-values at n=1.
- **How current handles them:** each lane is loaded separately (`Read10X`), QC-filtered,
  tagged with a `$lane` metadata field, then `merge()`d within a condition (v5: `JoinLayers`).
  The `$lane` field lets pseudobulk treat the two lanes as two (technical) pseudo-replicates.
- **Known defect (documented, not yet fixed):** because each lane was counted by a *separate*
  PIPseeker run and then merged, ~every cell appears **twice at half its true depth** → cell
  counts ~2× inflated; the `nFeature ≥ 1500` cutoff was applied to half-depth cells. Standard
  fix: pool both lanes' FASTQs and run PIPseeker **once per library**
  (`scripts/recount_combine_lanes.sh`, runs on the cluster). See REPLICATES.md.
- **Previous:** the original run **did not track lane of origin** (no `$lane`), so it couldn't
  do per-lane pseudobulk and the double-counting was invisible.

## 3. Tags — there are none (no cell hashing)

This is a droplet **PIP-seq** protocol (Fluent PIPseeker), **not** a multiplexed cell-hashing
experiment. **There is no HTO / antibody-tag demultiplexing step** in either pipeline. The only
"tags" present are:
- **cell barcodes** (PIPseeker cell calling, from `filtered_matrix/sensitivity_5/`), and
- **sample labels** attached as metadata: `orig.ident` (P0WT/P0KO/P7WT/P7KO), `timepoint`,
  `genotype`, `lane`.

If someone says "look at the tags," they mean these sample/lane labels, not antibody hashtags.

## 4. Scrublet — doublet detection

- **What it does:** simulates artificial doublets by combining random pairs of real
  transcriptomes, then scores each real cell by similarity to those simulated doublets.
  High score → a droplet that likely captured two cells.
- **Current** (`scripts/scrublets.py`): run **once per lane** (8 runs) on the raw matrix.
  Parameters: `expected_doublet_rate=0.08`, `min_counts=2`, `min_cells=3`,
  `min_gene_variability_pctl=85`, `n_prin_comps=30`. Writes `predicted_doublets.csv`
  (`doublet_scores`, `predicted_doublets`) + a score histogram. Flags are attached to each
  Seurat object as `is_doublet` and **doublets are actively removed** at the per-lane QC step
  (`subset(... is_doublet == FALSE)`). The rerun also cross-checks recall with **scDblFinder**
  (`rerun/analysis/qc_doublets.R`) — Scrublet's ~1–2% called rate may under-call vs the ~8% expected.
- **Previous:** Scrublet was run, but **doublets were flagged and never removed** — one of the
  headline bugs the rerun fixed.

## 5. Normalization & QC

| Step | **Current (`rerun/`)** | **Previous — original run** | **Previous — MCW Quarto "upgrade" (intended)** |
|---|---|---|---|
| Counts | PIPseeker `filtered_matrix/sensitivity_5` | same | `filtered_feature_bc_matrix` |
| Mito QC | **`^mt-` (mouse), `percent.mt ≤ 20`** | **none** (disabled line used human `^MT-` → matched nothing) | `percent.mt ≤ 10` (placeholder) |
| nFeature | ≥1500 **and** ≤99.5th pctile/lane | ≥1500 only | 500–6000 (placeholder) |
| Doublets | Scrublet, **removed** + scDblFinder check | Scrublet, **flagged only** | scDblFinder, removed |
| Normalization | `SCTransform(method="glmGamPoi")`, **no `vars.to.regress`** | `SCTransform(glmGamPoi)` (similar) | `SCTransform(vars.to.regress="percent.mt")`, 3000 integration features |
| Integration | 4-group **Harmony** on `orig.ident`, dims 1:30, **res 0.8** | per-condition `merge()`, **no batch correction** | Harmony on `sample`, **res 0.6** |
| Annotation | marker-module argmax (SingleR attempted; celldex blocked) | MAESTRO TabulaMuris | SingleR `MouseRNAseqData` |
| Cell cycle | `CellCycleScoring`, human `cc.genes.updated.2019` → mouse via `babelgene` | original `Mus_musculus.csv` list | `CellCycleScoring` |
| KO-vs-WT DE | pseudobulk DESeq2 + apeglm `lfcShrink`, **labeled descriptive(n=1)** | single-cell Wilcoxon (pseudoreplication) | DESeq2 pseudobulk, but skips at n<2 → no output |

No ribosomal-content filter and no `nCount_RNA` minimum are applied (complexity is captured by `nFeature_RNA`).

### The one normalization choice worth flagging: `percent.mt` regression
The MCW pipeline **regresses out `percent.mt`**; the current pipeline deliberately **does not**,
and instead caps it generously at 20%. In heart tissue this is intentional — cardiomyocytes are
*genuinely* mitochondria-rich, so much of `percent.mt` variance is real CM biology, and regressing
it can blunt true CM signal. Both choices are defensible; at n=1 it changes no biological
conclusion. (Full discussion in METHODS_COMPARISON.md.)

## 6. Downstream of normalization

- **Integration:** SCTransform → PCA → `RunHarmony(group.by.vars = "orig.ident")` →
  `FindNeighbors`/`FindClusters(resolution = 0.8)` (dims 1:30) → UMAP. Harmony is for the
  **embedding/annotation only**; all DE stays on the RNA assay / pseudobulk.
- **Annotation:** canonical-marker module-score argmax per cluster → `$celltype`
  (`rerun/analysis/_common.R::CELLTYPE_MARKERS`). SingleR is kept as an optional cross-check
  column but celldex's dependency stack would not install.
- **Cell cycle:** `CellCycleScoring` using Seurat's `cc.genes.updated.2019` mapped human→mouse
  via `babelgene` (title-case fallback). Adds `S.Score`, `G2M.Score`, `Phase`.

## 7. Caveats that gate all KO-vs-WT interpretation

1. **n = 1 per condition** (lanes are the same library) → DE is descriptive, no valid p-values.
2. **Sex confound:** KO = male, WT = female (Y-linked `Eif2s3y/Kdm5d/Uty/Ddx3y` top the KO-up list).
3. **KO allele:** `E2f7`/`E2f8` are **not** transcript-reduced in KO and `Gt(ROSA)26Sor` is KO-up
   → likely a ROSA26 conditional allele that 3′ scRNA can't see; confirm the construct with the lab.

See REPLICATES.md for the replicated-cohort spec (≥3 sex-matched animals/condition).

---

## 8. The PREVIOUS procedure (original run), described on its own

> "Previous" has two meanings here. This section describes **the project's own original
> analysis** — the first Seurat-Rmd run that the rerun replaced (the gaps it had are recorded
> in README.md "Changes made"). A *separate* MCW Quarto "upgrade" attempt (mostly a
> non-executed scaffold) is the other reference point; its distinct choices
> (`vars.to.regress="percent.mt"`, 3000 integration features, scDblFinder, SingleR,
> Harmony at res 0.6) are detailed in METHODS_COMPARISON.md. Everything below is the
> **original run**, so you can read it against §2–§6 above.

**One line:** a straightforward per-lane → `merge()` Seurat-Rmd pipeline with **no batch
correction** and several QC steps that were either missing or broken — later fixed in the rerun.

- **Upstream (unchanged):** same PIPseeker `filtered_matrix/sensitivity_5` matrices. The
  difference is entirely in the R/Seurat steps, not in cell calling.
- **Lanes:** lanes were `merge()`d within a condition but **lane of origin was NOT tracked**
  (no `$lane` field). No batch correction across lanes/conditions, and no awareness of the
  97–100% barcode double-counting.
- **Tags:** none (same — PIP-seq, no HTO).
- **Doublets:** Scrublet **was run**, but the calls were **flagged only and never removed**
  (`is_doublet` computed, but no `subset()` to drop them). No scDblFinder cross-check.
- **Normalization (core method):** **`SCTransform` — essentially the same call as the rerun.**
  The normalization *method* is not where the versions differ; the differences are in QC,
  doublet handling, lane tracking, integration, cardiac selection, and DE (below).
- **QC filtering:** **`nFeature_RNA ≥ 1500` only.** **No mitochondrial filter** — the mito line
  was disabled *and* used the human `^MT-` regex (which matches nothing in mouse), so
  `percent.mt` was effectively unused. **No upper `nFeature` cap.**
- **Clustering:** per-condition `merge()` (KO+WT) → `SCTransform` → PCA → `FindNeighbors` /
  `FindClusters` (dims 1:30) → UMAP. **No Harmony / no integration.**
- **Cardiomyocyte selection:** a **hardcoded UMAP region** (cells within radius 9 of point
  (4, 0)) — not reproducible and lane/embedding-dependent.
- **Cell-cycle genes:** the original `Mus_musculus.csv` ENSEMBL list (via `bitr`/`org.Mm.eg.db`),
  rather than the rerun's `babelgene` human→mouse mapping of `cc.genes.updated.2019`.
- **Annotation:** **MAESTRO** (`mouse.all.droplet.TabulaMuris`, Heart_and_Aorta signatures).
- **KO-vs-WT DE:** **single-cell Wilcoxon `FindMarkers`** (KO vs WT) + `EnhancedVolcano`
  (pCutoff 0.01, FCcutoff 2); GO via `enrichGO` (BP), P0 only. This treats individual cells as
  replicates (**pseudoreplication**) — not valid for condition-level inference.

### What the rerun changed (the deltas to compare)
| Aspect | Original run (previous) | Rerun (current) |
|---|---|---|
| Doublets | flagged, **not removed** | **removed** + scDblFinder check |
| Mito QC | **none** (broken `^MT-`) | `^mt-`, `percent.mt ≤ 20` |
| Upper nFeature | none | ≤ 99.5th pctile/lane |
| Lane tracking | none | `$lane` added (enables pseudobulk) |
| Batch correction | none | 4-group **Harmony** |
| Cardiac selection | hardcoded UMAP circle | marker-module score |
| Cell-cycle genes | `Mus_musculus.csv` | `babelgene` ortholog mapping |
| Annotation | MAESTRO | marker-module argmax |
| KO-vs-WT DE | single-cell Wilcoxon | pseudobulk DESeq2 + apeglm, **descriptive(n=1)** |
| Normalization call | `SCTransform` | `SCTransform(glmGamPoi)` — **essentially unchanged** |

**Bottom line:** the previous run's *normalization* (SCTransform) is basically the same as the
rerun's; the meaningful differences are the **missing mito QC, un-removed doublets, untracked
lanes, no batch correction, a hardcoded cardiac gate, and a pseudoreplicated DE test** — all
addressed in the rerun.

---
*Generated as a reference summary; sources are the scripts and docs cited above.*
