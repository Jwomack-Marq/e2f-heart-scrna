# Project overview — scRNA-seq of E2F7/8 Knockout Mouse Heart (Han, 2025)

> **Navigation:** the top-level [`README.md`](../../README.md) explains the
> `original_Han_analysis/` vs `our_analysis/` split; [`RUN_ORDER.md`](RUN_ORDER.md)
> is the step-by-step reproduction guide. This file is the biological + methods
> background.

Single-cell RNA-seq analysis of **E2F7/8 knockout in the developing mouse heart**.
Eight libraries: **2 timepoints (P0, P7) × 2 conditions (WT, KO) × 2 sequencing lanes (lane1, lane6)**.

## Biological purpose

E2F7 and E2F8 are atypical *repressor* E2F transcription factors implicated in
cardiomyocyte cell-cycle exit and polyploidization. This study asks how their
loss (KO vs WT) changes cardiac cell populations and gene expression across the
perinatal window — **P0 (birth)** and **P7 (postnatal day 7)**. The central
readout is **KO-vs-WT differential expression within cardiac cells**, followed
by GO enrichment, at each timepoint.

## Folder layout & run order

The project was reorganized (2026-06-23) into two top-level folders —
`original_Han_analysis/` (received/upstream work) and `our_analysis/` (our
re-analysis, organized into numbered step folders `00_DOCS … 06_outputs`). See
the top-level [`README.md`](../../README.md) for the full tree.

**To reproduce, follow [`RUN_ORDER.md`](RUN_ORDER.md)** — it documents every step
(environment → lane merge → condition merge → integrate/annotate → analyses →
report) with the exact command, inputs and outputs. Per-script detail is in
[`SCRIPTS.md`](SCRIPTS.md); the full process diagram is in
[`PIPELINE_FLOWCHART.md`](PIPELINE_FLOWCHART.md).

Scripts no longer run "from the project root": each anchors on the
`our_analysis/.projroot` sentinel (walking up from its own location), reads input
from `01_input/`, and writes to `processing/` and `results/`.

## Environment

- **R / Bioconductor:** Seurat, SeuratObject, glmGamPoi, dplyr, tibble, patchwork,
  dittoSeq, MAESTRO, clusterProfiler, enrichplot, org.Mm.eg.db, EnhancedVolcano, DT.
- **Python:** scrublet, scipy, numpy, pandas, matplotlib.
- The condition-merge steps load large objects; allow several GB of RAM
  (the P0WT object alone is ~300 MB on disk).

## QC assessment (first PIP-seq run for the lab)

The method choices are sound and standard — PIPseeker cell-calling → Scrublet
doublets → Seurat SCTransform normalization is a legitimate, documented PIP-seq
workflow; nothing needs to be invented. Preprocessing metrics are healthy
(confidently-mapped ~89–91%; median genes/cell ~2,700–4,500; cells recovered
P0 ~11–12k/lane, P7 ~6–8k/lane). Watch points: sequencing saturation is modest
(40–60%, lowest in P0WT_lane1 ~40%), and P7 has higher mitochondrial content
(~7–10%) than P0 (~5–7%), suggesting more stressed cells at the later timepoint.

The original R workflow had several fixable gaps, now addressed (see below):
doublets were flagged but never removed; there was no mitochondrial/upper-bound
QC (and the disabled line used the human `^MT-` pattern, not mouse `^mt-`); the
"cardiac cells" were chosen by a hardcoded UMAP region; and KO-vs-WT used a
single-cell Wilcoxon test (pseudoreplication) rather than pseudobulk.

> **Replicate structure — RESOLVED (see [`REPLICATES.md`](REPLICATES.md)):**
> `lane1`/`lane6` share the identical Novogene library ID per condition
> (e.g. `P0WT_CKDL250003755` in both `Raw_data/` and `Raw_data_Lane6/`), so they
> are the **same library on two flow-cell lanes** — **n=1 animal per condition.**
> KO-vs-WT DE is therefore **descriptive only**; valid p-values need a replicated
> cohort (≥3 animals/condition). This is why the pseudobulk DE shows thousands of
> genes at `padj≈0`.

### Reference guidance

- PIP-seq specifics: the PIPseeker / PIPseq User Guides (`docs/`).
- [Best practices for single-cell analysis across modalities (Nat Rev Genet 2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10066026/) and the Bioconductor OSCA book.
- [10x best-practices analysis guide](https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data).
- [HBC pseudobulk DESeq2 tutorial](https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html); [pseudoreplication pitfalls](https://pmc.ncbi.nlm.nih.gov/articles/PMC10695556/).

## Changes made (reorganization 2026-06-16, then QC fixes)

Reorganization:
- Flat folder reorganized into `scripts/`, `processing/`, `results/`, `docs/`;
  the 8 raw sample folders left untouched in place.
- `scRNA_P0WT.Rmd` renamed to `scRNA_lane_merge.Rmd` (reusable template).
- Hardcoded macOS paths removed (`setwd("~/Downloads/...")`, `/Users/ji51349/...`)
  → portable root-relative paths; outputs to `processing/` and `results/`.
- Marker-CSV save chunk re-enabled; final merged objects now saved.

Phase-1 QC correctness fixes:
- **Doublets are now removed** at the per-lane filter step (`scRNA_lane_merge.Rmd`).
- **Mouse mito QC added** (`^mt-`, generous `cutoff.mt = 20` because cardiomyocytes
  are mito-rich) plus a conservative per-lane upper `nFeature` cap; `percent.mt`
  added to the QC violins.
- **Lane of origin tracked** (`$lane`) so downstream pseudobulk can aggregate per replicate.
- **Cardiac-cell selection** replaced with annotation/marker-based selection
  (`Tnnt2`/`Myh6`/`Actc1`/`Nppa` fallback) in both merge scripts.
- **Pseudobulk DESeq2 path added** for KO-vs-WT (gated on ≥2 replicates/condition);
  the single-cell Wilcoxon volcano is kept but labeled *exploratory only*.
- **Cell-cycle genes** mapped via `babelgene` orthologs (title-case fallback).
- **Lane-effect DimPlot** added to flag whether integration is needed.

## Known issues / future steps (phase 2)

- **Replicate structure — RESOLVED:** n=1 per condition (technical lanes); see
  [`REPLICATES.md`](REPLICATES.md). All DE is descriptive only.
- **In progress (descriptive analyses on the pilot, written to `rerun/analysis/`):**
  4-group Harmony integration; full cell-type annotation + composition; E2F7/8
  hypothesis readouts (KO verification, target de-repression, cycling-CM fraction,
  ploidy surrogates); P7 GO parity + GSEA; cross-timepoint (P0 vs P7); scDblFinder
  doublet recall.
- **Ambient RNA correction** (SoupX or CellBender) before QC.
- **Doublet caller:** consider `scDblFinder`/`DoubletFinder` (current Scrublet
  called rate of 1–2% may be under-calling vs. the 8% expected; eyeball the
  score histograms).
- **Integration:** if the lane-effect DimPlot shows lane-driven clusters, add
  Harmony (or CCA/RPCA) for joint clustering/annotation (keep DE on RNA/pseudobulk).
- **Sensitivity:** revisit `sensitivity_4` vs `_5` per sample via knee plots.
  Note only `sensitivity_5` filtered matrices are saved locally; `sensitivity_4`
  would require re-running PIPseeker from FASTQ (not stored here).
- **Analyses not yet done:** cross-timepoint (P0 vs P7) and combined 4-group
  comparison; add the GO block to `scRNA_mergeP7.Rmd` (only P0 has it).
- The original `Mus_musculus.csv` cell-cycle list and the raw FASTQ are not here.
- New package dependencies introduced by the fixes: `DESeq2`, `babelgene`.
