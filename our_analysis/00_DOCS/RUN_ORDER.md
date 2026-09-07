# RUN ORDER — how to reproduce our re-analysis

This folder (`our_analysis/`) is organized into **numbered steps** that match the
pipeline stages. Run them top to bottom. Each step's scripts carry a header
comment describing exactly what they read and write; this file is the map.

> **Read the caveats first** (`00_DOCS/REPLICATES.md`, and the "Critical caveats"
> slide in the report): this is an **n = 1 animal per condition** pilot
> (lane1/lane6 are the *same* library on two flow-cell lanes), **genotype is
> confounded with sex** (KO male / WT female), and the KO is likely a ROSA26
> conditional allele invisible to 3′ scRNA. **All KO-vs-WT differential
> expression is descriptive / hypothesis-generating only — there are no valid
> p-values.**

## How paths work (so this runs anywhere)
Every script anchors itself on the `our_analysis/.projroot` sentinel file by
walking up from its own location, so you can launch scripts from any working
directory. Inputs come from `01_input/`, intermediate objects from `processing/`,
and outputs go to `results/` (`figures/`, `markers/`, `tables/`) and `06_outputs/`.
Set `$rs` to your Rscript once:

```powershell
$rs = "C:\Program Files\R\R-4.5.2\bin\Rscript.exe"   # adjust to your R install
```

## Step 0 — environment (one time)
Install the R/Bioconductor stack and **`presto`** (prevents the FindAllMarkers
memory blow-up). Full instructions: [`HANDOFF_setup_and_run.md`](HANDOFF_setup_and_run.md).
Helper install scripts live in [`../_setup/`](../_setup/). Run **one sample at a
time** on a machine with ≥25 GB free RAM; the runner forces a sequential `future`
plan with a 3 GiB brake (never `Inf`).

Upstream preprocessing (PIPseeker v3.3.0 → count matrices, then Scrublet doublets)
was **done by others on the cluster** — we received its outputs. Those reference
scripts are in [`../../original_Han_analysis/preprocessing_scripts/`](../../original_Han_analysis/preprocessing_scripts/);
they are **not re-run** here (the raw FASTQ are not stored locally).

## Step 01 — input
`01_input/` holds the 8 raw PIPseeker sample folders
(`P0KO_lane1 … P7WT_lane6`, each with `filtered_matrix/sensitivity_5/` +
`predicted_doublets.csv`). This is the copy our pipeline reads. Confirm all 8 are
present and fully downloaded (not OneDrive cloud placeholders) before starting.

## Step 02 — QC + lane merge  (`02_qc_lane_merge/`)
Per condition: load both lanes → QC filter (nFeature ≥1500 & ≤99.5pctl,
percent.mt ≤20, drop Scrublet doublets) → merge lanes → SCTransform → cluster →
UMAP → markers. Run all four:

```powershell
& $rs ..\run_pipeline.R lane P0WT
& $rs ..\run_pipeline.R lane P0KO
& $rs ..\run_pipeline.R lane P7WT
& $rs ..\run_pipeline.R lane P7KO
```
→ writes `processing/merge.lanes.<group>.rds`, `results/figures/<group>cluster.umap.png`,
`results/markers/<group>.{allmarkers,cluster.top30.markers}.csv`.

## Step 03 — condition merge + pseudobulk DE  (`03_condition_merge/`)
Combine WT+KO per timepoint → cluster/annotate → cardiac selection →
pseudobulk DESeq2 (descriptive). Needs the 4 objects from Step 02.

```powershell
& $rs ..\run_pipeline.R merge P0
& $rs ..\run_pipeline.R merge P7
```
→ writes `processing/seurat.{P0,P7}.merge.rds`,
`results/markers/{P0,P7}.cardiac.pseudobulk.DESeq2.csv`.

## Step 04 — integrate + annotate  (`04_integrate_annotate/`)
Builds the shared objects the Step 05 analyses depend on. Run **in this order**:

```powershell
& $rs combined.R            # 4-group Harmony embedding   -> processing/seurat.combined.rds
& $rs annotate.R            # marker-based cell types     -> processing/seurat.combined.annotated.rds
& $rs annotate_singler.R    # optional SingleR cross-check (skips if celldex absent)
& $rs cm_subcluster_build.R # cardiomyocyte subclustering -> processing/seurat.cm.subclustered.rds
```

## Step 05 — analyses  (`05_analyses/`)
Each reads an object from Step 03/04 and writes to `results/{tables,figures}/`.
They are **independent of each other** — run any/all, in any order (one isolated
Rscript process each). See each file's header for its exact inputs/outputs.

| Script | What it does |
|---|---|
| `de_descriptive.R` | Cardiomyocyte KO-vs-WT pseudobulk DE (descriptive) |
| `per_celltype_de.R` | KO-vs-WT DE within each annotated cell type |
| `cross_timepoint.R` | P0↔P7 cardiac DE within genotype |
| `pathway_msigdb.R` | GO/GSEA + fgsea over Hallmark/KEGG/E2F sets |
| `gene_set_provenance.R` | Audit, generated from the code, of where every gene set came from. The headline it exists to surface rather than bury: 38 of 46 sets were typed into a script with no citation — not wrong, but unauditable |
| `gene_set_benchmark.R` | Follow-up to the audit: benchmarks each hand-curated panel against MSigDB Hallmark, since "are they right" has no ground truth but "do they agree with a published set" does |
| `e2f_readouts.R` | KO verification, E2F-target de-repression, cycling/maturation |
| `e2f_atlas.R` | E2F family expression atlas plots |
| `abundance_propeller.R` | Cell-type composition + propeller |
| `cell_cycle.R` | Phase composition; cycling fraction by cell type |
| `cellcycle_ridge.R` | Cell-cycle score ridge plots |
| `cellcycle_tricycle.R` | Independent cross-check of the phase calls: projects cells onto a fixed external reference cycle for a continuous position θ, plus depth and ambient-RNA controls. Agrees with `CellCycleScoring` (κ 0.87); both are depth-floored. Runs in `e2f-tricycle` (`Dockerfile.tricycle`) |
| `cellcycle_tricycle_depthmatched.R` | Re-asks KO-vs-WT after binomially thinning both genotypes to one depth distribution per cell type × timepoint, then re-projecting through tricycle. At P7 the WT is the deeper library, so the raw gap *understates*: P7 CM KO−WT goes +3.4 → +4.1 points. Runs in `e2f-tricycle` |
| `cm_cycling_investigate.R` | Advisor follow-up: are the two "cycling" CM subclusters two identities or one proliferating population split by phase? Marker Jaccard + phase composition per subcluster. Answer: one population — CM4 is 24 % S / 76 % G2M with no G1 |
| `cm_subtypes.R` | Cardiomyocyte subtype labeling/composition |
| `cm_immune_contamination.R` | Flags the ~99 leukocytes (about half mast cells) that `annotate.R` labelled "Cardiomyocyte" via `max.col` over nine marker scores. Not doublets — Scrublet clears them all. Feeds `shiny_app/build_immune_flag.R` |
| `cm_subcluster_analyze.R` | KO-vs-WT DE per CM subcluster (resolution sweep) |
| `cm_subcluster_celllevel_de.R` | Cell-level Wilcoxon KO-vs-WT for the CM subclusters pseudobulk cannot test (CM12 at res 0.2 has 99 cells, so only one sample clears the 20-cell floor). A **ranking by AUC, not a test** — feeds `shiny_app/build_celllevel_de.R` |
| `trajectory_slingshot.R` | Pseudotime / CM maturation ordering |
| `cell_state_classifier.R` | Portable glmnet cell-type + CM-stage (P0/P7) predictors; marker panels, held-out accuracy, `predict_cell_state()` applier (descriptive; genotype NOT predicted) |
| `baniol_maturation_validation.R` | The one **external** test this dataset supports: applies the CM maturation (P0-vs-P7) model to Baniol et al. 2021 (285 FACS-sorted CMs, Smart-seq2, ENA PRJEB47622). Validates the timepoint model only — genotype is not testable there |
| `tf_activity.R` | decoupleR TF-regulon activity over E2F regulons |
| `cellchat.R` | Cell–cell communication, KO vs WT per timepoint |
| `sex_check.R` | Xist / Y-gene sex calls per sample-lane |
| `qc_doublets.R` | scDblFinder vs Scrublet doublet recall (reads `01_input/`) |
| `barcode_overlap.R` | lane1∩lane6 barcode overlap (reads `01_input/`) |
| `marker_qc_plots.R` | Marker-QC plots from `results/markers/` |

## Step 06 — outputs / comparison  (`06_outputs/`)
```powershell
& $rs make_result_plots.R   # assembles result figures from results/tables
& $rs make_report.R         # self-contained HTML + PPTX deck (E2F_scRNA_pilot_report.*)
& $rs compare_results.R     # our outputs vs ../original_Han_analysis (cells, clusters, markers)
& $rs qc_compare.R          # QC filter impact per lane (reads 01_input) -> qc_comparison.csv
```

## Dependency summary
```
01_input ─▶ 02 (lane merges) ─▶ 03 (condition merges) ─┐
                                  └▶ 04 combined ─▶ annotate ─▶ cm_subcluster_build
                                                           └────────────┬─────────────┘
                                                                        ▼
                                                                  05 analyses ─▶ 06 outputs
```
Full process diagram: [`PIPELINE_FLOWCHART.md`](PIPELINE_FLOWCHART.md).
Per-script detail: [`SCRIPTS.md`](SCRIPTS.md). Methods rationale & fixes vs the
original run: [`METHODS_COMPARISON.md`](METHODS_COMPARISON.md), [`NORMALIZATION.md`](NORMALIZATION.md).
