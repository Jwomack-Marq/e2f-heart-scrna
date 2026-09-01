# original_Han_analysis — work we did NOT do

This folder isolates the **upstream and original analysis we received**, kept
separate from our re-analysis (which lives in `../our_analysis/`). Nothing here
was authored by us; it is preserved as-is for provenance and for the
`../our_analysis/06_outputs/compare_results.R` comparison.

## Contents

| Folder | What it is |
|---|---|
| `raw_pipseeker_outputs/` | The 8 raw PIPseeker v3.3.0 sample folders (`P0KO_lane1 … P7WT_lane6`), produced on the compute cluster: STAR alignment to mouse GRCm38 → cell-called `filtered_matrix/sensitivity_5/` count matrices, plus `metrics/`, `star/`, and the per-lane Scrublet output (`predicted_doublets.csv`). **The data everything started from.** (Raw FASTQ stayed on the cluster.) |
| `preprocessing_scripts/` | Reference commands for the upstream steps: `pipseeker_scripts.txt` (PIPseeker cell-calling, all 8 libraries) and `scrublets.py` (Scrublet doublet detection). Documentation of already-completed work — not re-run. |
| `processing/` | The original run's intermediate Seurat objects (`merge.lanes.{P0KO,P0WT,P7KO,P7WT}.rds`). |
| `results/` | The original run's outputs: `figures/`, `markers/`, and rendered `reports/` (`.html`). |
| `docs/` | `PIPseeker-v3.3-User-Guide.pdf`. |

## Important note on the analysis scripts
The original R pipeline scripts (`scRNA_lane_merge.Rmd`, `scRNA_mergeP0.Rmd`,
`scRNA_mergeP7.Rmd`) were **edited by us** (QC correctness fixes — doublet removal,
mouse mito QC, marker-based cardiac selection, pseudobulk DE, portable paths) and
therefore now live with our re-analysis at
`../our_analysis/02_qc_lane_merge/` and `../our_analysis/03_condition_merge/`.
The `processing/` and `results/` here are the outputs of the **original,
un-patched** versions of those scripts. See
`../our_analysis/00_DOCS/METHODS_COMPARISON.md` for exactly what changed and why.

> A copy of `raw_pipseeker_outputs/` also lives at `../our_analysis/01_input/`
> because our pipeline reads it as input; the two copies are identical.
