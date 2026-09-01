# Methods comparison — our reanalysis vs. the prior MCW analysis

Compares our work (`rerun/` + `rerun/analysis/`) against the earlier analysis at
`…/OneDrive - mcw…/Projects/E2F 7_8/2025_scRNA-seq/`. Purpose: confirm we cover everything
they did, and document how our **normalization** (and methods generally) differ, plus what we added.

## The two things in their folder

1. **Original analyst's executed pipeline** — rendered `scRNA_*.html` reports + `*.allmarkers.csv`.
   Same Seurat-Rmd approach our project's originals come from: per-lane → `merge()` (no batch
   correction) → SCTransform → cluster → single-cell-Wilcoxon `FindAllMarkers` → MAESTRO annotation.
2. **A newer Quarto "upgrade" pipeline** `analysis/01_load_qc…08_e2f_atlas.qmd` (+ `_env.R`, a
   `00_tierA_csv_audit.qmd`). This is the substantive comparison target. **Important:** it is largely a
   *scaffold* — QC thresholds are explicitly `# placeholder`, several steps are "pending arrival of
   `.rds` files," and only the tier-A CSV audit (00/00b) actually ran (producing the `tierA_plots/`).
   Its 02–08 steps therefore describe their **intended** methods, mostly not executed.

## Side-by-side methods

| Stage | Prior MCW (Quarto intended) | Ours (executed) | Same? |
|---|---|---|---|
| Input matrices | `filtered_feature_bc_matrix` per sample; **no lane handling** (left as a TODO in `_env.R`) | PIPseeker `filtered_matrix/sensitivity_5`; lane1+lane6 `merge()`d (inherited) — then **we found this double-counts** | We went further |
| QC cutoffs | nFeature 500–6000, percent.mt ≤10 (**placeholders**) | nFeature ≥1500 & ≤99.5th pctile, percent.mt ≤20 | Differ (cutoffs) |
| Doublets | `scDblFinder`, removed | Scrublet (removed in filter) **+ scDblFinder recall** vs Scrublet | Both; we cross-checked |
| **Normalization** | `SCTransform(vars.to.regress="percent.mt")`, `SelectIntegrationFeatures(nfeatures=3000)` | `SCTransform(method="glmGamPoi")`, **no** percent.mt regression | **Differ — see below** |
| Integration | **Harmony** on `sample` (4 jointly), res 0.6 | per-condition `merge()` **+** 4-group **Harmony** on `orig.ident`, res 0.8 | Both Harmony ✓ |
| Annotation | SingleR `celldex::MouseRNAseqData` `label.fine` + marker curation + CellCycleScoring | marker-module argmax + CellCycleScoring; **SingleR attempted, celldex won't build** | We matched except SingleR (blocked) |
| Cluster markers | `FindAllMarkers` (Wilcoxon), min.pct 0.25, logfc 0.25 | `FindAllMarkers` (presto) per condition | Same ✓ |
| KO-vs-WT DE | DESeq2 pseudobulk `~genotype` by `celltype×sample`; **n<2 warning**; `if (sum(keep) < 4) return(NULL)` → **skips everything at n=1** → no DE output; **no lfcShrink** | DESeq2 pseudobulk `~condition`, lanes as 2 pseudoreps → runs; **apeglm `lfcShrink`**; cardiac + per-cell-type + cross-timepoint; labeled descriptive | Same tool; we produce a descriptive ranking, they produce none |
| Pathway | `fgsea` over **Hallmark + KEGG_LEGACY + E2F-target(TFT:GTRD)** + `enrichGO` BP | `gseGO` BP + `enrichGO` BP — **now + `fgsea` Hallmark/KEGG/E2F-TFT** (`pathway_msigdb.R`) | Now matched ✓ |
| Abundance | `speckle::propeller` per timepoint | descriptive proportions — **now + `propeller`** (`abundance_propeller.R`) | Now matched ✓ |
| E2F / KO check | E2f7/8 FeaturePlot (WT) + VlnPlot (KO vs WT) + E2F-target module score; tier-A flags positive E2f7/8 in KO | KO verification + E2F-target module + **cycling-CM fraction** + **ploidy surrogates** + concluded the KO-not-reduced cause | We went further |
| Replicate structure | **Unresolved TODO**; pseudobulk warns "p-values not meaningful" | **Resolved**: identical CKDL library IDs + 99% lane barcode overlap → n=1; `recount_combine_lanes.sh` | We resolved |
| Sex | **Never checked** | `sex_check.R`: KO=male, WT=female (Xist + Y genes) | We added |

## Normalization (the headline difference)

**Theirs** (`02_normalize_integrate.qmd`):
```r
SCTransform(obj, vars.to.regress = "percent.mt")
features <- SelectIntegrationFeatures(samples, nfeatures = 3000)
RunHarmony(merged, group.by.vars = "sample", assay.use = "SCT", dims.use = 1:30)
FindClusters(resolution = 0.6)
```
**Ours** (`scRNA_*merge*.Rmd`, `rerun/analysis/combined.R`):
```r
SCTransform(obj, method = "glmGamPoi", conserve.memory = TRUE)   # no vars.to.regress
RunHarmony(comb, group.by.vars = "orig.ident")                    # combined object
FindClusters(resolution = 0.8)
```

Both use **SCTransform → Harmony**, so the integration *engine* is the same. Three differences:

1. **`percent.mt` regression — they do, we don't.** This is the most consequential normalization
   choice. Regressing `percent.mt` removes mito-content as a covariate. In heart tissue this is a
   double-edged sword: cardiomyocytes are *genuinely* mito-rich, so a large part of `percent.mt`
   variance is real CM biology, not a technical artifact — regressing it can blunt true cardiomyocyte
   signal. We deliberately omitted it (and instead capped `percent.mt ≤ 20`, generous on purpose).
   Their choice is defensible too (removes a stress/quality gradient); it's a judgment call, not a
   right/wrong. **At n=1 it does not change any biological conclusion.**
2. **glmGamPoi backend (ours) vs default (theirs).** `method="glmGamPoi"` is a faster, numerically
   robust negative-binomial fit; results are essentially equivalent to the default for downstream use.
3. **Feature selection & resolution.** They fix 3000 integration features and cluster at res 0.6;
   we use SCTransform's default variable features and res 0.8 (→ a few more clusters). Minor, and
   re-clustering granularity doesn't affect the pseudobulk/descriptive DE.

Per the agreed scope this is an **analytical** comparison (no head-to-head re-run); given n=1 the
gene-level descriptive rankings are robust to these choices.

## Coverage check — do we now do everything they did?

| Their analysis | Our status |
|---|---|
| QC + doublet removal | ✓ (scDblFinder also cross-checked) |
| SCTransform + Harmony integration | ✓ (`combined.R`) |
| Cell-type annotation | ✓ marker-based; **SingleR blocked** (celldex/ExperimentHub won't install — `singler_blocked.txt`; marker-based is the heart-appropriate fallback) |
| Cluster markers | ✓ |
| Pseudobulk KO-vs-WT DE | ✓ and beyond (they produce none at n=1; we give a descriptive lfcShrink ranking) |
| fgsea Hallmark/KEGG/E2F-target | ✓ **added** (`pathway_msigdb.R`; 303 terms at padj<0.1) |
| enrichGO BP ORA | ✓ |
| propeller abundance | ✓ **added** (`abundance_propeller.R`) |
| E2F atlas / KO verification | ✓ and beyond |

Residual gaps are cosmetic only (their tier-A `specificity_scatter`, `top5_heatmap` — visualization
audits, not analyses). The one true gap is **SingleR**, blocked by celldex's dependency stack.

## What we did in addition (not in their analysis)

- **Resolved the replicate structure** they left open: matching `CKDL` library IDs + 97–100% lane
  barcode overlap → **n=1 per condition**; wrote the standard re-count fix (`recount_combine_lanes.sh`).
- **Sex confound** (they never checked): `KO=male, WT=female` from Xist + Y genes — genotype is
  perfectly confounded with sex.
- **Concluded the KO-not-reduced finding** their own tier-A data showed but didn't interpret
  (autorepression / likely ROSA26 conditional allele; `Gt(ROSA)26Sor` KO-up).
- **apeglm lfcShrink** effect-size ranking (they used raw LFC); **cross-timepoint** P0-vs-P7 DE;
  **cycling-CM fraction** and **ploidy/maturation surrogates**; explicit **descriptive(n=1)** framing
  + `REPLICATES.md` + replicated-cohort spec; scDblFinder-vs-Scrublet comparison; self-contained
  **HTML report + PPTX deck**.

## Cross-pipeline agreement (reassuring)

- Both independently flagged the **n=1 replicate problem** (they as a warning/TODO, we resolved it).
- Their tier-A audit and our DE **both show E2f7/E2f8 NOT depleted in KO**.
- The newly-added **propeller** reproduces our descriptive composition shift (**P7 KO: more
  cardiomyocytes, fewer fibroblasts/pericytes**), and the newly-added **Hallmark E2F-target fgsea**
  quantifies **E2F-target de-repression in KO** (positive NES in CM/fibroblast at P7, CM at P0) — the
  predicted direction for loss of repressor E2Fs. All descriptive (n=1, sex-confounded).

Outputs: `rerun/analysis/tables/{pathway_fgsea_*, propeller_results, singler_blocked}.csv/.txt`,
`rerun/analysis/figures/{pathway_fgsea_dotplot, abundance_propeller}.png`.
