#!/usr/bin/env Rscript
# Shared report CONTENT for the our_analysis outputs: the ordered `sections` list
# (narrative + figures + captions + table previews) plus the small helpers that
# build it. This is the SINGLE SOURCE OF TRUTH consumed by both:
#   * make_report.R              -> static self-contained HTML + PPTX deck
#   * interactive/report.Rmd     -> interactive self-contained HTML report
#
# Requires the path globals OUTFIG / OUTTAB (and OUT) to already be defined by the
# caller (make_report.R sources _common.R; the Rmd defines them via a light
# .projroot walk). Sourcing this file READS the result CSVs in OUTTAB, so those
# paths must exist first.
#
# All KO-vs-WT content is framed DESCRIPTIVE (n=1, sex-confounded).

stopifnot(exists("OUTFIG"), exists("OUTTAB"))

# ---- helpers ---------------------------------------------------------------
FIGDIRS <- c(OUTFIG)   # all figures consolidated under results/figures
fig <- function(name) { for (d in FIGDIRS) { p <- file.path(d, name); if (file.exists(p)) return(p) }; NA_character_ }
tab <- function(name) { p <- file.path(OUTTAB, name); if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL }
png_dims <- function(f) {            # read width/height from PNG IHDR (no deps)
  con <- file(f, "rb"); on.exit(close(con)); readBin(con, "raw", 8)
  readBin(con, "integer", 1, 4, endian = "big"); readBin(con, "raw", 4)
  c(w = readBin(con, "integer", 1, 4, endian = "big"), h = readBin(con, "integer", 1, 4, endian = "big"))
}
CONFOUND <- c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")
top_de <- function(name, n = 12) {
  d <- tab(name); if (is.null(d)) return(NULL)
  d <- d[!d$gene %in% CONFOUND & is.finite(d$log2FoldChange), ]
  d <- head(d[order(-abs(d$log2FoldChange)), c("gene","log2FoldChange","baseMean")], n)
  d$log2FoldChange <- round(d$log2FoldChange, 2); d$baseMean <- round(d$baseMean); d
}
rnd <- function(df) { if (is.null(df)) return(df); for (j in seq_along(df)) if (is.numeric(df[[j]])) df[[j]] <- round(df[[j]], 3); df }

# ---- shared content: ordered sections -------------------------------------
shared <- readLines(file.path(OUTTAB, "shared_KO_up_P0_and_P7.txt"))
shared_bio <- setdiff(shared, CONFOUND)
sections <- list(
  list(title = "E2F7/8 KO mouse-heart scRNA-seq — pilot analysis",
       narr = c("Mouse heart single-cell RNA-seq: E2F7/8 knockout vs wild-type at P0 and P7.",
                "This deck reports the corrected re-analysis and descriptive next-steps work.",
                "READ THE CAVEATS SLIDE: this is an n=1 pilot — KO-vs-WT differences are descriptive only."),
       figs = NULL, tabs = NULL),
  list(title = "Critical caveats (read first)",
       narr = c("n = 1 animal per condition: lane1/lane6 are the SAME library sequenced on two flow-cell lanes (identical CKDL IDs), not biological replicates.",
                "Processing: ~99% of barcodes are shared between lanes, so the current objects double-count nearly every cell at half depth (fix: re-count with both lanes pooled — 00_DOCS/recount_combine_lanes.sh).",
                "Sex confound: the KO and WT animals are different sexes (Y genes Eif2s3y/Kdm5d/Uty/Ddx3y top the KO-up list).",
                "KO allele: E2f7/E2f8 are NOT reduced in KO and Gt(ROSA)26Sor is KO-up -> likely a ROSA26 conditional allele invisible to 3' scRNA; confirm with the lab.",
                "=> All KO-vs-WT DE is descriptive/hypothesis-generating. Valid inference needs a replicated, sex-matched cohort (>=3 animals/condition)."),
       figs = NULL, tabs = list(list(df = tab("lane_barcode_overlap.csv"), cap = "Lane1 vs lane6 barcode overlap per library"))),
  list(title = "QC — doublets (Scrublet under-called)",
       narr = c("scDblFinder calls ~5-8% doublets vs Scrublet's ~1-2% (expected ~8%)."),
       figs = list(c(fig("QC_doublet_comparison.png"), "scDblFinder vs Scrublet per lane",
                     "Grouped bar chart: for each lane, the fraction of cells flagged as doublets by two methods. Doublets are partitions that captured two cells and masquerade as fake 'hybrid' cell types, so they must be removed before analysis. The original Scrublet calls (~1-2%) sit well below the ~8% expected for this loading, while scDblFinder recovers ~5-8% — evidence the original run under-removed doublets.")),
       tabs = list(list(df = tab("doublet_comparison.csv"), cap = "Doublet rate by lane"))),
  list(title = "Integration (Harmony, 4 groups)",
       narr = c("4-group object (58,917 cells) Harmony-integrated for joint annotation; DE kept on RNA/pseudobulk."),
       figs = list(c(fig("combined_harmony_before_after.png"), "Before vs after Harmony, colored by library",
                     "Two UMAP plots. A UMAP is a 2-D map of all cells where points close together have similar expression. Cells are colored by library (the 4 samples). LEFT (before correction) cells split by library — a technical batch effect. RIGHT (after Harmony integration) the libraries intermix while biological structure is preserved, so downstream clusters reflect cell type rather than which sample a cell came from."),
                   c(fig("combined_umap_facets.png"), "Clusters / timepoint / genotype",
                     "The same integrated UMAP shown three ways: colored by cluster, by timepoint (P0 vs P7), and by genotype (KO vs WT). Comparing the panels shows which clusters are shared across all conditions versus driven by timepoint or genotype.")),
       tabs = NULL),
  list(title = "Cell-type annotation",
       narr = c("All 28 clusters labeled (marker-based): cardiomyocyte, fibroblast, endothelial, mural/pericyte, immune/myeloid, epicardial, RBC."),
       figs = list(c(fig("annotation_umap.png"), "Cell types (UMAP)",
                     "The integrated UMAP with each cell colored by its assigned cell type. Each colored island is a distinct cell population in the heart (cardiomyocytes, fibroblasts, endothelial cells, mural/pericytes, immune/myeloid, epicardial, red blood cells)."),
                   c(fig("annotation_marker_dotplot.png"), "Canonical markers by cell type",
                     "Dot plot: known marker genes (columns) across cell types (rows). Dot SIZE = the percentage of cells in that type expressing the gene; dot COLOR = average expression level. The block-diagonal pattern (each cell type lighting up its own canonical markers, e.g. cardiomyocytes for Tnnt2/Myh6) is the evidence supporting the labels in the UMAP.")),
       tabs = NULL),
  list(title = "Composition KO vs WT (descriptive)",
       narr = c("Proportions only (no testing at n=1). P7 KO trends to more cardiomyocytes, fewer fibroblasts/pericytes."),
       figs = list(c(fig("composition_stacked.png"), "Cell-type proportions per group",
                     "Stacked bar chart: one bar per group (P0/P7 x KO/WT), each summing to 100%, with segment heights giving the fraction of each cell type. Use it to spot shifts in tissue makeup — e.g. P7 KO trends toward more cardiomyocytes and fewer fibroblasts/pericytes. With n=1 per condition this is descriptive only (no statistical test).")),
       tabs = list(list(df = tab("composition.csv"), cap = "Proportions + KO/WT log2FC"))),
  list(title = "E2F7/8 readouts (descriptive)",
       narr = c("KO verification: E2f7/E2f8 NOT reduced in KO (see caveats).",
                "Cardiomyocyte cell-cycle / E2F-target / maturation signals are modest and mixed."),
       figs = list(c(fig("e2f_P0_KO_verification.png"), "E2f7/E2f8 by genotype (P0)",
                     "Violin/box plots of E2f7 and E2f8 expression split by genotype at P0. This is meant as a knockout sanity check — you would expect the targeted genes to be LOWER in KO. They are not clearly reduced here, which is flagged in the caveats (the KO is likely a conditional allele that a 3' assay can't see)."),
                   c(fig("e2f_P7_KO_verification.png"), "E2f7/E2f8 by genotype (P7)",
                     "Same E2f7/E2f8 knockout-verification plot at P7. As at P0, the genes are not convincingly reduced in KO — see the caveats slide."),
                   c(fig("e2f_featureplot_WT.png"), "E2f7/E2f8 expression in WT (FeaturePlot, by timepoint)",
                     "FeaturePlots: the UMAP with cells shaded by E2f7 / E2f8 expression level (grey = low, color = high) in wild-type only, split by timepoint. Shows where in the tissue and at which stage these transcription factors are normally expressed."),
                   c(fig("e2f_family_dotplot.png"), "E2F family + targets by cell type (KO vs WT)",
                     "Dot plot of the broader E2F transcription-factor family plus representative E2F target genes (columns) across cell types split by genotype (rows). Dot size = % expressing, color = mean expression. Lets you scan whether E2F-target programs differ between KO and WT in any cell type.")),
       tabs = list(list(df = tab("e2f_ko_verification.csv"), cap = "E2f7/E2f8 expression"),
                   list(df = tab("e2f_cm_readouts.csv"), cap = "Cardiomyocyte readouts"))),
  list(title = "Cell cycle: proliferation & cell-cycle exit",
       narr = c("E2F7/8 drive cardiomyocyte cell-cycle EXIT; loss is predicted to keep cells cycling.",
                "P7 KO cardiomyocytes show a higher S/G2M (cycling) fraction than WT (31.6% vs 25.6%), and the cycling-CM compartment is ~1.41x larger in P7 KO (propeller). P0 shows NO difference -- the effect is specific to the P7 maturation window.",
                "Strongest hypothesis-consistent signal so far, but DESCRIPTIVE (n=1; KO=male/WT=female; transcriptional S/G2M score, not EdU/Ki67). Top priority to confirm in the replicated, sex-matched cohort."),
       figs = list(c(fig("cellcycle_fraction_by_celltype.png"), "Cycling fraction KO vs WT by cell type x timepoint",
                     "Bar chart of the fraction of actively cycling cells (in S or G2/M phase, inferred from cell-cycle gene scores) for KO vs WT, split by cell type and timepoint. The key result is in the cardiomyocytes at P7: a higher cycling fraction in KO (31.6% vs 25.6%) — consistent with E2F7/8 loss letting cardiomyocytes keep dividing instead of exiting the cell cycle. P0 shows no difference."),
                   c(fig("cellcycle_phase_composition.png"), "Cell-cycle phase composition per group",
                     "Stacked bars showing the proportion of cells in each cell-cycle phase (G1 = resting/non-dividing, S = DNA synthesis, G2/M = dividing) for each group. A larger S + G2/M slice means more proliferation."),
                   c(fig("cellcycle_ridge_P7.png"), "P7 cell-cycle score ridges by cluster",
                     "Ridge plot: stacked density curves of the S-phase and G2/M scores for each cluster at P7. Curves shifted to the right indicate clusters with more cycling cells; it shows which populations carry the proliferative signal."),
                   c(fig("cellcycle_phase_umap_P7.png"), "P7 phase UMAP (overall + by genotype)",
                     "The P7 UMAP colored by assigned cell-cycle phase, shown overall and split KO vs WT. Cycling cells (S, G2/M) cluster together; comparing the KO and WT panels shows where the extra cycling cells in KO sit.")),
       tabs = list(list(df = tab("cellcycle_cyclingCM_summary.csv"), cap = "Cycling cardiomyocytes KO vs WT"))),
  list(title = "Cardiomyocyte subtypes",
       narr = c("Cardiomyocytes subtyped by canonical panels (ventricular/atrial/trabecular/compact/cycling), matching the prior pipeline's CM subtypes.",
                "At this resolution CMs are predominantly ventricular (+ a little atrial) -- expected for the immature P0/P7 heart; the cycling state is a cross-cutting program captured by phase. Composition and cycling-by-subtype shown KO vs WT (descriptive, n=1)."),
       figs = list(c(fig("cm_subtypes_umap.png"), "CM subtypes (UMAP)",
                     "UMAP of the cardiomyocytes only, re-embedded and colored by subtype (ventricular, atrial, trabecular, compact, cycling). Shows how the cardiomyocyte compartment breaks down into finer states."),
                   c(fig("cm_subtypes_dotplot.png"), "CM subtype markers",
                     "Dot plot of the canonical marker panels used to call each cardiomyocyte subtype (dot size = % expressing, color = mean expression). Confirms each subtype label is driven by its expected markers."),
                   c(fig("cm_subtype_composition.png"), "CM subtype composition per group",
                     "Stacked bars of cardiomyocyte subtype proportions per group (P0/P7 x KO/WT). At this resolution the CMs are mostly ventricular (+ some atrial), as expected for the immature P0/P7 heart; cycling is a cross-cutting state. Descriptive (n=1).")),
       tabs = list(list(df = tab("cm_subtype_cycling.csv"), cap = "Cycling fraction by CM subtype"))),
  list(title = "Cardiomyocyte SUBCLUSTERING — marker heatmaps & resolution choice",
       narr = c("True re-clustering of the cardiomyocyte compartment (42,416 CMs re-embedded: SCTransform -> Harmony -> graph-based clustering), distinct from the marker-module CM subtypes above.",
                "Shown at two resolutions: res 0.1 (8 subclusters, clean/parsimonious) and res 0.2 (13 subclusters, which additionally splits S-phase from G2M cardiomyocytes).",
                "Reading the heatmap: clean diagonal blocks = distinct subclusters; smeared/overlapping blocks = over-clustered. At res 0.3 (not shown) blocks clearly merged, so 0.1/0.2 are kept.",
                "Caveat: a few 'subclusters' are non-CM contaminants swept in by the coarse annotation (endothelial: Car4/Cldn5/Esam; mast: Cma1/Cpa3) and should be filtered before CM-only interpretation. Descriptive (n=1, sex-confounded)."),
       figs = list(c(fig("cm_subcluster_marker_heatmap_res0.1.png"), "Subcluster markers, res 0.1 (8 subclusters)",
                     "Heatmap of the top marker genes (rows) averaged per cardiomyocyte subcluster (columns), z-scored across subclusters; the colored strip on the RIGHT labels which subcluster each block of markers defines. Sharp red blocks down the diagonal mean the 8 subclusters are transcriptionally distinct (metabolic, stress, cycling, atrial, plus endothelial/mast contaminants)."),
                   c(fig("cm_subcluster_marker_heatmap_res0.2.png"), "Subcluster markers, res 0.2 (13 subclusters)",
                     "Same marker-block heatmap at the finer resolution 0.2. It keeps the 8 groups but additionally separates S-phase from G2/M cardiomyocytes and resolves two stress states — useful for the cell-cycle-exit question, at the cost of a couple of partially-overlapping blocks.")),
       tabs = list(list(df = { d <- tab("cm_subcluster_res0.1_KOvsWT_summary.csv"); if (!is.null(d)) d[, c("subcluster","n_cells","status","n_DE_absLFC_gt1","top_KO_up")] else NULL }, cap = "res 0.1 — per-subcluster KO-vs-WT DE summary (descriptive)"),
                   list(df = { d <- tab("cm_subcluster_res0.2_KOvsWT_summary.csv"); if (!is.null(d)) d[, c("subcluster","n_cells","status","n_DE_absLFC_gt1","top_KO_up")] else NULL }, cap = "res 0.2 — per-subcluster KO-vs-WT DE summary (descriptive)"))),
  list(title = "Cardiomyocyte subclusters — cell cycle (S/G2M)",
       narr = c("Cell-cycle phase (G1/S/G2M) per cardiomyocyte subcluster, KO vs WT, at res 0.2.",
                "One subcluster is ~100% cycling and another is G2/M-dominated; these proliferative subclusters are where the failed-cell-cycle-exit signal concentrates. KO trends to a larger cycling subcluster (descriptive).",
                "Per-subcluster KO-vs-WT DE recovers the same biological candidates as the bulk-CM contrast (Gabbr2, Tcf4, Adamts9, Ralyl). Descriptive (n=1, sex-confounded)."),
       figs = list(c(fig("cm_subcluster_res0.2_cellcycle_phase_composition.png"), "Phase composition by subcluster (res 0.2)",
                     "Stacked bars of the fraction of cells in each cell-cycle phase (G1 resting, S = DNA synthesis, G2/M = dividing) for every cardiomyocyte subcluster, split KO vs WT. Tall S+G2/M stacks flag the proliferative subclusters."),
                   c(fig("cm_subcluster_res0.2_cellcycle_marker_heatmap.png"), "S/G2M phase markers by subcluster (res 0.2)",
                     "Heatmap of canonical S-phase and G2/M marker genes (rows, grouped and labeled by phase on the RIGHT) averaged per subcluster. Confirms which subclusters are cycling and whether they sit in S or G2/M."),
                   c(fig("cm_subcluster_res0.2_umap.png"), "Subcluster + phase UMAP (res 0.2)",
                     "Cardiomyocyte UMAP colored by subcluster (left) and by cell-cycle phase (right); the cycling subclusters coincide with the S/G2M-phase cells.")),
       tabs = list(list(df = tab("cm_subcluster_res0.2_cellcycle.csv"), cap = "res 0.2 — cycling fraction per subcluster KO vs WT"))),
  list(title = "Differential expression — cardiomyocytes (descriptive)",
       narr = c("KO-vs-WT in cardiomyocytes; sex/construct genes greyed; p-axis for ranking only (NOT valid, n=1).",
                paste("Top KO-up candidates (sex/construct removed):", paste(head(shared_bio, 10), collapse = ", "))),
       figs = list(c(fig("DE_P0_cardiac_volcano_MA.png"), "P0 volcano + MA",
                     "Two paired views of KO-vs-WT expression in P0 cardiomyocytes. VOLCANO (left): each dot is a gene; x = log2 fold change (right = up in KO, left = up in WT), y = statistical ranking. MA (right): x = average expression, y = log2 fold change. Sex/construct genes are greyed out. IMPORTANT: with n=1 the y-axis is for ranking candidates only — these are not valid p-values."),
                   c(fig("DE_P7_cardiac_volcano_MA.png"), "P7 volcano + MA",
                     "Same volcano + MA pair for P7 cardiomyocytes (KO vs WT). Read identically to the P0 panel; ranking only, not valid significance (n=1).")),
       tabs = list(list(df = top_de("P0.cardiac.descriptive.DE.csv"), cap = "P0 top |log2FC| (biological)"),
                   list(df = top_de("P7.cardiac.descriptive.DE.csv"), cap = "P7 top |log2FC| (biological)"))),
  list(title = "DE across all cell types + pathways",
       narr = c("Per-lineage KO-vs-WT extends beyond cardiomyocytes. Gabbr2 recurs across lineages/timepoints.",
                "GSEA (GO:BP) summarizes the P0 cardiomyocyte ranked list."),
       figs = list(c(fig("DE_percelltype_counts.png"), "DE-gene counts by cell type",
                     "Bar chart of how many genes change between KO and WT (|log2FC| > 1) in each cell type and timepoint. Taller bars = lineages where the genotype has the largest transcriptional footprint, helping prioritize where to look beyond cardiomyocytes."),
                   c(fig("GSEA_P0_cardiac_BP.png"), "P0 cardiomyocyte GSEA",
                     "Gene-set enrichment analysis (GSEA) of the P0 cardiomyocyte KO-vs-WT ranked gene list against GO Biological Process terms. Bars are biological pathways enriched among the most up- or down-shifted genes — a pathway-level summary rather than single genes. Descriptive at n=1."),
                   c(fig("candidates_Gabbr2_dotplot.png"), "Candidate genes by cell type/genotype",
                     "Dot plot tracking the recurring candidate genes (e.g. Gabbr2) across cell types split by genotype (dot size = % expressing, color = mean expression). Shows whether a candidate's KO-vs-WT shift is cardiomyocyte-specific or shared across lineages.")),
       tabs = list(list(df = tab("percelltype_KOvsWT_summary.csv")[, c("timepoint","celltype","n_cells","n_DE_absLFC_gt1")], cap = "Per-cell-type DE summary"))),
  list(title = "Trajectory: cardiomyocyte maturation pseudotime",
       narr = c("Slingshot orders cardiomyocytes along an immature->mature axis (validated: pseudotime correlates +0.55 with maturation markers, negative with cycling).",
                "P7 cells are more mature than P0 (expected), but KO and WT pseudotime distributions overlap within each timepoint -- no clear maturation delay in KO along this axis. Descriptive (n=1, sex-confounded).",
                "RNA velocity was not possible here (no spliced/unspliced counts; would need the cluster BAMs)."),
       figs = list(c(fig("traj_umap_pseudotime.png"), "CM maturation pseudotime (UMAP)",
                     "Cardiomyocyte UMAP colored by pseudotime — a computed ordering of cells along an inferred immature-to-mature trajectory (Slingshot). It is a relative ordering, not real clock time: darker/lighter shading = earlier/later along the maturation axis."),
                   c(fig("traj_pseudotime_KOvsWT.png"), "Pseudotime KO vs WT by timepoint",
                     "Distributions (density/violin) of cardiomyocyte pseudotime for KO vs WT, separated by timepoint. P7 cells sit later (more mature) than P0 as expected, but KO and WT largely overlap within each timepoint — i.e. no obvious maturation delay in KO along this axis. Descriptive (n=1)."),
                   c(fig("traj_marker_validation.png"), "Direction check vs maturation markers",
                     "Validation that pseudotime points the right way: pseudotime correlates positively (~+0.55) with maturation markers and negatively with cycling genes. This confirms 'higher pseudotime = more mature' before interpreting the KO-vs-WT comparison.")),
       tabs = list(list(df = tab("pseudotime_KOvsWT.csv"), cap = "CM pseudotime summary KO vs WT"))),
  list(title = "TF regulon activity (E2F)",
       narr = c("decoupleR over MSigDB E2F-target regulons (robust substitute for full SCENIC, which needs ~1GB cisTarget databases).",
                "E2F-target regulon activity is modestly HIGHER in KO cardiomyocytes (E2F2/E2F5 targets) -- consistent with de-repression on loss of the repressor E2F7/8; mixed in other cell types. Descriptive (n=1, sex-confounded)."),
       figs = list(c(fig("tf_activity_E2F_heatmap.png"), "E2F-family regulon activity by cell type/genotype",
                     "Heatmap of inferred E2F-family transcription-factor ACTIVITY (estimated from how strongly each TF's target genes are expressed, via decoupleR) across cell types and genotype. Color = activity (red higher / blue lower). E2F-target activity is modestly higher in KO cardiomyocytes — consistent with de-repression when the repressors E2F7/8 are lost. Descriptive (n=1)."),
                   c(fig("tf_activity_top_KOvsWT.png"), "Top cardiomyocyte TF activity differences (KO-WT)",
                     "Bar chart of the transcription factors whose inferred activity differs most between KO and WT cardiomyocytes (KO minus WT). Bars to the right = more active in KO. A regulon-level complement to the gene-by-gene DE.")),
       tabs = list(list(df = { d <- tab("e2f_regulon_activity.csv"); if (!is.null(d)) d[d$celltype == "Cardiomyocyte", ] else NULL },
                        cap = "E2F regulon activity, cardiomyocytes (KO vs WT)"))),
  list(title = "Cell-cell communication (CellChat)",
       narr = c("CellChat ligand-receptor signaling, KO vs WT per timepoint (descriptive; aggregates over cells, no replicate stats at n=1).",
                "Differential interaction maps (KO - WT) and the most-changed signaling pathways."),
       figs = list(c(fig("cellchat_P0_diff_heatmap.png"), "P0 differential interactions (KO - WT)",
                     "Heatmap of the DIFFERENCE in cell-cell signaling between KO and WT at P0 (KO minus WT), inferred by CellChat from ligand-receptor co-expression. Rows = signal-sending cell types, columns = receiving cell types; red = more signaling in KO, blue = more in WT. Highlights which communication routes are rewired by genotype. Descriptive — no replicate statistics at n=1."),
                   c(fig("cellchat_P7_diff_heatmap.png"), "P7 differential interactions (KO - WT)",
                     "Same KO-minus-WT differential interaction heatmap at P7. Read the same way (rows = sender, columns = receiver; red = up in KO). Compare against P0 to see signaling changes specific to the later timepoint.")),
       tabs = list(list(df = { d <- tab("cellchat_P7_pathway_diff.csv"); if (!is.null(d)) head(d, 12) else NULL },
                        cap = "P7 signaling pathways changed KO vs WT (top by |diff|)"))),
  list(title = "Can we answer the biological question?",
       narr = c("Question: does E2F7/8 loss disrupt cardiomyocyte cell-cycle exit, maturation, and binucleation (P0->P7)?",
                "Cell-cycle EXIT: addressable (descriptive) -- P7 KO retains an expanded cycling-CM compartment; strongest lead.",
                "MATURATION: addressable (descriptive) -- pseudotime + maturation scores; KO-vs-WT approximately null along the axis.",
                "BINUCLEATION: NOT measurable by 3' scRNA (a nuclei/ploidy phenotype) -- needs imaging or flow cytometry; neither pipeline can address it.",
                "All gated by n=1 + sex confound + unconfirmed KO -> a replicated, sex-matched cohort + KO confirmation + a binucleation assay are required. Full detail in BIOLOGICAL_QUESTION.md."),
       figs = NULL, tabs = NULL),
  list(title = "Conclusions & next steps",
       narr = c("Pilot conclusions are descriptive: candidate KO-up genes Gabbr2, Tcf4, Adamts9 (consistent across P0/P7).",
                "Top biological lead: P7 KO cardiomyocytes retain an expanded cycling fraction (failed cell-cycle exit) -- confirm with EdU/Ki67 in the replicated cohort.",
                "1) Re-count with lanes pooled (00_DOCS/recount_combine_lanes.sh) for accurate per-cell metrics.",
                "2) Confirm the KO allele / Cre driver and animal sexes with the lab.",
                "3) Run a replicated, sex-matched cohort (>=3 animals/condition) for valid KO-vs-WT statistics.",
                "Full detail: REPLICATES.md; tables in results/tables/."),
       figs = NULL, tabs = NULL)
)

# Stable keys so consumers (interactive/report.Rmd) can reference sections by name
# instead of fragile positional index. Order is unchanged, so make_report.R's
# positional `for (i in seq_along(sections))` loop is unaffected.
names(sections) <- c(
  "title", "caveats", "qc", "integration", "annotation", "composition",
  "e2f", "cellcycle", "cm_subtypes", "cm_subcluster", "cm_subcluster_cc",
  "de_cardiac", "de_all", "trajectory", "tf", "cellchat", "biological_q",
  "conclusions"
)
