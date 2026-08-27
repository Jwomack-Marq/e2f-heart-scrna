#!/usr/bin/env Rscript
# build_linked_results.R
# ---------------------------------------------------------------------------
# Link every result the upstream pipeline already COMPUTED into the app bundle, so
# the browser reads them instead of recomputing, approximating, or simply not
# offering them.
#
#   Rscript shiny_app/build_linked_results.R --probe   # inventory + sizes, writes nothing
#   Rscript shiny_app/build_linked_results.R           # ingest + save
#
# WHY THIS EXISTS. our_analysis/results/ holds ~100 computed tables. The bundle
# carried twelve of them. Everything else was either invisible in the app or, worse,
# re-approximated inside it -- the Cell-cell signalling tab scores 20 curated ligand
# pairs by hand while a full CellChat run with permutation tests sat on disk unused;
# the Composition tab shows proportions with no test while propeller results sat on
# disk unused. Linking beats recomputing on both correctness and load.
#
# It also removes a real data-completeness problem. build_app_data.R trims stored DE
# tables at baseMean >= 3 (DE_MIN_BASEMEAN), which drops ~4 % of genes -- the lowest
# expressed, but a gene a reader looks up and cannot find reads as "not measured".
# This script re-reads the complete CSVs and keeps every row.
#
# Two things it deliberately does NOT ingest:
#   * results/markers/*.allmarkers.csv (64 MB) -- the per-cluster top markers already
#     in the bundle answer the same question at 1/100th the size.
#   * anything requiring the multi-GB Seurat objects. This script only reads CSVs.
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))

ARGS  <- commandArgs(trailingOnly = TRUE)
PROBE <- "--probe" %in% ARGS
argval <- function(f, d = "") { h <- grep(paste0("^", f, "="), ARGS, value = TRUE)
                                if (length(h)) sub(paste0("^", f, "="), "", h[1]) else d }

if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

OUR <- argval("--our", "../our_analysis")
TAB <- file.path(OUR, "results", "tables")
MRK <- file.path(OUR, "results", "markers")
stopifnot("our_analysis/results/tables not found -- pass --our=<path>" = dir.exists(TAB))

# ---- what to link, and what each thing is ---------------------------------
# `key` is the slot name the app uses; `file` is relative to results/; `label` and
# `note` are what the app shows, because a table browser full of raw filenames is not
# a linked result, it is a directory listing.
LINK <- function(key, file, label, note, group)
  list(key = key, file = file, label = label, note = note, group = group)

WANT <- list(
  # ---- the gap that matters most: full-data temporal contrasts -------------
  LINK("crosstime_WT_P7_vs_P0", "tables/crosstime_WT_P7vsP0.cardiac.descriptive.DE.csv",
       "WT: P7 vs P0 — all cardiomyocytes, full data",
       paste("The wild-type maturation contrast computed on the FULL object, every gene,",
             "no row gate. The Four-group tab's version of this contrast is computed from",
             "the downsampled bundle and keeps ~3 % of genes; this is the complete answer."),
       "Differential expression (full data)"),
  LINK("crosstime_KO_P7_vs_P0", "tables/crosstime_KO_P7vsP0.cardiac.descriptive.DE.csv",
       "KO: P7 vs P0 — all cardiomyocytes, full data",
       "As above, within the knockout.",
       "Differential expression (full data)"),
  LINK("pseudobulk_P0", "markers/P0.cardiac.pseudobulk.DESeq2.csv",
       "P0 cardiomyocytes — pseudobulk DESeq2",
       paste("DESeq2 on lane-level pseudobulk rather than cell-level Wilcoxon. With two lanes",
             "per library these p-values are still not biologically replicated, but the test",
             "does not pretend cells are independent, so the ranking is less inflated."),
       "Differential expression (full data)"),
  LINK("pseudobulk_P7", "markers/P7.cardiac.pseudobulk.DESeq2.csv",
       "P7 cardiomyocytes — pseudobulk DESeq2", "As above, at P7.",
       "Differential expression (full data)"),

  # ---- cell-cell communication: the real CellChat run ---------------------
  LINK("cellchat_P0_pathway", "tables/cellchat_P0_pathway_diff.csv",
       "CellChat P0 — pathway-level KO vs WT",
       paste("A full CellChat run with permutation testing. The app's Cell-cell signalling",
             "tab scores 20 curated ligand-receptor pairs by hand with no significance test;",
             "this is the database-driven, tested version of the same question."),
       "Cell-cell communication"),
  LINK("cellchat_P7_pathway", "tables/cellchat_P7_pathway_diff.csv",
       "CellChat P7 — pathway-level KO vs WT", "As above, at P7.", "Cell-cell communication"),
  LINK("cellchat_P0_weight", "tables/cellchat_P0_diff_weight.csv",
       "CellChat P0 — interaction strength, KO − WT", "Sender x receiver interaction weight.",
       "Cell-cell communication"),
  LINK("cellchat_P7_weight", "tables/cellchat_P7_diff_weight.csv",
       "CellChat P7 — interaction strength, KO − WT", "Sender x receiver interaction weight.",
       "Cell-cell communication"),
  LINK("cellchat_P0_count", "tables/cellchat_P0_diff_count.csv",
       "CellChat P0 — interaction count, KO − WT", "Number of inferred interactions.",
       "Cell-cell communication"),
  LINK("cellchat_P7_count", "tables/cellchat_P7_diff_count.csv",
       "CellChat P7 — interaction count, KO − WT", "Number of inferred interactions.",
       "Cell-cell communication"),

  # ---- composition, with a test ------------------------------------------
  LINK("propeller", "tables/propeller_results.csv",
       "Cell-type abundance — propeller test",
       paste("propeller tests whether cell-type proportions differ, on the logit scale with",
             "an empirical-Bayes moderated statistic. The Composition tab draws proportions",
             "with no test at all; this is the test."),
       "Composition and cell cycle"),
  LINK("propeller_cellcycle", "tables/cellcycle_propeller_CMsplit.csv",
       "Cell-cycle phase abundance — propeller test",
       "The same test applied to phase composition, split within cardiomyocytes.",
       "Composition and cell cycle"),
  LINK("cellcycle_phase_composition", "tables/cellcycle_phase_composition.csv",
       "Cell-cycle phase composition", "Phase fractions per group.", "Composition and cell cycle"),
  LINK("cellcycle_cycling_cm", "tables/cellcycle_cyclingCM_summary.csv",
       "Cycling cardiomyocyte summary", "Cycling fraction within the CM compartment.",
       "Composition and cell cycle"),
  LINK("cm_subtype_composition", "tables/cm_subtype_composition.csv",
       "CM subtype composition", "Subtype fractions per group.", "Composition and cell cycle"),
  LINK("cm_subtype_cycling", "tables/cm_subtype_cycling.csv",
       "CM subtype cycling fractions", "Cycling fraction per CM subtype.",
       "Composition and cell cycle"),

  # ---- the confounds, quantified upstream --------------------------------
  LINK("sex_calls", "tables/sex_calls.csv",
       "Sex call per sample",
       paste("The sex confound as the upstream pipeline measured it. KO and WT animals are",
             "different sexes; this is the evidence, not the assertion."),
       "Confounds and QC"),
  LINK("sex_markers", "tables/sex_markers_by_sample_lane.csv",
       "Sex-marker expression per sample and lane", "Y-linked and Xist expression per lane.",
       "Confounds and QC"),
  LINK("e2f_ko_verification", "tables/e2f_ko_verification.csv",
       "E2f7/E2f8 knockout verification",
       paste("The upstream check on whether the knockout is visible in the transcript.",
             "It is not -- see the E2F focus tab and the methods book."),
       "Confounds and QC"),
  LINK("e2f_cm_readouts", "tables/e2f_cm_readouts.csv",
       "E2F cardiomyocyte readouts", "Per-group E2F-axis readouts.", "Confounds and QC"),
  LINK("lane_barcode_overlap", "tables/lane_barcode_overlap.csv",
       "Lane barcode overlap",
       paste("How far the two lanes of each library share barcodes -- the check that they are",
             "the same library sequenced twice rather than independent samples."),
       "Confounds and QC"),
  LINK("doublet_comparison", "tables/doublet_comparison.csv",
       "Doublet calls — Scrublet vs scDblFinder", "Per-lane doublet rates from both callers.",
       "Confounds and QC"),

  # ---- annotation provenance ---------------------------------------------
  LINK("annotation_cluster_labels", "tables/annotation_cluster_labels.csv",
       "Cluster to cell-type mapping",
       paste("The mapping the book previously listed as unrecoverable: which cluster became",
             "which cell type."),
       "Annotation"),
  LINK("celltype_marker_panel", "tables/celltype_marker_panel.csv",
       "Cell-type marker panel used", "The marker sets behind the annotation.", "Annotation"),
  LINK("celltype_confusion", "tables/celltype_confusion.csv",
       "Annotation confusion matrix", "Assigned vs reference-predicted labels.", "Annotation"),
  LINK("celltype_recall", "tables/celltype_recall.csv",
       "Annotation recall per cell type", "Per-type recall against the reference.", "Annotation"),

  # ---- trajectory --------------------------------------------------------
  LINK("pseudotime_validation", "tables/pseudotime_validation.csv",
       "Pseudotime validation", "Checks on the trajectory fit.", "Trajectory"),
  LINK("cm_stage_genes", "tables/cm_stage_genes.csv",
       "CM stage genes", "Genes marking maturation stage along the trajectory.", "Trajectory")
)

read_one <- function(rel) {
  p <- file.path(OUR, "results", rel)
  if (!file.exists(p)) return(NULL)
  d <- tryCatch(utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  # an unnamed first column is DESeq2's rowname dump; it is the gene
  if (names(d)[1] %in% c("", "X") && !"gene" %in% names(d)) names(d)[1] <- "gene"
  for (nm in names(d)) if (is.numeric(d[[nm]])) {
    d[[nm]] <- if (all(is.na(d[[nm]]) | abs(d[[nm]]) < 1e-4 | d[[nm]] == 0, na.rm = TRUE))
      signif(d[[nm]], 4) else round(d[[nm]], 4)
  }
  d
}

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
CONF <- app$confound

# ---- 1. link the computed tables -----------------------------------------
cat("\n== linking computed results ==\n")
linked <- list(); man <- list(); miss <- character(0)
for (w in WANT) {
  d <- read_one(w$file)
  if (is.null(d)) { miss <- c(miss, w$file); cat(sprintf("  %-30s MISSING (%s)\n", w$key, w$file)); next }
  if ("gene" %in% names(d)) d$confounder <- d$gene %in% CONF
  linked[[w$key]] <- d
  man[[length(man) + 1L]] <- data.frame(
    key = w$key, label = w$label, group = w$group, rows = nrow(d),
    cols = paste(names(d), collapse = ", "), note = w$note,
    source = w$file, stringsAsFactors = FALSE)
  cat(sprintf("  %-30s %6d rows  %s\n", w$key, nrow(d), w$file))
}
manifest <- do.call(rbind, man)

# ---- 2. un-trim the stored DE tables -------------------------------------
# build_app_data.R applied baseMean >= 3. Re-read the complete CSVs and keep every
# row, in the same schema the app already expects.
cat("\n== re-reading the pooled DE tables WITHOUT the baseMean filter ==\n")
KEEP <- c("gene","log2FoldChange","baseMean","pvalue","padj","confounder")
untrim <- function(path) {
  d <- read_one(path); if (is.null(d)) return(NULL)
  if (!"gene" %in% names(d)) return(NULL)
  d$confounder <- d$gene %in% CONF
  d <- d[, intersect(KEEP, names(d)), drop = FALSE]
  d <- d[is.finite(d$log2FoldChange), , drop = FALSE]
  d$log2FoldChange <- round(d$log2FoldChange, 3)
  if ("baseMean" %in% names(d)) d$baseMean <- round(d$baseMean, 1)
  if ("pvalue" %in% names(d)) d$pvalue <- signif(d$pvalue, 3)
  if ("padj" %in% names(d)) d$padj <- signif(d$padj, 3)
  d[order(-abs(d$log2FoldChange)), , drop = FALSE]
}
n_before <- sum(vapply(app$tables$ct_DE, nrow, 0L)); n_after <- n_before; swapped <- 0L
for (k in names(app$tables$ct_DE)) {
  parts <- strsplit(k, "_", fixed = TRUE)[[1]]
  tp <- parts[1]; ct <- paste(parts[-1], collapse = "_")
  d <- untrim(sprintf("tables/percelltype_%s_%s_KOvsWT.descriptive.DE.csv", tp, ct))
  if (!is.null(d) && nrow(d) >= nrow(app$tables$ct_DE[[k]])) {
    cat(sprintf("  %-24s %6d -> %6d rows\n", k, nrow(app$tables$ct_DE[[k]]), nrow(d)))
    app$tables$ct_DE[[k]] <- d; swapped <- swapped + 1L
  }
}
n_after <- sum(vapply(app$tables$ct_DE, nrow, 0L))
cat(sprintf("  %d of %d tables replaced; %s -> %s rows total (+%.1f%%)\n", swapped,
            length(app$tables$ct_DE), format(n_before, big.mark=","),
            format(n_after, big.mark=","), 100*(n_after/n_before - 1)))

app$linked <- linked
app$linked_manifest <- manifest
app$linked_built <- Sys.Date()

cat(sprintf("\n== summary ==\n  linked tables : %d\n  missing       : %d\n  linked size   : %.1f MB\n",
            length(linked), length(miss), as.numeric(object.size(linked))/1024^2))
cat(sprintf("  bundle size   : %.0f MB in memory\n", as.numeric(object.size(app))/1024^2))
if (length(miss)) cat("  MISSING:\n", paste0("    - ", miss, collapse = "\n"), "\n")

if (PROBE) { cat("\n--probe: nothing written.\n"); quit(save = "no") }

cat("\nBacking up -> app_data.pre_linked.bak.rds\n")
file.copy("app_data.rds", "app_data.pre_linked.bak.rds", overwrite = TRUE)
cat("Saving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
