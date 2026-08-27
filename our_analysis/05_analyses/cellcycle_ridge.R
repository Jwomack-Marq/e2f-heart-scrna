#!/usr/bin/env Rscript
# Reproduce the prior analysis's cell-cycle "sorting" plot: RidgePlot of S.Score and
# G2M.Score grouped by cluster (their scRNA_merge{P0,P7}.html), per timepoint.
# Plus companion Phase UMAP and a KO-vs-WT phase split. Descriptive.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(ggplot2)))

for (tp in TIMEPOINTS) {
  obj <- readRDS(merged_path(tp))
  obj$genotype <- genotype_of(obj$orig.ident)
  obj <- SetIdent(obj, value = "seurat_clusters")

  png(file.path(OUTFIG, paste0("cellcycle_ridge_", tp, ".png")), 1100, 700)
  print(RidgePlot(obj, features = c("S.Score", "G2M.Score"), ncol = 2) +
          patchwork::plot_annotation(title = paste0(tp, ": cell-cycle scores by cluster (RidgePlot)")))
  dev.off()

  if ("umap" %in% names(obj@reductions)) {
    png(file.path(OUTFIG, paste0("cellcycle_phase_umap_", tp, ".png")), 1200, 520)
    print((DimPlot(obj, group.by = "Phase", reduction = "umap") + ggtitle(paste0(tp, " phase"))) |
          (DimPlot(obj, group.by = "Phase", reduction = "umap", split.by = "genotype") + ggtitle("by genotype")))
    dev.off()
  }
  cat(sprintf("  %s: ridge + phase UMAP written; phase x genotype:\n", tp))
  print(table(obj$Phase, obj$genotype))
  rm(obj); gc(verbose = FALSE)
}
cat("=== DONE cellcycle_ridge ===\n")
