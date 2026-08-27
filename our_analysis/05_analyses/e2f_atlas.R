#!/usr/bin/env Rscript
# E2F atlas (parity with their 08_e2f_atlas.qmd): E2f7/E2f8 spatial expression in WT
# (FeaturePlot split by timepoint) + E2F-family & canonical-target DotPlot by cell
# type, split by genotype. Descriptive.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(ggplot2)))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(comb) <- "RNA"; comb <- NormalizeData(comb, verbose = FALSE)
comb$genotype <- genotype_of(comb$orig.ident)

## E2f7/E2f8 FeaturePlot in WT, split by timepoint ---------------------------
wt <- comb[, comb$genotype == "WT"]
e2f_present <- intersect(c("E2f7","E2f8"), rownames(wt))
if (length(e2f_present) && "umap" %in% names(wt@reductions)) {
  png(file.path(OUTFIG, "e2f_featureplot_WT.png"), 1100, 520 * length(e2f_present))
  print(FeaturePlot(wt, features = e2f_present, split.by = "timepoint", reduction = "umap", order = TRUE))
  dev.off()
  cat("  E2f7/E2f8 WT FeaturePlot written for:", paste(e2f_present, collapse=", "), "\n")
}

## E2F family + canonical targets DotPlot, by cell type, split by genotype ----
e2f_family <- intersect(paste0("E2f", 1:8), rownames(comb))
feats <- unique(c(e2f_family, intersect(E2F_TARGETS, rownames(comb))))
png(file.path(OUTFIG, "e2f_family_dotplot.png"), 1300, 700)
print(DotPlot(comb, features = feats, group.by = "celltype", split.by = "genotype",
              cols = c("blue","red")) + RotatedAxis() +
        ggtitle("E2F family + canonical targets by cell type (KO vs WT; descriptive)"))
dev.off()
cat("  E2F-family dotplot written;", length(e2f_family), "E2F genes,", length(feats), "features total\n")
cat("=== DONE e2f_atlas ===\n")
