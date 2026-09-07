#!/usr/bin/env Rscript
# Phase B -- annotate ALL major cardiac cell types on the Harmony-integrated object
# (not just cardiomyocytes), then DESCRIBE KO-vs-WT composition (n=1: no testing).
#   primary  : SingleR cluster-level vs celldex MouseRNAseqData (broad lineages)
#   validate : canonical-marker module argmax per cluster (heart-specific)
# Reads processing/seurat.combined.rds; writes labeled object + dotplot +
# composition table/stacked bar.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(ggplot2); library(tidyr) }))
have_singler <- requireNamespace("SingleR", quietly = TRUE) &&
                requireNamespace("celldex", quietly = TRUE)

comb <- readRDS(file.path(PROC, "seurat.combined.rds"))
DefaultAssay(comb) <- "RNA"
comb <- NormalizeData(comb, verbose = FALSE)

## --- SingleR (cluster-level, broad lineages) -- optional cross-check ---------
# Final labels are marker-based (below); SingleR/celldex is only a validation column
# and is skipped if celldex's ExperimentHub dependency stack isn't installed.
cl2lab <- NULL
if (have_singler) {
  ref  <- tryCatch(celldex::MouseRNAseqData(), error = function(e) { message("celldex ref failed: ", conditionMessage(e)); NULL })
  if (!is.null(ref)) {
    mat  <- GetAssayData(comb, assay = "RNA", layer = "data")
    pred <- SingleR::SingleR(test = mat, ref = ref, labels = ref$label.main,
                             clusters = comb$seurat_clusters)
    cl2lab <- setNames(pred$labels, rownames(pred))
    comb$celltype_singler <- unname(cl2lab[as.character(comb$seurat_clusters)])
  }
}
if (is.null(cl2lab)) { comb$celltype_singler <- NA_character_; message("SingleR cross-check skipped (celldex unavailable).") }

## --- canonical-marker module argmax per cluster (validation) ----------------
ms.cols <- c()
for (ct in names(CELLTYPE_MARKERS)) {
  mk <- intersect(CELLTYPE_MARKERS[[ct]], rownames(comb))
  if (length(mk) < 2) next
  comb <- AddModuleScore(comb, features = list(mk), name = paste0("ms_", ct))
  ms.cols[ct] <- paste0("ms_", ct, "1")
}
cl.means <- sapply(ms.cols, function(col) tapply(comb[[col]][,1], comb$seurat_clusters, mean))
marker.lab <- colnames(cl.means)[max.col(cl.means, ties.method = "first")]
names(marker.lab) <- rownames(cl.means)
comb$celltype_marker <- unname(marker.lab[as.character(comb$seurat_clusters)])

## final label: marker-based heart-specific call (SingleR kept as a column too)
comb$celltype <- comb$celltype_marker
saveRDS(comb, file.path(PROC, "seurat.combined.annotated.rds"))

## --- validation dotplot ------------------------------------------------------
flat <- unique(unlist(lapply(CELLTYPE_MARKERS, function(g) intersect(g, rownames(comb)))))
png(file.path(OUTFIG, "annotation_marker_dotplot.png"), 1300, 700)
print(DotPlot(comb, features = flat, group.by = "celltype") +
        RotatedAxis() + ggtitle("Canonical markers by assigned cell type"))
dev.off()
png(file.path(OUTFIG, "annotation_umap.png"), 1100, 600)
print(DimPlot(comb, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("Cell types (marker-based)"))
dev.off()

## --- per-cluster label reconciliation table ---------------------------------
recon <- data.frame(cluster = rownames(cl.means),
                    singleR = if (is.null(cl2lab)) NA_character_ else cl2lab[rownames(cl.means)],
                    marker  = marker.lab[rownames(cl.means)])
write.csv(recon, file.path(OUTTAB, "annotation_cluster_labels.csv"), row.names = FALSE)

## --- composition (DESCRIPTIVE, n=1) -----------------------------------------
tab <- as.data.frame.matrix(table(comb$celltype, comb$orig.ident))
prop <- sweep(tab, 2, colSums(tab), "/")
comp <- data.frame(celltype = rownames(prop), round(prop, 4), check.names = FALSE)
for (tp in TIMEPOINTS) {
  wt <- paste0(tp, "WT"); ko <- paste0(tp, "KO")
  if (all(c(wt, ko) %in% colnames(prop)))
    comp[[paste0(tp, "_log2FC_KOvsWT")]] <-
      round(log2((prop[[ko]] + 1e-4) / (prop[[wt]] + 1e-4)), 3)
}
write.csv(comp, file.path(OUTTAB, "composition.csv"), row.names = FALSE)

longp <- prop; longp$celltype <- rownames(longp)
longp <- pivot_longer(longp, -celltype, names_to = "group", values_to = "prop")
png(file.path(OUTFIG, "composition_stacked.png"), 800, 600)
print(ggplot(longp, aes(group, prop, fill = celltype)) + geom_col() +
        ylab("proportion") + ggtitle("Cell-type composition (descriptive, n=1)") +
        theme_minimal())
dev.off()

cat("\n--- cluster labels (SingleR vs marker) ---\n"); print(recon, row.names = FALSE)
cat("\n--- composition (proportions + KO/WT log2FC; DESCRIPTIVE) ---\n"); print(comp, row.names = FALSE)
cat("\n=== DONE annotate ===\n")
