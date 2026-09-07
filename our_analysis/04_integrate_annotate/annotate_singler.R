#!/usr/bin/env Rscript
# Parity with their 03_annotate.qmd: SingleR vs celldex::MouseRNAseqData (label.fine
# + label.main), cluster-level, on the combined object. Compares to our marker-based
# labels. Degrades gracefully if celldex (ExperimentHub/alabaster stack) won't install.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))

ok <- requireNamespace("SingleR", quietly = TRUE) && requireNamespace("celldex", quietly = TRUE)
if (!ok) {
  msg <- "celldex/SingleR unavailable (ExperimentHub/alabaster dependency stack failed to install). SingleR annotation skipped; our marker-based annotation (annotate.R) is the fallback and is heart-appropriate (celldex MouseRNAseqData is a coarse bulk reference)."
  writeLines(msg, file.path(OUTTAB, "singler_blocked.txt")); cat(msg, "\n=== DONE annotate_singler ===\n"); quit(save = "no")
}
suppressWarnings(suppressMessages({ library(SingleR); library(celldex) }))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(comb) <- "RNA"; comb <- NormalizeData(comb, verbose = FALSE)
ref <- tryCatch(celldex::MouseRNAseqData(), error = function(e) NULL)
if (is.null(ref)) { writeLines(paste("celldex reference download failed."), file.path(OUTTAB, "singler_blocked.txt"))
  cat("celldex ref failed.\n=== DONE annotate_singler ===\n"); quit(save = "no") }

mat <- GetAssayData(comb, assay = "RNA", layer = "data")
for (lab in c("label.fine", "label.main")) {
  pred <- SingleR::SingleR(test = mat, ref = ref, labels = ref[[lab]], clusters = comb$seurat_clusters)
  comb[[paste0("singler_", sub("label.", "", lab, fixed = TRUE))]] <-
    unname(setNames(pred$labels, rownames(pred))[as.character(comb$seurat_clusters)])
}
# agreement vs our marker-based label
agree <- as.data.frame.matrix(table(comb$celltype, comb$singler_main))
write.csv(agree, file.path(OUTTAB, "singler_vs_marker_confusion.csv"))
cluster_tab <- unique(data.frame(cluster = comb$seurat_clusters,
                                 marker = comb$celltype,
                                 singler_main = comb$singler_main,
                                 singler_fine = comb$singler_fine))
cluster_tab <- cluster_tab[order(as.numeric(as.character(cluster_tab$cluster))), ]
write.csv(cluster_tab, file.path(OUTTAB, "singler_cluster_labels.csv"), row.names = FALSE)
saveRDS(comb, file.path(PROC, "seurat.combined.annotated.rds"))   # add singler columns

cat("\n--- cluster labels: marker vs SingleR ---\n"); print(cluster_tab, row.names = FALSE)
cat("\n--- confusion (marker rows x SingleR.main cols) ---\n"); print(agree)
cat("=== DONE annotate_singler ===\n")
