#!/usr/bin/env Rscript
# Phase C -- E2F7/8 hypothesis readouts (DESCRIPTIVE, n=1 per condition).
#   1) KO verification: E2f7/E2f8 expression KO vs WT (sanity check the knockout)
#   2) E2F/cell-cycle target DE-REPRESSION in cardiomyocytes (KO should be higher)
#   3) cycling-CM fraction (S/G2M) KO vs WT  -- the cell-cycle-exit readout
#   4) ploidy/maturation SURROGATES (indirect; scRNA cannot measure DNA content)
# Reads processing/seurat.{P0,P7}.merge.rds; writes results/{tables,figures}/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(ggplot2)))

mean_pct <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj[["RNA"]]))
  dat <- FetchData(obj, vars = genes, layer = "data")          # normalized RNA
  cnt <- FetchData(obj, vars = genes, layer = "counts")
  data.frame(gene = genes,
             mean_logexpr = round(colMeans(dat), 4),
             pct_expressing = round(100 * colMeans(cnt > 0), 1),
             row.names = NULL)
}

ko_rows <- list(); cm_rows <- list()
for (tp in TIMEPOINTS) {
  cat("\n==============", tp, "==============\n")
  obj <- readRDS(merged_path(tp))
  obj$genotype <- genotype_of(obj$orig.ident)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)

  ## 1) KO verification (E2f7/E2f8), global ----------------------------------
  for (g in c("WT","KO")) {
    sub <- obj[, obj$genotype == g]
    mp <- mean_pct(sub, c("E2f7","E2f8"))
    mp$timepoint <- tp; mp$genotype <- g; mp$n_cells <- ncol(sub)
    ko_rows[[paste(tp,g)]] <- mp
    cat(sprintf("  %s %s  E2f7=%.3f(%.0f%%)  E2f8=%.3f(%.0f%%)\n", tp, g,
        mp$mean_logexpr[mp$gene=="E2f7"], mp$pct_expressing[mp$gene=="E2f7"],
        mp$mean_logexpr[mp$gene=="E2f8"], mp$pct_expressing[mp$gene=="E2f8"]))
  }
  png(file.path(OUTFIG, paste0("e2f_", tp, "_KO_verification.png")), 900, 450)
  print(VlnPlot(obj, features = c("E2f7","E2f8"), group.by = "genotype",
                pt.size = 0) + patchwork::plot_annotation(title = paste(tp, "KO verification")))
  dev.off()

  ## 2-4) cardiomyocyte-restricted readouts ----------------------------------
  is.cm <- select_cardiac(obj)
  cm <- obj[, is.cm]
  cat(sprintf("  cardiomyocytes: %d (%.1f%%)\n", ncol(cm), 100*ncol(cm)/ncol(obj)))

  e2f.tg <- intersect(E2F_TARGETS, rownames(cm))
  cm <- AddModuleScore(cm, features = list(e2f.tg),  name = "E2Ftarget")
  cm <- AddModuleScore(cm, features = list(intersect(CM_MATURE, rownames(cm))),   name = "Mature")
  cm <- AddModuleScore(cm, features = list(intersect(CM_IMMATURE, rownames(cm))), name = "Immature")

  for (g in c("WT","KO")) {
    m <- cm$genotype == g
    cm_rows[[paste(tp,g)]] <- data.frame(
      timepoint = tp, genotype = g, n_CM = sum(m),
      cycling_fraction = round(mean(cm$Phase[m] %in% c("S","G2M")), 3),
      E2F_target_score = round(median(cm$E2Ftarget1[m]), 4),
      mature_score     = round(median(cm$Mature1[m]), 4),
      immature_score   = round(median(cm$Immature1[m]), 4),
      median_nCount    = round(median(cm$nCount_RNA[m])),
      median_nFeature  = round(median(cm$nFeature_RNA[m])))
  }
  png(file.path(OUTFIG, paste0("e2f_", tp, "_CM_cellcycle_targets.png")), 1000, 450)
  print((VlnPlot(cm, "E2Ftarget1", group.by = "genotype", pt.size = 0) + ggtitle("E2F-target score (CM)")) |
        (VlnPlot(cm, "S.Score",   group.by = "genotype", pt.size = 0) + ggtitle("S score (CM)")))
  dev.off()
  rm(obj, cm); gc(verbose = FALSE)
}

ko <- do.call(rbind, ko_rows); cmt <- do.call(rbind, cm_rows)
write.csv(ko,  file.path(OUTTAB, "e2f_ko_verification.csv"), row.names = FALSE)
write.csv(cmt, file.path(OUTTAB, "e2f_cm_readouts.csv"),     row.names = FALSE)
cat("\n--- KO verification ---\n"); print(ko, row.names = FALSE)
cat("\n--- CM readouts (descriptive, n=1) ---\n"); print(cmt, row.names = FALSE)
cat("\n=== DONE e2f_readouts ===\n")
