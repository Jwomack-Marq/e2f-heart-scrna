#!/usr/bin/env Rscript
# Phase D -- DESCRIPTIVE KO-vs-WT DE in cardiomyocytes (n=1: NO valid p-values).
# Ranks genes by apeglm-SHRUNKEN log2FC (effect size), runs rank-based GSEA (more
# robust than threshold ORA at n=1), and intersects P0 & P7 for consistent signal.
# Reads seurat.{P0,P7}.merge.rds; writes results/tables/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({
  library(DESeq2); library(clusterProfiler); library(org.Mm.eg.db)
}))

up_lists <- list()
for (tp in TIMEPOINTS) {
  cat("\n==============", tp, "==============\n")
  obj <- readRDS(merged_path(tp))
  cm <- obj[, select_cardiac(obj)]          # uses default SCT assay (has data layer)
  cm$pb.sample <- paste(cm$orig.ident, cm$lane, sep = "_")
  pb <- AggregateExpression(cm, assays = "RNA", group.by = "pb.sample", slot = "counts")$RNA

  cond <- ifelse(grepl(paste0(tp, "KO"), colnames(pb)), paste0(tp, "KO"), paste0(tp, "WT"))
  coldata <- data.frame(row.names = colnames(pb), condition = relevel(factor(cond), ref = paste0(tp, "WT")))
  dds <- DESeqDataSetFromMatrix(round(as.matrix(pb)), colData = coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)
  coef <- grep(paste0("condition_", tp, "KO_vs_", tp, "WT"), resultsNames(dds), value = TRUE)
  res  <- as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm"))   # shrunken effect size
  res$gene <- rownames(res)
  res <- res[order(-abs(res$log2FoldChange)), ]
  res$NOTE <- "descriptive_n1_no_valid_pvalues"
  write.csv(res, file.path(OUTTAB, paste0(tp, ".cardiac.descriptive.DE.csv")), row.names = FALSE)
  cat(sprintf("  %d genes; top KO-up by shrunken LFC: %s\n", nrow(res),
              paste(head(res$gene[res$log2FoldChange > 0], 8), collapse = ", ")))
  up_lists[[tp]] <- res$gene[res$log2FoldChange > 1]

  ## rank-based GSEA on shrunken LFC (guarded) --------------------------------
  gl <- sort(setNames(res$log2FoldChange, res$gene), decreasing = TRUE)
  gl <- gl[is.finite(gl)]
  gsea <- tryCatch(
    clusterProfiler::gseGO(gl, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                           minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.25, verbose = FALSE),
    error = function(e) { message("GSEA failed for ", tp, ": ", conditionMessage(e)); NULL })
  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    write.csv(as.data.frame(gsea), file.path(OUTTAB, paste0(tp, ".cardiac.GSEA_BP.csv")), row.names = FALSE)
    cat(sprintf("  GSEA: %d enriched BP terms\n", nrow(as.data.frame(gsea))))
  } else cat("  GSEA: no terms\n")
  rm(obj, cm); gc(verbose = FALSE)
}

## shared cross-timepoint KO-up signal ----------------------------------------
if (length(up_lists) == 2) {
  shared <- intersect(up_lists[[1]], up_lists[[2]])
  writeLines(shared, file.path(OUTTAB, "shared_KO_up_P0_and_P7.txt"))
  cat(sprintf("\nShared KO-up (shrunkenLFC>1 at BOTH P0 & P7): %d genes\n  %s\n",
              length(shared), paste(head(shared, 25), collapse = ", ")))
}
cat("\n=== DONE de_descriptive ===\n")
