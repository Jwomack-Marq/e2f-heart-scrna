#!/usr/bin/env Rscript
# Phase D (cross-timepoint) -- DESCRIPTIVE P0-vs-P7 maturation axis within each
# genotype, in cardiomyocytes. Lanes act as the (technical) replicates, so this is
# descriptive too -- but the P0->P7 maturation signal is large and biologically real.
# Reads seurat.{P0,P7}.merge.rds; writes results/tables/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(DESeq2)))

# Build CM pseudobulk (per orig.ident_lane) for each timepoint, then combine.
pb_list <- list()
for (tp in TIMEPOINTS) {
  obj <- readRDS(merged_path(tp))
  cm <- obj[, select_cardiac(obj)]          # uses default SCT assay (has data layer)
  cm$pb.sample <- paste(cm$orig.ident, cm$lane, sep = "_")
  m <- AggregateExpression(cm, assays = "RNA", group.by = "pb.sample", slot = "counts")$RNA
  pb_list[[tp]] <- as.matrix(m)
  rm(obj, cm); gc(verbose = FALSE)
}
common <- Reduce(intersect, lapply(pb_list, rownames))
pb <- do.call(cbind, lapply(pb_list, function(m) m[common, ]))
colnames(pb) <- gsub("-", "_", colnames(pb))   # AggregateExpression rewrites _ -> -

for (gt in c("WT", "KO")) {
  cols <- grep(gt, colnames(pb), value = TRUE)        # P0WT_lane1.. P7WT_lane6 (4 cols)
  tpf  <- ifelse(grepl("^P0", cols), "P0", "P7")
  coldata <- data.frame(row.names = cols, timepoint = relevel(factor(tpf), ref = "P0"))
  dds <- DESeqDataSetFromMatrix(round(pb[, cols]), colData = coldata, design = ~ timepoint)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)
  coef <- grep("timepoint_P7_vs_P0", resultsNames(dds), value = TRUE)
  res <- as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm"))
  res$gene <- rownames(res); res$NOTE <- "descriptive_P7vsP0_lanes_as_techreps"
  res <- res[order(-abs(res$log2FoldChange)), ]
  write.csv(res, file.path(OUTTAB, paste0("crosstime_", gt, "_P7vsP0.cardiac.descriptive.DE.csv")), row.names = FALSE)
  cat(sprintf("%s P7-vs-P0 CM: top |LFC| genes: %s\n", gt,
              paste(head(res$gene, 10), collapse = ", ")))
}
cat("=== DONE cross_timepoint ===\n")
