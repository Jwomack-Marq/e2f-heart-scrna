#!/usr/bin/env Rscript
# Per-cell-type KO-vs-WT DESCRIPTIVE DE (n=1; sex-confounded) -- extends the
# cardiomyocyte-only contrasts to every annotated lineage, per timepoint.
# Pseudobulk over the 2 lanes (technical) -> DESeq2 ~genotype -> apeglm shrinkage.
# Reads processing/seurat.combined.annotated.rds; writes results/tables/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(DESeq2)))

MIN_CELLS_PER_PB <- 20    # min cells per (genotype,lane) pseudobulk sample
# Known confounder genes (sex + the ROSA26 knock-in construct) -- flagged, not dropped.
CONFOUND <- c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(comb) <- "RNA"
comb$pb.sample <- paste(comb$orig.ident, comb$lane, sep = "_")

summary_rows <- list()
for (tp in TIMEPOINTS) {
  for (ct in sort(unique(comb$celltype))) {
    cells <- which(comb$timepoint == tp & comb$celltype == ct)
    if (length(cells) < 4 * MIN_CELLS_PER_PB) {
      summary_rows[[paste(tp,ct)]] <- data.frame(timepoint=tp, celltype=ct,
        n_cells=length(cells), status="skipped_too_few", n_DE_absLFC_gt1=NA,
        top_KO_up=NA, top_KO_down=NA); next
    }
    sub <- comb[, cells]
    counts_per <- table(sub$pb.sample)
    if (length(counts_per) < 4 || min(counts_per) < MIN_CELLS_PER_PB) {
      summary_rows[[paste(tp,ct)]] <- data.frame(timepoint=tp, celltype=ct,
        n_cells=length(cells), status="skipped_unbalanced", n_DE_absLFC_gt1=NA,
        top_KO_up=NA, top_KO_down=NA); next
    }
    pb <- AggregateExpression(sub, assays = "RNA", group.by = "pb.sample", slot = "counts")$RNA
    cond <- ifelse(grepl(paste0(tp,"KO"), colnames(pb)), paste0(tp,"KO"), paste0(tp,"WT"))
    cd <- data.frame(row.names = colnames(pb),
                     condition = relevel(factor(cond), ref = paste0(tp,"WT")))
    dds <- DESeqDataSetFromMatrix(round(as.matrix(pb)), colData = cd, design = ~ condition)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) NULL)
    if (is.null(dds)) {
      summary_rows[[paste(tp,ct)]] <- data.frame(timepoint=tp, celltype=ct,
        n_cells=length(cells), status="DESeq_failed", n_DE_absLFC_gt1=NA,
        top_KO_up=NA, top_KO_down=NA); next
    }
    coef <- grep(paste0("condition_", tp, "KO_vs_", tp, "WT"), resultsNames(dds), value = TRUE)
    res <- as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm"))
    res$gene <- rownames(res)
    res$confounder <- res$gene %in% CONFOUND
    res$NOTE <- "descriptive_n1_sex_confounded"
    res <- res[order(-abs(res$log2FoldChange)), ]
    ctf <- gsub("[^A-Za-z0-9]+", "_", ct)
    write.csv(res, file.path(OUTTAB, sprintf("percelltype_%s_%s_KOvsWT.descriptive.DE.csv", tp, ctf)), row.names = FALSE)

    bio <- res[!res$confounder, ]                       # exclude sex/construct genes
    up  <- head(bio$gene[bio$log2FoldChange >  1], 6)
    dn  <- head(bio$gene[bio$log2FoldChange < -1], 6)
    summary_rows[[paste(tp,ct)]] <- data.frame(timepoint=tp, celltype=ct,
      n_cells=length(cells), status="ok",
      n_DE_absLFC_gt1 = sum(abs(bio$log2FoldChange) > 1, na.rm=TRUE),
      top_KO_up = paste(up, collapse=","), top_KO_down = paste(dn, collapse=","))
    cat(sprintf("  %s %-16s n=%5d  KO-up: %s\n", tp, ct, length(cells), paste(up, collapse=", ")))
  }
}
summ <- do.call(rbind, summary_rows)
write.csv(summ, file.path(OUTTAB, "percelltype_KOvsWT_summary.csv"), row.names = FALSE)
cat("\n--- per-cell-type KO-vs-WT summary (biological genes; sex/construct excluded) ---\n")
print(summ[, c("timepoint","celltype","n_cells","status","n_DE_absLFC_gt1","top_KO_up")], row.names = FALSE)
cat("\n=== DONE per_celltype_de ===\n")
