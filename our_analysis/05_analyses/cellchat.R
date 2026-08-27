#!/usr/bin/env Rscript
# Cell-cell communication (DESCRIPTIVE, n=1, sex-confounded).
# CellChat per timepoint for KO and WT separately, then differential (KO vs WT)
# interaction counts/strength and pathway-level changes. No replicate stats at n=1.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
if (!requireNamespace("CellChat", quietly = TRUE)) {
  writeLines("CellChat not installed; communication analysis skipped.", file.path(OUTTAB, "cellchat_blocked.txt"))
  cat("CellChat missing.\n=== DONE cellchat ===\n"); quit(save = "no")
}
suppressWarnings(suppressMessages({ library(CellChat); library(ggplot2) }))
future::plan("sequential")

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds")); DefaultAssay(comb) <- "RNA"
comb <- NormalizeData(comb, verbose = FALSE)
comb$genotype <- genotype_of(comb$orig.ident)

build_cc <- function(obj) {
  data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
  meta <- data.frame(labels = as.character(obj$celltype), row.names = colnames(obj))
  cc <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")
  cc@DB <- CellChat::CellChatDB.mouse
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = "triMean")
  cc <- CellChat::filterCommunication(cc, min.cells = 10)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  cc
}

for (tp in TIMEPOINTS) tryCatch({
  cat("\n==============", tp, "==============\n")
  ccs <- list()
  for (g in c("WT","KO")) {
    sub <- comb[, comb$timepoint == tp & comb$genotype == g]
    sub$celltype <- droplevels(factor(sub$celltype))
    ccs[[g]] <- build_cc(sub)
    cat(sprintf("  %s %s: built CellChat (%d cells, %d cell types)\n", tp, g, ncol(sub), length(unique(sub$celltype))))
  }
  saveRDS(ccs, file.path(PROC, paste0("cellchat_", tp, "_objs.rds")))   # cache the slow computeCommunProb result
  # differential interaction matrices (align cell types)
  cts <- sort(union(rownames(ccs$WT@net$count), rownames(ccs$KO@net$count)))
  pad <- function(m) { full <- matrix(0, length(cts), length(cts), dimnames = list(cts, cts))
    full[rownames(m), colnames(m)] <- m; full }
  dcount <- pad(ccs$KO@net$count)  - pad(ccs$WT@net$count)
  dweight<- pad(ccs$KO@net$weight) - pad(ccs$WT@net$weight)
  write.csv(dcount,  file.path(OUTTAB, paste0("cellchat_", tp, "_diff_count.csv")))
  write.csv(dweight, file.path(OUTTAB, paste0("cellchat_", tp, "_diff_weight.csv")))

  # pathway-level total strength KO vs WT
  pw <- function(cc) { p <- cc@netP$prob; if (is.null(p)) return(setNames(numeric(0), character(0)))
    apply(p, 3, sum) }
  wt <- pw(ccs$WT); ko <- pw(ccs$KO); paths <- union(names(wt), names(ko))
  pwt <- data.frame(pathway = paths, WT = unname(round(wt[paths], 4)), KO = unname(round(ko[paths], 4)))
  pwt[is.na(pwt)] <- 0; pwt$diff_KO_WT <- round(pwt$KO - pwt$WT, 4)
  pwt <- pwt[order(-abs(pwt$diff_KO_WT)), ]
  write.csv(pwt, file.path(OUTTAB, paste0("cellchat_", tp, "_pathway_diff.csv")), row.names = FALSE)

  # differential interaction-count heatmap
  dd <- as.data.frame(as.table(dcount)); names(dd) <- c("source","target","diff")
  png(file.path(OUTFIG, paste0("cellchat_", tp, "_diff_heatmap.png")), 800, 700)
  print(ggplot(dd, aes(target, source, fill = diff)) + geom_tile() +
          scale_fill_gradient2(low="steelblue", mid="white", high="firebrick") +
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = paste0(tp, " differential interactions (KO - WT, count; descriptive)")))
  dev.off()
  cat(sprintf("  %s top pathways by |KO-WT|: %s\n", tp, paste(head(pwt$pathway, 8), collapse = ", ")))
}, error = function(e) message("CellChat ", tp, " failed: ", conditionMessage(e)))

cat("\n=== DONE cellchat ===\n")
