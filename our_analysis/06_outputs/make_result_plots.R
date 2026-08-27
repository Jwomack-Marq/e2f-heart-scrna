#!/usr/bin/env Rscript
# Result figures for the descriptive analyses (DE/GSEA/doublets/composition were
# table-only). All KO-vs-WT panels are labeled DESCRIPTIVE (n=1, sex-confounded);
# p-values are shown for ranking only, not inference.
# Reads results/tables/ (+ the annotated object for the candidate dotplot);
# writes results/figures/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(ggplot2); library(tidyr) }))
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
CONFOUND <- c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")
savep <- function(f, p, w=1000, h=750) { png(file.path(OUTFIG, f), w, h); print(p); dev.off() }

## 1) Volcano + MA per timepoint (cardiac KO-vs-WT, descriptive) --------------
for (tp in TIMEPOINTS) tryCatch({
  d <- read.csv(file.path(OUTTAB, paste0(tp, ".cardiac.descriptive.DE.csv")))
  d <- d[is.finite(d$log2FoldChange), ]
  d$class <- "ns"
  d$class[abs(d$log2FoldChange) > 1] <- "|LFC|>1 (biological)"
  d$class[d$gene %in% CONFOUND]      <- "sex/construct (confounder)"
  d$neglog <- -log10(pmax(d$pvalue, 1e-300))   # add BEFORE building `lab` (the label layer inherits it)
  cols <- c("ns"="grey80","|LFC|>1 (biological)"="firebrick","sex/construct (confounder)"="grey35")
  bio <- d[d$class == "|LFC|>1 (biological)" & is.finite(d$neglog), ]
  lab <- head(bio[order(-abs(bio$log2FoldChange)), ], 14)
  v <- ggplot(d, aes(log2FoldChange, neglog, color = class)) +
    geom_point(size = 0.6, alpha = 0.6) + scale_color_manual(values = cols) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    labs(title = paste0(tp, " cardiomyocytes — KO vs WT (DESCRIPTIVE, n=1)"),
         subtitle = "apeglm-shrunken log2FC; p-value axis for RANKING only — NOT valid (n=1, sex-confounded)",
         x = "shrunken log2FC (KO/WT)", y = "-log10(p) [not valid]") + theme_bw()
  if (has_repel && nrow(lab)) v <- v + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = 20, color = "black")

  ma <- ggplot(d, aes(log10(baseMean + 1), log2FoldChange, color = class)) +
    geom_point(size = 0.6, alpha = 0.6) + scale_color_manual(values = cols) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    labs(title = paste0(tp, " — MA (effect size vs expression)"),
         x = "log10(mean expression)", y = "shrunken log2FC (KO/WT)") + theme_bw()
  if (has_repel && nrow(lab)) ma <- ma + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = 20, color = "black")

  savep(paste0("DE_", tp, "_cardiac_volcano_MA.png"), v | ma, 1500, 700)
  cat("  wrote volcano/MA for", tp, "\n")
}, error = function(e) message("volcano ", tp, " failed: ", conditionMessage(e)))

## 2) GSEA dotplot (P0) -------------------------------------------------------
tryCatch({
  g <- read.csv(file.path(OUTTAB, "P0.cardiac.GSEA_BP.csv"))
  g <- head(g[order(-abs(g$NES)), ], 15)
  p <- ggplot(g, aes(NES, reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point() + scale_color_gradient(low = "firebrick", high = "blue") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(title = "P0 cardiomyocyte GSEA (GO:BP), KO-vs-WT ranked", y = NULL) + theme_bw()
  savep("GSEA_P0_cardiac_BP.png", p, 1100, 650)
  cat("  wrote GSEA dotplot\n")
}, error = function(e) message("GSEA plot failed: ", conditionMessage(e)))

## 3) Doublet comparison (scDblFinder vs Scrublet) ---------------------------
tryCatch({
  dd <- read.csv(file.path(OUTTAB, "doublet_comparison.csv"))
  long <- pivot_longer(dd, c(scrublet_pct, scDblFinder_pct), names_to = "method", values_to = "pct")
  p <- ggplot(long, aes(lane, pct, fill = method)) +
    geom_col(position = "dodge") + coord_flip() +
    labs(title = "Doublet rate: Scrublet under-called vs scDblFinder", y = "% doublets", x = NULL) + theme_bw()
  savep("QC_doublet_comparison.png", p, 900, 600)
  cat("  wrote doublet barplot\n")
}, error = function(e) message("doublet plot failed: ", conditionMessage(e)))

## 4) Per-cell-type DE-gene counts -------------------------------------------
tryCatch({
  s <- read.csv(file.path(OUTTAB, "percelltype_KOvsWT_summary.csv"))
  s <- s[s$status == "ok" & !is.na(s$n_DE_absLFC_gt1), ]
  p <- ggplot(s, aes(reorder(celltype, n_DE_absLFC_gt1), n_DE_absLFC_gt1, fill = timepoint)) +
    geom_col(position = "dodge") + coord_flip() +
    labs(title = "KO-vs-WT genes |shrunken LFC|>1 by cell type (DESCRIPTIVE, n=1)",
         subtitle = "sex/construct genes excluded", x = NULL, y = "# genes |LFC|>1") + theme_bw()
  savep("DE_percelltype_counts.png", p, 950, 600)
  cat("  wrote per-cell-type DE-count barplot\n")
}, error = function(e) message("percelltype plot failed: ", conditionMessage(e)))

## 5) Candidate genes across cell types, split by genotype -------------------
tryCatch({
  comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
  DefaultAssay(comb) <- "RNA"; comb <- NormalizeData(comb, verbose = FALSE)
  comb$genotype <- genotype_of(comb$orig.ident)
  cand <- intersect(c("Gabbr2","Tcf4","Adamts9","Ralyl"), rownames(comb))
  p <- DotPlot(comb, features = cand, group.by = "celltype", split.by = "genotype",
               cols = c("blue","red")) + RotatedAxis() +
    labs(title = "Recurring KO-up candidates by cell type (KO vs WT; descriptive)")
  savep("candidates_Gabbr2_dotplot.png", p, 1000, 700)
  # focused Gabbr2 violin in cardiomyocytes
  cm <- comb[, comb$celltype == "Cardiomyocyte"]
  savep("candidates_Gabbr2_CM_violin.png",
        VlnPlot(cm, "Gabbr2", group.by = "genotype", split.by = "timepoint", pt.size = 0) +
          ggtitle("Gabbr2 in cardiomyocytes (descriptive)"), 800, 500)
  cat("  wrote candidate dotplot + Gabbr2 violin\n")
}, error = function(e) message("candidate plot failed: ", conditionMessage(e)))

cat("=== DONE make_result_plots ===\n")
