#!/usr/bin/env Rscript
# Parity with their 06_pathway.qmd: fgsea over Hallmark + KEGG(LEGACY) + E2F-target
# (TFT:GTRD) MSigDB sets, plus enrichGO BP ORA on KO-up genes. Run per (celltype,
# timepoint) using OUR descriptive DE tables, ranked by apeglm-shrunken log2FC.
# Sex/ROSA26 CONFOUNDERS are removed from the ranking (primary) since they would
# otherwise distort enrichment (their pipeline predates the sex finding).
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(fgsea); library(msigdbr); library(ggplot2) }))
CONFOUND <- c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")

## --- MSigDB gene sets (mouse) -- match their collections ---------------------
get_sets <- function() {
  m_h    <- msigdbr(species = "Mus musculus", collection = "H")
  m_kegg <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  m_tft  <- msigdbr(species = "Mus musculus", collection = "C3", subcollection = "TFT:GTRD")
  e2f    <- m_tft[grepl("E2F", m_tft$gs_name), ]
  ps <- c(split(m_h$gene_symbol, m_h$gs_name),
          split(m_kegg$gene_symbol, m_kegg$gs_name),
          split(e2f$gene_symbol, e2f$gs_name))
  lapply(ps, unique)
}
pathways <- tryCatch(get_sets(), error = function(e) { message("msigdbr fetch failed: ", conditionMessage(e)); NULL })
if (is.null(pathways)) { cat("Cannot build gene sets; aborting.\n=== DONE pathway_msigdb ===\n"); quit(save = "no") }
cat(sprintf("Loaded %d gene sets (Hallmark + KEGG_LEGACY + E2F TFT)\n", length(pathways)))

## --- collect our descriptive DE tables (celltype x timepoint) ---------------
de_files <- c(file.path(OUTTAB, paste0(TIMEPOINTS, ".cardiac.descriptive.DE.csv")),
              list.files(OUTTAB, pattern = "^percelltype_.*KOvsWT\\.descriptive\\.DE\\.csv$", full.names = TRUE))
de_files <- de_files[file.exists(de_files)]

parse_meta <- function(f) {
  b <- basename(f)
  if (grepl("^percelltype_", b)) {
    tp <- sub("^percelltype_(P[07])_.*$", "\\1", b)
    ct <- sub("^percelltype_P[07]_(.*)_KOvsWT.*$", "\\1", b)
  } else { tp <- substr(b, 1, 2); ct <- "Cardiomyocyte(cardiac-subset)" }
  list(tp = tp, ct = ct)
}

gsea_all <- list(); go_all <- list()
suppressWarnings(suppressMessages({ library(clusterProfiler); library(org.Mm.eg.db) }))
for (f in de_files) {
  m <- parse_meta(f); d <- read.csv(f)
  d <- d[!is.na(d$log2FoldChange) & !(d$gene %in% CONFOUND), ]
  ranks <- sort(setNames(d$log2FoldChange, d$gene), decreasing = TRUE)
  ranks <- ranks[!duplicated(names(ranks)) & is.finite(ranks)]
  if (length(ranks) >= 100) {
    r <- tryCatch(fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500),
                  error = function(e) NULL)
    if (!is.null(r) && nrow(r)) {
      r$leadingEdge <- sapply(r$leadingEdge, paste, collapse = ",")
      r$celltype <- m$ct; r$timepoint <- m$tp
      gsea_all[[f]] <- r
    }
  }
  up <- d$gene[d$padj < 0.05 & d$log2FoldChange > 1]; up <- up[!is.na(up)]
  if (length(up) >= 10) {
    ego <- tryCatch(enrichGO(up, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                             pAdjustMethod = "BH", qvalueCutoff = 0.1), error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego))) {
      g <- as.data.frame(ego); g$celltype <- m$ct; g$timepoint <- m$tp; go_all[[f]] <- g
    }
  }
  cat(sprintf("  %s %s: %d ranked genes\n", m$tp, m$ct, length(ranks)))
}

gsea <- do.call(rbind, gsea_all)
if (!is.null(gsea)) {
  gsea <- gsea[order(gsea$padj), ]
  write.csv(gsea, file.path(OUTTAB, "pathway_fgsea_hallmark_kegg_e2f.csv"), row.names = FALSE)
  write.csv(gsea[grepl("E2F", gsea$pathway, ignore.case = TRUE), ],
            file.path(OUTTAB, "pathway_fgsea_E2F_focus.csv"), row.names = FALSE)
  sig <- gsea[gsea$padj < 0.1, ]
  cat(sprintf("\nfgsea: %d total terms, %d at padj<0.1\n", nrow(gsea), nrow(sig)))
  if (nrow(sig)) {
    top <- head(sig[order(sig$padj), ], 25)
    top$grp <- paste(top$celltype, top$timepoint, sep = "@")
    p <- ggplot(top, aes(grp, reorder(pathway, NES), fill = NES, size = -log10(padj))) +
      geom_point(shape = 21) + scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "fgsea (Hallmark/KEGG/E2F) — top terms (descriptive, n=1)", x = NULL, y = NULL)
    ggsave(file.path(OUTFIG, "pathway_fgsea_dotplot.png"), p, width = 12, height = 9, dpi = 120)
  }
}
go <- do.call(rbind, go_all)
if (!is.null(go)) write.csv(go, file.path(OUTTAB, "pathway_GO_BP_up_in_KO.csv"), row.names = FALSE)
cat("=== DONE pathway_msigdb ===\n")
