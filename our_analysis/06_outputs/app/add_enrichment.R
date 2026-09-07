#!/usr/bin/env Rscript
# One-time augment: bundle the precomputed enrichment tables (GSEA / GO BP / TF
# activity) into app_data.rds as app$enrich, so the app can show them with no
# runtime network calls or heavy Bioconductor packages. Re-runnable.
#
#   Rscript add_enrichment.R [path/to/app_data.rds]

args <- commandArgs(trailingOnly = TRUE)
RDS  <- if (length(args)) args[1] else
  "C:/Users/Justi/OneDrive/Documents/GitHub/e2f-heart-scrna/shiny_app/app_data.rds"
TAB  <- "C:/Users/Justi/OneDrive - Marquette University/Personal/E2F 7_8/Han_scRNA_2025/our_analysis/results/tables"
stopifnot(file.exists(RDS), dir.exists(TAB))

rd <- function(f) {
  p <- file.path(TAB, f)
  if (!file.exists(p)) { cat("  MISSING:", f, "\n"); return(NULL) }
  d <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("  %s: %d rows x %d cols\n", f, nrow(d), ncol(d))); d
}

app <- readRDS(RDS)

cat("Reading enrichment CSVs ...\n")
gsea <- rd("pathway_fgsea_hallmark_kegg_e2f.csv")
go   <- rd("pathway_GO_BP_up_in_KO.csv")
tf   <- rd("tf_activity_by_celltype_genotype.csv")

# light trims for display size (keep leadingEdge / geneID strings intact)
if (!is.null(gsea)) {
  for (c in c("NES","ES","padj","pval","log2err")) if (c %in% names(gsea)) gsea[[c]] <- signif(gsea[[c]], 3)
}
if (!is.null(go)) {
  for (c in c("FoldEnrichment","RichFactor","zScore","pvalue","p.adjust","qvalue"))
    if (c %in% names(go)) go[[c]] <- signif(go[[c]], 3)
}
if (!is.null(tf)) tf$mean_activity <- signif(tf$mean_activity, 4)

app$enrich <- list(gsea = gsea, go = go, tf = tf)

cat("\nsummary:\n")
cat("  gsea celltypes:", paste(sort(unique(gsea$celltype)), collapse=", "), "\n")
cat("  gsea timepoints:", paste(sort(unique(gsea$timepoint)), collapse=", "), "\n")
cat("  go celltypes:", paste(sort(unique(go$celltype)), collapse=", "), "\n")
cat("  tf celltypes:", paste(sort(unique(tf$celltype)), collapse=", "), "\n")

bak <- sub("\\.rds$", ".pre_enrich.bak.rds", RDS)
if (!file.exists(bak)) file.copy(RDS, bak)
saveRDS(app, RDS, compress = "xz")
cat("saved:", RDS, "\n=== DONE add_enrichment ===\n")
