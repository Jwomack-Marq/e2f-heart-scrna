#!/usr/bin/env Rscript
# Compare the rerun outputs against the original run.
#  - cell counts (impact of the new QC: doublet + mito + upper-cap filtering)
#  - clustering granularity (# clusters) and marker-gene agreement per condition
#  - the new pseudobulk DESeq2 results (rerun only)
suppressWarnings(suppressMessages(library(Seurat)))

args0 <- commandArgs(FALSE); sp <- sub("^--file=", "", args0[grep("^--file=", args0)])
# This script lives in our_analysis/06_outputs/. our_analysis is one level up;
# the original run is the sibling folder original_Han_analysis/.
our  <- normalizePath(file.path(dirname(sp), ".."))
orig <- normalizePath(file.path(our, "..", "original_Han_analysis"))
rer  <- our                                    # our (re-run) outputs live here
conds <- c("P0WT","P0KO","P7WT","P7KO")

ncol_rds <- function(f) if (file.exists(f)) ncol(readRDS(f)) else NA_integer_

cat("==================================================================\n")
cat("1) CELLS PER CONDITION (merged, post-QC)   original -> rerun\n")
cat("==================================================================\n")
cell_rows <- list()
for (cc in conds) {
  o <- ncol_rds(file.path(orig, "processing",          paste0("merge.lanes.", cc, ".rds")))
  r <- ncol_rds(file.path(rer,  "processing",          paste0("merge.lanes.", cc, ".rds")))
  pct <- if (!is.na(o) && o > 0) round(100*(o-r)/o, 1) else NA
  cell_rows[[cc]] <- data.frame(condition=cc, cells_original=o, cells_rerun=r,
                                removed=o-r, pct_removed=pct)
  cat(sprintf("  %-5s  %6s -> %6s   (-%s cells, -%s%%)\n", cc,
              o, r, o-r, pct))
}
cells <- do.call(rbind, cell_rows)

read_markers <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path)
  names(df)[1] <- "row"
  df
}
sig_genes <- function(df, padj = 0.05) {
  if (is.null(df) || !"p_val_adj" %in% names(df)) return(character(0))
  unique(df$gene[df$p_val_adj < padj])
}
jacc <- function(a,b) { u <- length(union(a,b)); if (u==0) NA else round(length(intersect(a,b))/u,3) }

cat("\n==================================================================\n")
cat("2) CLUSTERS & MARKER AGREEMENT per condition\n")
cat("   (#clusters; #sig marker genes padj<0.05; gene-set Jaccard; top30 overlap)\n")
cat("==================================================================\n")
mk_rows <- list()
for (cc in conds) {
  od <- read_markers(file.path(orig, "results","markers", paste0(cc, ".allmarkers.csv")))
  rd <- read_markers(file.path(rer,  "results","markers", paste0(cc, ".allmarkers.csv")))
  o.cl <- if (!is.null(od)) length(unique(od$cluster)) else NA
  r.cl <- if (!is.null(rd)) length(unique(rd$cluster)) else NA
  o.g <- sig_genes(od); r.g <- sig_genes(rd)
  ot <- read_markers(file.path(orig, "results","markers", paste0(cc, ".cluster.top30.markers.csv")))
  rt <- read_markers(file.path(rer,  "results","markers", paste0(cc, ".cluster.top30.markers.csv")))
  o.tg <- if (!is.null(ot)) unique(ot$gene) else character(0)
  r.tg <- if (!is.null(rt)) unique(rt$gene) else character(0)
  mk_rows[[cc]] <- data.frame(condition=cc,
    clusters_original=o.cl, clusters_rerun=r.cl,
    sigGenes_original=length(o.g), sigGenes_rerun=length(r.g),
    sigGene_jaccard=jacc(o.g,r.g),
    top30_overlap=length(intersect(o.tg,r.tg)),
    top30_union=length(union(o.tg,r.tg)))
  cat(sprintf("  %-5s  clusters %2s->%2s | sigGenes %4d->%4d (Jacc %s) | top30 overlap %d/%d\n",
      cc, o.cl, r.cl, length(o.g), length(r.g), jacc(o.g,r.g),
      length(intersect(o.tg,r.tg)), length(union(o.tg,r.tg))))
}
markers <- do.call(rbind, mk_rows)

cat("\n==================================================================\n")
cat("3) NEW pseudobulk DESeq2 (rerun only) -- cardiac KO vs WT\n")
cat("==================================================================\n")
for (ts in c("P0","P7")) {
  f <- file.path(rer, "results","markers", paste0(ts, ".cardiac.pseudobulk.DESeq2.csv"))
  if (!file.exists(f)) { cat("  ", ts, "pseudobulk: MISSING\n"); next }
  d <- read.csv(f); names(d)[1] <- "gene"
  d <- d[!is.na(d$padj), ]
  nsig05 <- sum(d$padj < 0.05); nsig10 <- sum(d$padj < 0.10)
  cat(sprintf("\n  %s: %d genes tested, %d sig (padj<0.05), %d (padj<0.10)\n",
              ts, nrow(d), nsig05, nsig10))
  top <- head(d[order(d$padj), c("gene","log2FoldChange","pvalue","padj")], 15)
  top$log2FoldChange <- round(top$log2FoldChange,2)
  top$pvalue <- signif(top$pvalue,3); top$padj <- signif(top$padj,3)
  print(top, row.names = FALSE)
}

outdir <- file.path(our, "06_outputs")
write.csv(cells,   file.path(outdir, "comparison_cells.csv"),   row.names = FALSE)
write.csv(markers, file.path(outdir, "comparison_markers.csv"), row.names = FALSE)
cat("\nWrote 06_outputs/comparison_cells.csv and 06_outputs/comparison_markers.csv\n")
cat("=== COMPARE DONE ===\n")
