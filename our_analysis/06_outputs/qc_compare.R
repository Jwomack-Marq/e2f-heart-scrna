#!/usr/bin/env Rscript
# QC comparison: original filter (nFeature >= 1500 only) vs the fixed filter
# (nFeature >= 1500 & <= per-lane 99.5th pctl & percent.mt <= 20 & not doublet).
# Quantifies how many cells the new QC removes, per lane. Needs only Seurat.

suppressMessages(library(Seurat))

# our_analysis root = parent of this script's folder (06_outputs/)
args0 <- commandArgs(trailingOnly = FALSE)
script.path <- sub("^--file=", "", args0[grep("^--file=", args0)])
our <- normalizePath(file.path(dirname(script.path), ".."))
input <- file.path(our, "01_input")
cat("our_analysis root:", our, "\n")

cutoff.nFeature <- 1500
cutoff.mt <- 20
upper.pctl <- 0.995

lanes <- c("P0KO_lane1","P0KO_lane6","P0WT_lane1","P0WT_lane6",
           "P7KO_lane1","P7KO_lane6","P7WT_lane1","P7WT_lane6")

rows <- list()
for (s in lanes) {
  mtx.dir <- file.path(input, s, "filtered_matrix", "sensitivity_5")
  counts <- Read10X(data.dir = mtx.dir)
  obj <- CreateSeuratObject(counts = counts)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

  # doublet flags (predicted_doublets.csv rows are in matrix-column order)
  dbl <- read.csv(file.path(input, s, "predicted_doublets.csv"))
  is.dbl <- as.logical(as.character(dbl$predicted_doublets))
  if (length(is.dbl) == ncol(obj)) obj$is_doublet <- is.dbl
  else { warning(sprintf("%s: doublet rows (%d) != cells (%d); skipping doublet filter",
                         s, length(is.dbl), ncol(obj))); obj$is_doublet <- FALSE }

  nf <- obj$nFeature_RNA
  mt <- obj$percent.mt
  db <- obj$is_doublet
  max.nf <- as.numeric(quantile(nf, upper.pctl))

  pass.min   <- nf >= cutoff.nFeature                 # original filter
  pass.max   <- nf <= max.nf
  pass.mt    <- mt <= cutoff.mt
  pass.dbl   <- !db
  pass.new   <- pass.min & pass.max & pass.mt & pass.dbl

  rows[[s]] <- data.frame(
    sample = s,
    barcodes = ncol(obj),
    kept_original = sum(pass.min),
    kept_new = sum(pass.new),
    removed_doublet = sum(pass.min & db),                       # doublets among original-kept
    removed_highmito = sum(pass.min & !pass.mt),                # high-mito among original-kept
    removed_uppercap = sum(pass.min & pass.max == FALSE),       # upper-cap among original-kept
    median_genes = round(median(nf)),
    median_mt = round(median(mt), 2),
    upper_nFeature_cap = round(max.nf)
  )
  cat(sprintf("%-12s barcodes=%5d  orig_kept=%5d  new_kept=%5d  (-%d doublet, -%d highmito)\n",
              s, ncol(obj), sum(pass.min), sum(pass.new),
              sum(pass.min & db), sum(pass.min & !pass.mt)))
  rm(counts, obj); gc(verbose = FALSE)
}

res <- do.call(rbind, rows)
res$pct_dropped_vs_original <- round(100 * (res$kept_original - res$kept_new) / res$kept_original, 1)
out <- file.path(our, "06_outputs", "qc_comparison.csv")
write.csv(res, out, row.names = FALSE)
cat("\nWrote:", out, "\n")
print(res, row.names = FALSE)
