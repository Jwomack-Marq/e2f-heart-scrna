#!/usr/bin/env Rscript
# Phase E -- recall doublets with scDblFinder and compare to the existing Scrublet
# calls (predicted_doublets.csv). The README suspects Scrublet under-called
# (1-2% vs the 8% expected_doublet_rate). One isolated process per lane is cheap.
# Reads 01_input/<lane>/filtered_matrix/sensitivity_5/ + 01_input/<lane>/predicted_doublets.csv;
# writes results/tables/doublet_comparison.csv.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(scDblFinder); library(SingleCellExperiment) }))

lanes <- c("P0KO_lane1","P0KO_lane6","P0WT_lane1","P0WT_lane6",
           "P7KO_lane1","P7KO_lane6","P7WT_lane1","P7WT_lane6")
rows <- list()
for (s in lanes) {
  mtx <- file.path(INPUT, s, "filtered_matrix", "sensitivity_5")
  counts <- Read10X(data.dir = mtx)
  sce <- SingleCellExperiment(assays = list(counts = counts))
  set.seed(1)
  sce <- scDblFinder(sce)
  scd_rate <- mean(sce$scDblFinder.class == "doublet")

  scr <- read.csv(file.path(INPUT, s, "predicted_doublets.csv"))
  scr_rate <- mean(as.logical(as.character(scr$predicted_doublets)), na.rm = TRUE)

  rows[[s]] <- data.frame(lane = s, n_cells = ncol(sce),
                          scrublet_pct = round(100 * scr_rate, 2),
                          scDblFinder_pct = round(100 * scd_rate, 2))
  cat(sprintf("  %-12s scrublet=%.2f%%  scDblFinder=%.2f%%\n", s,
              100*scr_rate, 100*scd_rate))
  rm(counts, sce); gc(verbose = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUTTAB, "doublet_comparison.csv"), row.names = FALSE)
cat("\n--- doublet rate comparison ---\n"); print(res, row.names = FALSE)
cat("=== DONE qc_doublets ===\n")
