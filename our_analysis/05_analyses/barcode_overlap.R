#!/usr/bin/env Rscript
# Check, per library, how many cell barcodes are shared between lane1 and lane6
# (same library sequenced on two lanes). High overlap => the lane-merge step
# double-counts those cells. Reads only barcodes.tsv.gz; fast.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))

libs <- c("P0WT","P0KO","P7WT","P7KO")
read_bc <- function(lib, lane)
  readLines(gzfile(file.path(INPUT, paste0(lib, "_", lane),
                             "filtered_matrix", "sensitivity_5", "barcodes.tsv.gz")))
rows <- list()
for (lib in libs) {
  b1 <- read_bc(lib, "lane1"); b6 <- read_bc(lib, "lane6")
  inter <- length(intersect(b1, b6)); uni <- length(union(b1, b6))
  rows[[lib]] <- data.frame(library = lib, lane1 = length(b1), lane6 = length(b6),
    shared = inter, union = uni,
    pct_lane1_shared = round(100*inter/length(b1), 1),
    pct_lane6_shared = round(100*inter/length(b6), 1),
    jaccard = round(inter/uni, 3))
  cat(sprintf("  %-5s lane1=%d lane6=%d shared=%d (%.1f%% of lane1, %.1f%% of lane6)\n",
              lib, length(b1), length(b6), inter, 100*inter/length(b1), 100*inter/length(b6)))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUTTAB, "lane_barcode_overlap.csv"), row.names = FALSE)
cat("\n"); print(res, row.names = FALSE)
cat("=== DONE barcode_overlap ===\n")
