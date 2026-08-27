#!/usr/bin/env Rscript
# Runner for steps 02-03 (lane merges + condition merges). Pandoc isn't required:
# instead of rendering the .Rmd reports we extract their code with knitr::purl and
# source it. Reads raw matrices from 01_input/ and writes into our_analysis/
# (processing/ + results/).
#
# Usage (run from anywhere; paths are anchored on the our_analysis/.projroot file):
#   Rscript our_analysis/run_pipeline.R lane  P0WT   # step 02: one lane-merge group
#   Rscript our_analysis/run_pipeline.R merge P0     # step 03: one condition merge (P0/P7)

suppressWarnings(suppressMessages(library(knitr)))

args <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
key   <- args[2]

# our_analysis root = the folder containing the .projroot sentinel. Walk up from
# this script's location so it works regardless of the caller's working directory.
a0 <- commandArgs(FALSE); sp <- sub("^--file=", "", a0[grep("^--file=", a0)])
our_root <- normalizePath(if (length(sp) && nzchar(sp)) dirname(sp) else getwd())
while (!file.exists(file.path(our_root, ".projroot"))) {
  up <- dirname(our_root)
  if (identical(up, our_root)) stop("our_analysis/.projroot not found above ", our_root)
  our_root <- up
}
setwd(our_root)
Sys.setenv(SCRNA_OUT_DIR = ".")          # outputs -> our_analysis/{processing,results}
Sys.setenv(SCRNA_IN_DIR  = "01_input")   # raw PIPseeker matrices (copied into our_analysis)

# Run future sequentially so Seurat's FindMarkers/SCTransform never try to ship
# multi-GB globals to parallel workers (the cause of the maxSizeOfObjects errors).
suppressWarnings(suppressMessages(library(future)))
future::plan("sequential")
# 3 GiB safety brake (NOT Inf -- Inf is what caused the original crash). The Rmd
# setup chunks set the same value; this keeps the runner consistent with it.
options(future.globals.maxSize = 3 * 1024^3)

if (!dir.exists("logs")) dir.create("logs")
run_rmd <- function(rmd) {
  rfile <- tempfile(fileext = ".R")
  knitr::purl(rmd, output = rfile, quiet = TRUE, documentation = 0)
  pdf(file.path("logs", "throwaway_plots.pdf"))   # capture stray on-screen plots
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  source(rfile, local = new.env())
  invisible(NULL)
}

t0 <- Sys.time()
if (identical(stage, "lane")) {
  Sys.setenv(SCRNA_SAMPLE = key)
  cat("=== lane_merge:", key, "===\n")
  run_rmd("02_qc_lane_merge/scRNA_lane_merge.Rmd")
} else if (identical(stage, "merge")) {
  rmd <- if (key == "P0") "03_condition_merge/scRNA_mergeP0.Rmd"
         else            "03_condition_merge/scRNA_mergeP7.Rmd"
  cat("=== merge:", key, "===\n")
  run_rmd(rmd)
} else {
  stop("Unknown stage: ", stage)
}
cat(sprintf("=== DONE %s %s in %.1f min ===\n", stage, key,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
