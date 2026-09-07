#!/usr/bin/env Rscript
# Inject the QC / normalization figures (as base64 data-URIs) + the doublet-rate
# table into an EXISTING app_data.rds, without re-loading the multi-GB Seurat
# objects. (The same embedding is also done in build_app_data.R for fresh builds.)
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))   # OUTFIG, OUTTAB, OUT
suppressWarnings(suppressMessages(library(base64enc)))

fig_uri <- function(name) { p <- file.path(OUTFIG, name)
  if (file.exists(p)) base64enc::dataURI(file = p, mime = "image/png") else NA_character_ }

rds <- file.path(OUT, "app", "app_data.rds"); app <- readRDS(rds)
app$figs <- list(
  filtering  = fig_uri("norm_filtering_summary.png"),
  qc_violins = fig_uri("norm_qc_violins.png"),
  doublet    = fig_uri("QC_doublet_comparison.png"),
  hvg        = fig_uri("norm_hvg.png"),
  harmony    = fig_uri("combined_harmony_before_after.png"))
dcsv <- file.path(OUTTAB, "doublet_comparison.csv")
if (file.exists(dcsv)) app$tables$doublet <- read.csv(dcsv, check.names = FALSE)
saveRDS(app, rds, compress = "gzip")
cat(sprintf("Patched %s (%.1f MB) | figs: %s\n", rds, file.info(rds)$size/1024^2,
            paste(names(Filter(function(x) !is.na(x), app$figs)), collapse = ", ")))
