# build_immune_flag.R
# ---------------------------------------------------------------------------
# Carry the immune-contamination flag into app_data.rds so the app can keep those cells
# off the statistics they distort.
#
# The CM compartment contains leukocytes that annotate.R labelled "Cardiomyocyte" (see
# docs/05-cm-deepdive.qmd#cm12). This does NOT re-annotate or remove anything -- the real
# fix is in annotate.R and belongs with the regeneration the lane double-counting already
# requires. It marks them so a cycling fraction computed across CM subclusters is not
# quietly inflated by immune cells, which genuinely cycle.
#
# Why the flag is per-cell and score-based rather than "drop CM12": the measured
# concordance is that the isolating cluster holds only 55-63 % of flagged cells at
# dims 30/50, and at dims 10 no cluster is more than half flagged at all -- the rest are
# dispersed through otherwise-real CM subclusters. Dropping a cluster id would miss them.
#
#   Reads  results/tables/cm_immune_contamination.csv   (cm_immune_contamination.R)
#   Writes app$cm$meta$immune_contam + app$immune_contam summary
#
#   Rscript shiny_app/build_immune_flag.R [--tables=<dir>]
# ---------------------------------------------------------------------------

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))
TABLES <- argval("--tables", "../our_analysis/results/tables")

f <- file.path(TABLES, "cm_immune_contamination.csv")
if (!file.exists(f))
  stop("no flag table at ", f,
       "\n  run our_analysis/05_analyses/cm_immune_contamination.R first")
fl <- read.csv(f, stringsAsFactors = FALSE)
stopifnot(all(c("cell", "immune_contam", "n_immune_genes", "n_mast_genes") %in% names(fl)))

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_immuneflag.bak.rds\n")
saveRDS(app, "app_data.pre_immuneflag.bak.rds", compress = "gzip")

lut <- setNames(as.logical(fl$immune_contam), fl$cell)

# Attach to the CM meta the app plots from. Cells absent from the flag table are FALSE
# rather than NA: the flag drives an exclusion, and an NA would propagate into every
# downstream mean rather than simply leaving the cell in.
cmm <- app$cm$meta
cmm$immune_contam <- unname(lut[cmm$cell]); cmm$immune_contam[is.na(cmm$immune_contam)] <- FALSE
app$cm$meta <- cmm
n_cm <- sum(cmm$immune_contam)

# Same on the all-cell meta, so a CM-compartment statistic computed there agrees.
m <- app$meta
m$immune_contam <- unname(lut[m$cell]); m$immune_contam[is.na(m$immune_contam)] <- FALSE
app$meta <- m

# Cycling is the statistic this actually protects, so quantify the correction rather than
# asserting it matters.
cyc <- function(d) if (!"cycling" %in% names(d)) NA_real_ else 100 * mean(as.logical(d$cycling), na.rm = TRUE)
with_all  <- cyc(cmm)
with_none <- cyc(cmm[!cmm$immune_contam, , drop = FALSE])
only_flag <- cyc(cmm[cmm$immune_contam, , drop = FALSE])

app$immune_contam <- list(
  n_flagged_bundle = n_cm,
  n_flagged_total  = sum(fl$immune_contam),
  n_total          = nrow(fl),
  pct_mast         = 100 * mean(fl$n_mast_genes[fl$immune_contam] >= 2),
  cycling_all      = with_all,
  cycling_excl     = with_none,
  cycling_flagged  = only_flag,
  note = paste(
    "Cells labelled Cardiomyocyte by annotate.R that carry pan-leukocyte markers",
    "(>= 2 of Ptprc, Cd52, Laptm5, Coro1a, Cd3e, Trbc1, Cd79a). They are not doublets --",
    "Scrublet clears them, and they carry roughly half the RNA and a third the",
    "mitochondrial content of a cardiomyocyte. They were mislabelled because annotate.R",
    "takes max.col over nine marker panels per whole-heart cluster and none of the nine",
    "describes a mast cell or a lymphocyte, so ambient cardiomyocyte RNA wins by default.",
    "Flagged rather than removed: re-annotating is the real fix and belongs with the",
    "regeneration the lane double-counting already requires."),
  caveat = paste(
    "The flag is deliberately broad and errs toward exclusion. It is not a re-annotation:",
    "some dispersed flagged cells may be genuine cardiomyocytes carrying ambient immune",
    "RNA. At 0.38 % of the compartment, over-excluding a few costs nothing, while leaving",
    "immune cells in a cycling fraction does not."))

cat(sprintf("\n== Summary ==\n"))
cat(sprintf("  flagged: %d of %d CM cells in the bundle (%s of %s genome-wide, %.2f%%)\n",
            n_cm, nrow(cmm), format(sum(fl$immune_contam), big.mark = ","),
            format(nrow(fl), big.mark = ","), 100 * mean(fl$immune_contam)))
cat(sprintf("  %.0f%% of flagged cells are mast-marker positive\n", app$immune_contam$pct_mast))
cat(sprintf("  cycling fraction: all %.1f%% | excluding flagged %.1f%% | flagged themselves %.1f%%\n",
            with_all, with_none, only_flag))
cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
