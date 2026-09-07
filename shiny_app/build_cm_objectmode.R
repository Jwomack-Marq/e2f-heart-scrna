# build_cm_objectmode.R
# ---------------------------------------------------------------------------
# Carry the CM object-mode diagnostic into app_data.rds as app$cmtest, for the
# "Object-mode test" tab in the cardiomyocyte deep-dive.
#
# The diagnostic itself is our_analysis/04_integrate_annotate/cm_objectmode_check.R.
# It builds the SAME cardiomyocytes three ways -- filter / roundtrip / rebuilt -- and
# measures whether the embedding changes. This script does no computation of its own;
# it only reads the two CSVs that script writes:
#
#   results/tables/cm_objectmode_percell.csv.gz    cell, variant, UMAP1, UMAP2, clusters
#   results/tables/cm_objectmode_comparison.csv    pair, metric, value, note
#
# So it needs no Seurat and runs anywhere the bundle does.
#
#   Rscript shiny_app/build_cm_objectmode.R [--tables=<dir>] [--max-cells=N]
# ---------------------------------------------------------------------------

# The CM compartment is 42,416 cells, so this does NOT fire at the default -- and that is
# deliberate. build_app_data.R downsamples app$cm$meta to 30k to keep the interactive
# plotly map responsive, but this tab renders a static ggplot server-side, where the extra
# cells cost a few MB in the bundle and nothing at render time. Keeping every cell means
# the three panels are the whole comparison rather than a sample of it. Lower it with
# --max-cells if bundle size ever becomes the binding constraint.
MAXCELLS <- 50000L

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

TABLES <- argval("--tables", "../our_analysis/results/tables")
MAXCELLS <- as.integer(argval("--max-cells", MAXCELLS))
f_pc  <- file.path(TABLES, "cm_objectmode_percell.csv.gz")
f_cmp <- file.path(TABLES, "cm_objectmode_comparison.csv")
if (!file.exists(f_pc) || !file.exists(f_cmp))
  stop("diagnostic output not found under ", TABLES,
       "\n  run our_analysis/04_integrate_annotate/cm_objectmode_check.R first")

cat("Reading diagnostic tables ...\n")
pc  <- read.csv(f_pc,  stringsAsFactors = FALSE, check.names = FALSE)
cmp <- read.csv(f_cmp, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("cell","variant","UMAP1","UMAP2") %in% names(pc)),
          all(c("pair","metric","value") %in% names(cmp)))

variants <- unique(pc$variant)
res_cols <- grep("^SCT_snn_res\\.", names(pc), value = TRUE)
cat(sprintf("  variants   : %s\n", paste(variants, collapse = ", ")))
cat(sprintf("  resolutions: %s\n", paste(sub("^SCT_snn_res\\.", "", res_cols), collapse = ", ")))

# Every variant must cover the same cells -- the whole comparison rests on that, and a
# partial re-run that silently dropped some would otherwise look like a real difference.
by_var <- split(pc$cell, pc$variant)
stopifnot(length(unique(vapply(by_var, length, integer(1)))) == 1L)
for (v in variants[-1]) stopifnot(setequal(by_var[[1]], by_var[[v]]))
cells <- by_var[[1]]
cat(sprintf("  cells      : %s per variant\n", format(length(cells), big.mark = ",")))

if (length(cells) > MAXCELLS) {          # same cells across variants, so panels stay comparable
  set.seed(42); keep <- sample(cells, MAXCELLS)
  pc <- pc[pc$cell %in% keep, , drop = FALSE]
  cat(sprintf("  downsampled to %s cells per variant\n", format(MAXCELLS, big.mark = ",")))
}
for (rc in res_cols) pc[[rc]] <- as.character(pc[[rc]])

# The headline the tab shows without making anyone read the table. A vs B is the question
# that was actually asked; report it as a yes/no and keep the number that justifies it.
getv <- function(pair, metric) {
  v <- cmp$value[cmp$pair == pair & cmp$metric == metric]
  if (length(v)) v[1] else NA_real_
}
ab_diff <- getv("filter vs roundtrip", "umap_max_abs_diff")
ab_ari  <- getv("filter vs roundtrip", "ARI_SCT_snn_res.0.2")
ac_diff <- getv("filter vs rebuilt",   "umap_procrustes_rmsd")
ac_ari  <- getv("filter vs rebuilt",   "ARI_SCT_snn_res.0.2")
ab_same <- isTRUE(!is.na(ab_diff) && ab_diff == 0 && !is.na(ab_ari) && abs(ab_ari - 1) < 1e-9)

verdict <- if (ab_same) {
  paste0("Saving and reloading the subset changes nothing: identical UMAP coordinates ",
         "(max |difference| = 0) and identical clusters (ARI = 1). A Seurat subset holds ",
         "no reference to its parent, so there was no link for a round-trip to drop.")
} else {
  paste0("Round-trip did NOT reproduce the filtered object exactly (max |UMAP difference| = ",
         signif(ab_diff, 4), ", ARI = ", signif(ab_ari, 6), "). That points at something ",
         "non-deterministic in the pipeline and should be chased before reading the ",
         "rebuilt-object comparison.")
}

app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_objmode.bak.rds\n")
saveRDS(app, "app_data.pre_objmode.bak.rds", compress = "gzip")

app$cmtest <- list(
  percell  = pc,
  metrics  = cmp,
  variants = variants,
  res      = sub("^SCT_snn_res\\.", "", res_cols),
  ab_same  = ab_same,
  verdict  = verdict,
  labels   = c(filter    = "A — filter (production)",
               roundtrip = "B — round-trip (save + reload)",
               rebuilt   = "C — rebuilt (own object)"),
  blurb    = c(
    filter    = "cm <- comb[, comb$celltype == \"Cardiomyocyte\"] — exactly what cm_subcluster_build.R does.",
    roundtrip = "The same subset written to disk with saveRDS, dropped from memory, and read back before re-embedding.",
    rebuilt   = "CreateSeuratObject on the cardiomyocyte counts with zero-count genes dropped: no parent assays, reductions, graphs or SCT model carried over."),
  note     = paste("All three arms run the identical SCTransform -> PCA -> Harmony -> UMAP",
                   "-> clusters pipeline with the same seeds at every stage, on cells in the",
                   "same order. Rows tagged CONFOUNDED compare against the shipped production",
                   "object, which was built with different package versions -- a difference",
                   "there says nothing about how the object was constructed."))

cat("\n== Summary ==\n")
cat(sprintf("  per-cell rows : %s\n", format(nrow(pc), big.mark = ",")))
cat(sprintf("  metric rows   : %s\n", format(nrow(cmp), big.mark = ",")))
cat(sprintf("  A vs B identical : %s\n", if (ab_same) "YES" else "NO"))
cat(sprintf("  A vs C procrustes: %s   ARI(res 0.2): %s\n",
            signif(ac_diff, 4), signif(ac_ari, 6)))
cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
