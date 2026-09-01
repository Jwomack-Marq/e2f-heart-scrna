# build_clusterings.R
# ---------------------------------------------------------------------------
# Carry the clustering registry and its per-variant downstream results into
# app_data.rds as app$clusterings, so the CM panels can be driven by a chosen
# labelling instead of only the one that happens to be production.
#
# Reads only CSVs, so it needs no Seurat and runs in the app image.
#   results/tables/clustering_registry.csv              (pcdims_sweep.R)
#   results/tables/pcdims_cm_percell.csv.gz             (pcdims_sweep.R) -- labels
#   results/tables/cm_subcluster_<tag>_*.csv            (cm_subcluster_analyze.R)
#
# Each variant carries a bare cell -> cluster label vector, not a copy of the per-cell
# metadata: genotype, timepoint and Phase are identical across variants by construction,
# and the app already holds them once in app$meta.
#
# Labels are trimmed to cells the bundle can actually compute on. The sweep covers all
# 42,416 CMs but only 21,598 are in app$expr and 5,755 in app$deg_expr, so the rest cannot
# feed a marker test or a score no matter what the app does with them.
#
#   Rscript shiny_app/build_clusterings.R [--tables=<dir>]
# ---------------------------------------------------------------------------

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))
TABLES <- argval("--tables", "../our_analysis/results/tables")

regf <- file.path(TABLES, "clustering_registry.csv")
if (!file.exists(regf))
  stop("no clustering registry at ", regf,
       "\n  run our_analysis/04_integrate_annotate/pcdims_sweep.R --object=cm first")
reg <- read.csv(regf, stringsAsFactors = FALSE)
reg <- reg[reg$object == "cm", , drop = FALSE]          # only CM has a downstream pipeline
stopifnot(nrow(reg) > 0)

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_clusterings.bak.rds\n")
saveRDS(app, "app_data.pre_clusterings.bak.rds", compress = "gzip")

# Cells the bundle can compute on at all. app$expr is the curated panel over all CM cells
# in the bundle; deg_expr is the broad matrix over a smaller subset.
bundle_cells <- colnames(app$expr)
cat(sprintf("bundle cells: %s (broad matrix: %s)\n",
            format(length(bundle_cells), big.mark = ","),
            format(ncol(app$deg_expr), big.mark = ",")))

# ---- labels, from the sweep's per-cell table -------------------------------
pcf <- file.path(TABLES, "pcdims_cm_percell.csv.gz")
if (!file.exists(pcf)) stop("no per-cell labels at ", pcf)
pc <- read.csv(pcf, stringsAsFactors = FALSE, check.names = FALSE)

read_opt <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL

variants <- list()
for (i in seq_len(nrow(reg))) {
  vid <- reg$variant_id[i]
  tag <- sub("^cm_", "", vid)                       # analyze step's output namespace
  rescol <- paste0("SCT_snn_res.", reg$resolution[i])
  sl <- pc[pc$dims == reg$dims[i], , drop = FALSE]
  if (!rescol %in% names(sl)) { cat(sprintf("  %-22s SKIP (no %s)\n", vid, rescol)); next }

  lab <- setNames(as.character(sl[[rescol]]), sl$cell)
  lab <- lab[names(lab) %in% bundle_cells]

  # Per-subcluster KO-vs-WT DE: one file per cluster, keyed CM0..CMn.
  de_files <- list.files(TABLES, full.names = TRUE,
    pattern = sprintf("^cm_subcluster_%s_KOvsWT_CM[0-9]+\\.descriptive\\.DE\\.csv$",
                      gsub("\\.", "\\\\.", tag)))
  de <- setNames(lapply(de_files, read.csv, stringsAsFactors = FALSE),
                 sub(".*_KOvsWT_(CM[0-9]+)\\..*", "\\1", basename(de_files)))
  de <- de[order(as.integer(sub("CM", "", names(de))))]

  variants[[vid]] <- list(
    variant_id  = vid, dims = reg$dims[i], resolution = reg$resolution[i],
    label       = sprintf("dims 1:%d, res %s", reg$dims[i], reg$resolution[i]),
    is_production = isTRUE(as.logical(reg$is_production[i])),
    n_clusters  = reg$n_clusters[i],
    labels      = lab,
    de          = if (length(de)) de else NULL,
    de_summary  = read_opt(file.path(TABLES, sprintf("cm_subcluster_%s_KOvsWT_summary.csv", tag))),
    markers     = read_opt(file.path(TABLES, sprintf("cm_subcluster_top_markers_%s.csv", tag))),
    subtype_map = read_opt(file.path(TABLES, sprintf("cm_subcluster_subtype_map_%s.csv", tag))),
    cellcycle   = read_opt(file.path(TABLES, sprintf("cm_subcluster_%s_cellcycle.csv", tag))))
  cat(sprintf("  %-22s %2d clusters | %s labelled cells | DE %s | markers %s\n",
              vid, reg$n_clusters[i], format(length(lab), big.mark = ","),
              if (length(de)) sprintf("%d clusters", length(de)) else "MISSING",
              if (is.null(variants[[vid]]$markers)) "MISSING" else "ok"))
}
stopifnot(length(variants) > 0)

is_prod <- as.logical(reg$is_production); is_prod[is.na(is_prod)] <- FALSE
prod_id <- reg$variant_id[is_prod]
if (length(prod_id) != 1L)
  stop("registry must flag exactly one production variant; found ", length(prod_id),
       " (", paste(prod_id, collapse = ", "), ")")
if (!prod_id %in% names(variants))
  stop("production variant ", prod_id, " has no ingested results -- run ",
       "cm_subcluster_analyze.R --variant=", prod_id, " first")

for (vid in names(variants))
  variants[[vid]]$has_downstream <- !is.null(variants[[vid]]$de)

app$clusterings <- list(
  registry   = reg,
  production = prod_id,
  variants   = variants,
  note = paste("Every variant runs the identical SCTransform / PCA / Harmony; only the",
               "number of PCs carried into FindNeighbors and RunUMAP differs. Downstream",
               "results are keyed to the labelling that produced them, so a KO-vs-WT",
               "number read here belongs to the selected variant, not to production."))

# NOTE on a thing deliberately NOT done: app$pcdims stores per-cell metadata once per dims
# arm, and de-duplicating it looked like an easy ~20 MB win. Measured, it went 33.3 -> 35.7
# MB: splitting one data.frame into two copies the cell-barcode column, and R's global
# string cache already shares the repeated metadata strings that object.size counts three
# times over. The saving was an artefact of the measurement, so the code is not here.

cat(sprintf("\n== Summary ==\n  variants   : %d\n  production : %s\n",
            length(variants), prod_id))
cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
