# build_pcdims.R
# ---------------------------------------------------------------------------
# Carry the PC-dimension sweep into app_data.rds as app$pcdims, for the
# "PC dimensions" tab.
#
# The sweep itself is our_analysis/04_integrate_annotate/pcdims_sweep.R, run once per
# object. It holds SCTransform/PCA/Harmony fixed and varies only the dims cut carried
# into FindNeighbors and RunUMAP -- 10, 30, 50, where 30 is what both production
# embeddings use. This script only reads the CSVs it wrote, so it needs no Seurat.
#
#   results/tables/pcdims_<obj>_percell.csv.gz     cell, dims, UMAP1, UMAP2, clusters, meta
#   results/tables/pcdims_<obj>_comparison.csv     pair, metric, value, note
#
#   Rscript shiny_app/build_pcdims.R [--tables=<dir>] [--objects=cm,combined]
# ---------------------------------------------------------------------------

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

TABLES  <- argval("--tables", "../our_analysis/results/tables")
OBJECTS <- strsplit(argval("--objects", "cm,combined"), ",")[[1]]
LABELS  <- c(cm = "Cardiomyocytes", combined = "Whole heart")
COLBY   <- list(                     # what the map can be coloured by, per object
  cm       = c("Subcluster" = "cluster", "Genotype" = "genotype",
               "Timepoint" = "timepoint", "Cell-cycle phase" = "Phase"),
  combined = c("Cell type" = "celltype", "Cluster" = "cluster",
               "Genotype" = "genotype", "Timepoint" = "timepoint"))

one <- function(obj) {
  f_pc  <- file.path(TABLES, sprintf("pcdims_%s_percell.csv.gz", obj))
  f_cmp <- file.path(TABLES, sprintf("pcdims_%s_comparison.csv", obj))
  if (!file.exists(f_pc) || !file.exists(f_cmp)) {
    cat(sprintf("  %-9s SKIPPED (no sweep output under %s)\n", obj, TABLES)); return(NULL)
  }
  pc  <- read.csv(f_pc,  stringsAsFactors = FALSE, check.names = FALSE)
  cmp <- read.csv(f_cmp, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("cell","dims","UMAP1","UMAP2") %in% names(pc)),
            all(c("pair","metric","value") %in% names(cmp)))

  dims <- sort(unique(pc$dims))
  # Every dims arm must cover the same cells; that is the whole basis of the comparison,
  # and a partial re-run that dropped some would otherwise read as a real difference.
  by_d <- split(pc$cell, pc$dims)
  stopifnot(length(unique(vapply(by_d, length, integer(1)))) == 1L)
  for (i in seq_along(by_d)[-1]) stopifnot(setequal(by_d[[1]], by_d[[i]]))

  res_cols <- grep("^SCT_snn_res\\.", names(pc), value = TRUE)
  for (rc in res_cols) pc[[rc]] <- as.character(pc[[rc]])
  # The app colours by a single "cluster" column; pick the production resolution so the
  # tab shows the clusters people actually look at (app.R exposes res 0.2 for CM;
  # combined.R clusters once at 0.8).
  prod_rc <- if (obj == "cm") "SCT_snn_res.0.2" else "SCT_snn_res.0.8"
  if (!prod_rc %in% res_cols) prod_rc <- res_cols[1]
  pc$cluster <- pc[[prod_rc]]

  getv <- function(pair, metric) {
    v <- cmp$value[cmp$pair == pair & cmp$metric == metric]
    if (length(v)) v[1] else NA_real_
  }
  varpct <- c("10" = getv("(all)", "varpct_10"), "30" = getv("(all)", "varpct_30"),
              "50" = getv("(all)", "varpct_50"))
  nclust <- vapply(dims, function(d)
    getv(paste0("dims", d), paste0("nclust_", prod_rc)), numeric(1))
  names(nclust) <- as.character(dims)

  colby <- COLBY[[obj]]
  colby <- colby[colby %in% c("cluster", names(pc))]

  cat(sprintf("  %-9s %s cells x %d dims arms | %s | clusters %s\n", obj,
              format(length(by_d[[1]]), big.mark = ","), length(dims), prod_rc,
              paste(sprintf("%d:%g", dims, nclust), collapse = " ")))

  list(label = unname(LABELS[obj]),
       # CM subclusters are written CM0..CM10 everywhere else in the app; whole-heart
       # clusters have no such convention, so they stay bare.
       prefix = if (obj == "cm") "CM" else "C",
       percell = pc, metrics = cmp, dims = dims,
       res = sub("^SCT_snn_res\\.", "", res_cols), prod_res = prod_rc,
       varpct = varpct, nclust = nclust, colby = colby,
       # m2 and knn30, not umap_procrustes_rmsd: that one is n-dependent (two unrelated
       # embeddings score sqrt(2/n), not 1) and reads as "almost identical" when it is
       # nothing of the sort. See pcdims_shape_metrics.R.
       m2    = c("10v30" = getv("dims 10 vs 30", "procrustes_m2"),
                 "30v50" = getv("dims 30 vs 50", "procrustes_m2"),
                 "10v50" = getv("dims 10 vs 50", "procrustes_m2")),
       knn30 = c("10v30" = getv("dims 10 vs 30", "knn30_kept"),
                 "30v50" = getv("dims 30 vs 50", "knn30_kept"),
                 "10v50" = getv("dims 10 vs 50", "knn30_kept")),
       ari = c("10v30" = getv("dims 10 vs 30", paste0("ARI_", prod_rc)),
               "30v50" = getv("dims 30 vs 50", paste0("ARI_", prod_rc)),
               "10v50" = getv("dims 10 vs 50", paste0("ARI_", prod_rc))))
}

cat("Reading sweep tables ...\n")
got <- Filter(Negate(is.null), setNames(lapply(OBJECTS, one), OBJECTS))
if (!length(got)) stop("no sweep output found under ", TABLES,
                       "\n  run our_analysis/04_integrate_annotate/pcdims_sweep.R first")

app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_pcdims.bak.rds\n")
saveRDS(app, "app_data.pre_pcdims.bak.rds", compress = "gzip")

app$pcdims <- c(got, list(
  objects = names(got),
  note = paste("SCTransform, PCA (npcs = 50) and Harmony are computed ONCE and shared by",
               "all three arms -- only the dims cut carried into FindNeighbors and RunUMAP",
               "varies. Harmony corrects every PCA dimension and production then slices the",
               "corrected embedding, so \"carrying N PCs\" means taking the first N harmony",
               "dims, not re-running Harmony on a shorter PCA. Both production embeddings",
               "use dims 1:30 (combined.R:36-38, cm_subcluster_build.R:39-40); neither",
               "records why, which is what this sweep is for.")))

cat("\n== Summary ==\n")
for (o in names(got)) {
  g <- got[[o]]
  cat(sprintf("  %-9s variance 10/30/50 PCs: %.1f%% / %.1f%% / %.1f%%\n",
              o, g$varpct[["10"]], g$varpct[["30"]], g$varpct[["50"]]))
  cat(sprintf("  %-9s 30NN kept 10v30 %.1f%% | 30v50 %.1f%%   ARI 10v30 %.3f | 30v50 %.3f\n",
              "", 100 * g$knn30[["10v30"]], 100 * g$knn30[["30v50"]],
              g$ari[["10v30"]], g$ari[["30v50"]]))
}
cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
