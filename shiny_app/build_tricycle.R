# build_tricycle.R
# ---------------------------------------------------------------------------
# Carry the tricycle cell-cycle projection into app_data.rds.
#
# WHY. 05_analyses/cellcycle_tricycle.R gives every cell a CONTINUOUS position on the cell
# cycle -- an angle theta, plus the 2D cycle space it was measured in -- and an independent
# discrete stage call that is allowed to abstain. Nothing consumed any of it: the app knew
# only Seurat's three-way Phase label, which is an argmax and therefore cannot say "between
# S and G2" or "I can't tell". A position is the thing you can actually draw a cell cycle
# with, so the app needs it as data, not as two PNGs in a git-ignored folder.
#
# The app reads ONLY app_data.rds -- no CSV is opened at runtime -- so this is the step that
# makes the analysis reachable from the browser at all.
#
#   Reads  results/tables/cellcycle_tricycle_{percell,vs_seurat_by_celltype,confusion,
#                                             marker_peaks,depthmatched_cm,controls}.csv
#          results/tables/cellcycle_tricycle_matched_genotype.csv  (optional; from
#            cellcycle_tricycle_depthmatched.R -- the KO-vs-WT gap re-measured after
#            binomially thinning both genotypes to one depth distribution)
#   Writes app$tricycle           (all 58,917 cells + the summary tables)
#          app$meta$tricycle_*    (3 columns, joined onto the ~30k downsample)
#
#   Rscript shiny_app/build_tricycle.R [--tables=<dir>]
# ---------------------------------------------------------------------------

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))
TABLES <- argval("--tables", "../our_analysis/results/tables")

rd <- function(name, required = TRUE) {
  f <- file.path(TABLES, paste0("cellcycle_tricycle_", name, ".csv"))
  if (!file.exists(f)) {
    if (required) stop("missing ", f, "\n  run our_analysis/05_analyses/cellcycle_tricycle.R first")
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

percell  <- rd("percell")
by_ct    <- rd("vs_seurat_by_celltype")
confus   <- rd("confusion")
peaks    <- rd("marker_peaks")
depthm   <- rd("depthmatched_cm")
controls <- rd("controls")
# Optional: an un-run depth-matching step should leave the rest of the tab working.
matched_gt <- rd("matched_genotype", required = FALSE)   # NB: `matched` is taken below

# The embedding is the reason this build exists; without it the app can only draw cells on a
# ring of constant radius, which asserts a confidence the measurement does not have.
need <- c("tricycle_pc1", "tricycle_pc2", "ngene")
if (!all(need %in% names(percell)))
  stop("cellcycle_tricycle_percell.csv predates the cycle-space columns (missing: ",
       paste(setdiff(need, names(percell)), collapse = ", "),
       ")\n  re-run our_analysis/05_analyses/cellcycle_tricycle.R")

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_tricycle.bak.rds\n")
saveRDS(app, "app_data.pre_tricycle.bak.rds", compress = "gzip")

## ---- the per-cell frame ----------------------------------------------------
# Abstentions are given a LABEL rather than left as NA. NA reads as "missing data" and gets
# silently dropped or greyed by most of the app's plotting; the whole point of the Schwabe
# call is that refusing to stage a cell is an answer, so it has to survive as a category.
NOT_STAGED <- "(not staged)"
percell$tricycle_stage[is.na(percell$tricycle_stage) | !nzchar(percell$tricycle_stage)] <- NOT_STAGED

STAGE_LEVELS <- c("G1.S", "S", "G2", "G2.M", "M.G1", NOT_STAGED)   # cycle order, then abstain
stopifnot(all(percell$tricycle_stage %in% STAGE_LEVELS))

keep <- c("cell", "celltype", "timepoint", "genotype", "seurat_phase", "seurat_cycling",
          "tricycle_theta", "tricycle_pc1", "tricycle_pc2", "tricycle_stage",
          "tricycle_cycling", "ngene", "numi")
pc <- percell[, intersect(keep, names(percell)), drop = FALSE]
# Factors, not characters: 58,917 repeated strings across five columns is most of the size of
# this object, and the app needs a fixed level ORDER on stage and genotype anyway.
for (v in c("celltype", "timepoint", "seurat_phase")) pc[[v]] <- factor(pc[[v]])
pc$genotype       <- factor(pc$genotype, levels = c("WT", "KO"))   # WT first, as everywhere else
pc$tricycle_stage <- factor(pc$tricycle_stage, levels = STAGE_LEVELS)

app$tricycle <- list(
  percell      = pc,
  by_celltype  = by_ct,
  confusion    = confus,
  marker_peaks = peaks,
  depthmatched = depthm,
  controls     = controls,
  matched      = matched_gt,
  stage_levels = STAGE_LEVELS,
  not_staged   = NOT_STAGED,
  built        = as.character(Sys.time()))

## ---- join onto app$meta ----------------------------------------------------
# app$meta is a stratified downsample (~30k of 58,917), so this is a SUBSET join by design.
# Aggregate numbers in the app must come from app$tricycle$percell, never from meta, or they
# will not match the analysis they claim to show.
i <- match(app$meta$cell, pc$cell)
app$meta$tricycle_theta   <- pc$tricycle_theta[i]
app$meta$tricycle_stage   <- pc$tricycle_stage[i]
app$meta$tricycle_cycling <- pc$tricycle_cycling[i]
matched <- sum(!is.na(i))
cat(sprintf("\nJoined onto app$meta: %s of %s cells matched (%.1f%%)\n",
            format(matched, big.mark = ","), format(nrow(app$meta), big.mark = ","),
            100 * matched / nrow(app$meta)))
if (matched < 0.9 * nrow(app$meta))
  warning("under 90% of app$meta matched a tricycle cell -- check the barcode convention")

if (!is.null(app$cm$meta)) {
  j <- match(app$cm$meta$cell, pc$cell)
  app$cm$meta$tricycle_theta   <- pc$tricycle_theta[j]
  app$cm$meta$tricycle_stage   <- pc$tricycle_stage[j]
  app$cm$meta$tricycle_cycling <- pc$tricycle_cycling[j]
  cat(sprintf("Joined onto app$cm$meta: %s of %s cells matched\n",
              format(sum(!is.na(j)), big.mark = ","), format(nrow(app$cm$meta), big.mark = ",")))
}

## ---- summary ---------------------------------------------------------------
cat("\n== stage composition (all cells) ==\n")
print(table(pc$tricycle_stage))
cat("\n== cycling %, cardiomyocytes ==\n")
cm <- by_ct[by_ct$celltype == "Cardiomyocyte", c("timepoint","genotype","n",
                                                 "pct_cycling_seurat","pct_cycling_tricycle")]
print(cm, row.names = FALSE)
if (!is.null(controls)) {
  cat("\n== headline controls ==\n")
  print(controls[controls$metric %in% c("kappa","pct_concordance","r_group_fractions",
                                        "n_ref_genes_matched","pct_abstained"), ], row.names = FALSE)
}
if (!is.null(matched_gt)) {
  cat("\n== KO-WT gap, raw vs at matched depth ==\n")
  print(matched_gt[, c("celltype","timepoint","gap_raw","gap_matched",
                       "auc_depth_raw","auc_depth_matched")], row.names = FALSE)
} else cat("\n(no depth-matched table -- run cellcycle_tricycle_depthmatched.R for it)\n")

cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat(sprintf("Done. Bundle is now %.1f MB.\n", file.size("app_data.rds") / 1024^2))
