#!/usr/bin/env Rscript
# Cell-level KO-vs-WT within a CM subcluster -- for the subclusters pseudobulk cannot test.
#
# WHY THIS EXISTS, AND WHAT IT IS NOT. The pseudobulk pipeline
# (cm_subcluster_analyze.R) aggregates cells into one sample per library x lane, needs 20
# cells per sample and two samples per genotype, and skips a subcluster that cannot supply
# them -- CM12 at res 0.2 is the case in point: 99 cells spread 10/10/18/21/8/8/10/14, so
# exactly one sample clears the floor. That leaves the subcluster with no KO-vs-WT output
# at all, which reads as "nothing to see" when the honest statement is "not testable the
# way the rest of the pipeline tests things".
#
# This runs the comparison the other way: cell by cell, Wilcoxon, KO cells against WT
# cells. It is a RANKING, not a test.
#
#   The p-values are overconfident and should not be quoted. All the KO cells in a
#   subcluster come from two animals (one P0 KO, one P7 KO), so 36 cells are 36
#   measurements of two mice, not 36 independent observations. Wilcoxon does not know
#   that and will report significance proportional to cell count. This is the
#   pseudoreplication the pseudobulk design exists to avoid, accepted deliberately here
#   because a ranked gene list from a 99-cell cluster is still a useful thing to look at
#   as long as nobody calls it evidence.
#
# Reported alongside: AUC (the effect size that does not inflate with n -- 0.5 is no
# separation, 1.0 is perfect), the per-group cell counts, and a flag for the sex/construct
# confounders the rest of the project excludes.
#
#   Rscript cm_subcluster_celllevel_de.R [--res=0.2] [--subcluster=CM12] [--all-skipped]
#
# Writes results/tables/cm_subcluster_celllevel_<res>_<CL>.csv

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(presto)))
options(future.globals.maxSize = 16 * 1024^3)

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
RES         <- argval("--res", "0.2")
WANT        <- argval("--subcluster", "CM12")
ALL_SKIPPED <- "--all-skipped" %in% commandArgs(TRUE)
MIN_CELLS   <- 10L      # per genotype; below this even a ranking is not worth printing
CONFOUND <- c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")
# CM12 is immune contamination (docs/05-cm-deepdive.qmd#cm12). Flag the markers that say
# so, because a KO-vs-WT difference inside it is a difference between IMMUNE cells and
# must not be read as cardiomyocyte biology.
IMMUNE <- c("Ptprc","Cd52","Laptm5","Coro1a","Cd3e","Trbc1","Cd79a",
            "Cpa3","Kit","Cma1","Mcpt4","Tpsb2","Tpsab1","Ms4a2","Fcer1a","Mrgprb2","Il1rl1")

cm <- readRDS(file.path(PROC, "seurat.cm.subclustered.rds"))
cm$genotype <- genotype_of(cm$orig.ident)
res_col <- paste0("SCT_snn_res.", RES)
stopifnot(res_col %in% colnames(cm@meta.data))
cm$cl <- paste0("CM", cm[[res_col]][, 1])

targets <- if (ALL_SKIPPED) {
  f <- file.path(OUTTAB, sprintf("cm_subcluster_res%s_KOvsWT_summary.csv", RES))
  if (!file.exists(f)) stop("no summary at ", f, " -- run cm_subcluster_analyze.R first")
  s <- read.csv(f, stringsAsFactors = FALSE)
  s$subcluster[s$status != "ok"]
} else WANT
cat(sprintf("res %s | subclusters to test cell-level: %s\n\n", RES, paste(targets, collapse = ", ")))

DefaultAssay(cm) <- "RNA"
cm <- NormalizeData(cm, verbose = FALSE)
X  <- SeuratObject::GetAssayData(cm, assay = "RNA", layer = "data")

run_one <- function(cl, cells, tag) {
  g <- factor(cm$genotype[match(cells, colnames(cm))], levels = c("WT", "KO"))
  n <- table(g)
  if (min(n) < MIN_CELLS) {
    cat(sprintf("  %-6s %-10s WT=%d KO=%d -- below %d per genotype, skipped\n",
                cl, tag, n[["WT"]], n[["KO"]], MIN_CELLS)); return(NULL)
  }
  M <- X[, cells, drop = FALSE]
  M <- M[Matrix::rowSums(M > 0) >= 3, , drop = FALSE]     # a gene seen in <3 cells ranks nothing
  w <- suppressWarnings(presto::wilcoxauc(M, as.character(g)))
  w <- w[w$group == "KO", ]                                # KO vs rest == KO vs WT with two groups
  w$confounder   <- w$feature %in% CONFOUND
  w$immune_gene  <- w$feature %in% IMMUNE
  w$n_KO_cells   <- as.integer(n[["KO"]]); w$n_WT_cells <- as.integer(n[["WT"]])
  w$stratum      <- tag
  w$NOTE <- "descriptive_cell_level_pseudoreplicated_pvalues_not_valid"
  w <- w[order(-abs(w$auc - 0.5)), ]                       # rank by effect size, not p
  names(w)[names(w) == "feature"] <- "gene"
  keep <- c("gene","stratum","auc","logFC","pct_in","pct_out","pval","padj",
            "n_KO_cells","n_WT_cells","confounder","immune_gene","NOTE")
  w <- w[, intersect(keep, names(w))]
  cat(sprintf("  %-6s %-10s WT=%3d KO=%3d | %5d genes tested | AUC>=0.7 or <=0.3: %d (non-confounder %d)\n",
              cl, tag, n[["WT"]], n[["KO"]], nrow(w),
              sum(abs(w$auc - 0.5) >= 0.2), sum(abs(w$auc - 0.5) >= 0.2 & !w$confounder)))
  w
}

for (cl in targets) {
  cells <- colnames(cm)[cm$cl == cl]
  if (!length(cells)) { cat(sprintf("  %-6s not present at res %s\n", cl, RES)); next }
  out <- list(run_one(cl, cells, "pooled"))
  # per timepoint too: pooled KO-vs-WT can be a P0-vs-P7 difference in disguise if the
  # genotypes sit unevenly across timepoints
  for (tp in TIMEPOINTS) {
    sub <- cells[cm$timepoint[match(cells, colnames(cm))] == tp]
    if (length(sub)) out <- c(out, list(run_one(cl, sub, tp)))
  }
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) { cat(sprintf("  %-6s nothing testable\n", cl)); next }
  f <- file.path(OUTTAB, sprintf("cm_subcluster_celllevel_res%s_%s.csv", RES, cl))
  write.csv(out, f, row.names = FALSE)
  cat(sprintf("  -> %s\n", basename(f)))
  top <- out[out$stratum == "pooled" & !out$confounder, ]
  cat("\n  top 12 by |AUC-0.5| (pooled, confounders excluded):\n")
  print(head(top[, c("gene","auc","logFC","pct_in","pct_out","immune_gene")], 12), row.names = FALSE)
  cat("\n")
}
cat("=== DONE cm_subcluster_celllevel_de ===\n")
