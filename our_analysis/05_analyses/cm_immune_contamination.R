#!/usr/bin/env Rscript
# Flag immune cells sitting inside the cardiomyocyte compartment.
#
# The CM compartment contains ~99 leukocytes (about half of them mast cells) that
# annotate.R labelled "Cardiomyocyte". They are not doublets -- Scrublet clears all of
# them, and they carry half the RNA and a third the mitochondrial content of a real
# cardiomyocyte. They were mislabelled because annotate.R takes max.col over nine marker
# panels per whole-heart cluster and NONE of the nine describes a mast cell or a
# lymphocyte: the only immune panel is myeloid (Cd68/Lyz2/C1qa/Csf1r), so it scores 0.236
# while ambient CM RNA -- unavoidable in a 70 %-cardiomyocyte tissue -- scores 1.198. With
# no minimum score and no "Unassigned" option, Cardiomyocyte wins as the best of nine
# wrong answers. See docs/05-cm-deepdive.qmd#cm12.
#
# THIS SCRIPT DOES NOT RE-ANNOTATE ANYTHING. It writes a per-cell flag so the app can keep
# these cells off statistics they would distort. The real fix is in annotate.R (a mast /
# lymphocyte panel plus a minimum-score floor, and installing celldex so the SingleR
# cross-check actually runs) and belongs with the PIPseeker regeneration that the lane
# double-counting already requires -- doing it now means doing it twice.
#
# The flag is by MARKER SCORE, not by cluster id, on purpose: the same cells appear as
# CM12 in the shipped labelling, CM6/CM10/CM15 at dims 30 and CM6/CM11/CM15 at dims 50, and
# are not resolved at all at dims 10. A label-based flag would silently miss most of those.
#
#   Reads  processing/cm_variants.rds  (all 42,416 CM cells, full gene set)
#   Writes results/tables/cm_immune_contamination.csv
#
# Runs in e2f-seurat-full:latest.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
options(future.globals.maxSize = 16 * 1024^3)

# Pan-leukocyte, deliberately broad: the contamination is not only mast cells, so a
# mast-specific panel would miss the T cells in it. Ptprc alone is not enough -- ambient
# RNA puts a little of everything everywhere -- hence a >= 2 genes rule below.
IMMUNE <- c("Ptprc", "Cd52", "Laptm5", "Coro1a", "Cd3e", "Trbc1", "Cd79a")
MAST   <- c("Cpa3", "Kit", "Cma1", "Mcpt4", "Mcpt8", "Tpsb2", "Tpsab1", "Ms4a2",
            "Fcer1a", "Mrgprb2", "Il1rl1")
CM     <- c("Tnnt2", "Myh6", "Actc1", "Ttn", "Tnni3", "Myl2", "Myl3", "Actn2", "Des")
MIN_IMMUNE <- 2L        # genes detected; 1 would catch ambient Ptprc across the compartment

obj <- file.path(PROC, "cm_variants.rds")
if (!file.exists(obj)) stop("no CM object at ", obj,
                            "\n  run pcdims_sweep.R --object=cm first")
cat(sprintf("reading %s\n", basename(obj)))
cm <- readRDS(obj)
DefaultAssay(cm) <- "RNA"
cm <- NormalizeData(cm, verbose = FALSE)
X  <- SeuratObject::GetAssayData(cm, assay = "RNA", layer = "data")

have <- function(g) intersect(g, rownames(X))
ndet <- function(g) if (length(have(g))) Matrix::colSums(X[have(g), , drop = FALSE] > 0) else rep(0L, ncol(X))
scr  <- function(g) if (length(have(g))) Matrix::colMeans(X[have(g), , drop = FALSE]) else rep(0, ncol(X))

n_imm <- ndet(IMMUNE); n_mast <- ndet(MAST)
out <- data.frame(
  cell            = colnames(cm),
  n_immune_genes  = as.integer(n_imm),
  n_mast_genes    = as.integer(n_mast),
  immune_score    = round(scr(IMMUNE), 4),
  mast_score      = round(scr(MAST), 4),
  cm_score        = round(scr(CM), 4),
  nCount_RNA      = cm$nCount_RNA,
  percent_mt      = if ("percent.mt" %in% colnames(cm@meta.data)) cm$percent.mt else NA_real_,
  immune_contam   = n_imm >= MIN_IMMUNE,
  stringsAsFactors = FALSE)

f <- out$immune_contam
cat(sprintf("\nflagged %d of %d CM cells (%.2f%%) with >= %d pan-leukocyte genes detected\n",
            sum(f), nrow(out), 100 * mean(f), MIN_IMMUNE))
cat(sprintf("  of those, %.0f%% also detect >= 2 mast genes\n", 100 * mean(out$n_mast_genes[f] >= 2)))
cat(sprintf("  flagged   : CM score %.3f | %s counts | %.1f%% mito\n",
            mean(out$cm_score[f]), format(round(median(out$nCount_RNA[f])), big.mark = ","),
            median(out$percent_mt[f], na.rm = TRUE)))
cat(sprintf("  unflagged : CM score %.3f | %s counts | %.1f%% mito\n",
            mean(out$cm_score[!f]), format(round(median(out$nCount_RNA[!f])), big.mark = ","),
            median(out$percent_mt[!f], na.rm = TRUE)))

# Concordance with the clusters that isolate this population, as a check that the
# score-based flag finds what the label-based view found -- and, at dims 10, more.
for (vc in grep("^dims\\d+_res", colnames(cm@meta.data), value = TRUE)) {
  tab <- table(cm@meta.data[[vc]], f)
  if (ncol(tab) < 2) next
  frac <- tab[, "TRUE"] / rowSums(tab)
  hot  <- names(frac)[frac > 0.5]
  if (length(hot))
    cat(sprintf("  %-16s clusters >50%% flagged: %s (%s cells, %.0f%% of all flagged)\n",
                vc, paste0("CM", hot, collapse = ","),
                paste(tab[hot, "TRUE"], collapse = ","),
                100 * sum(tab[hot, "TRUE"]) / sum(f)))
  else
    cat(sprintf("  %-16s no cluster is >50%% flagged (they are dispersed at this cut)\n", vc))
}

write.csv(out[, setdiff(names(out), c("nCount_RNA", "percent_mt"))],
          file.path(OUTTAB, "cm_immune_contamination.csv"), row.names = FALSE)
cat(sprintf("\nwrote: %s\n", file.path(OUTTAB, "cm_immune_contamination.csv")))
cat("=== DONE cm_immune_contamination ===\n")
