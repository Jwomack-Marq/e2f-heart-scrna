#!/usr/bin/env Rscript
# Is the KO-vs-WT cell-cycle difference a depth artifact, or is depth HIDING it?
#
# WHY. cellcycle_tricycle.R established two things that have to be read together:
#
#   1. At P7 the knockout has more cycling cardiomyocytes than the wild type
#      (23.3 % vs 20.0 % on the theta arc; 31.6 % vs 25.6 % by Seurat).
#   2. The cycling call is strongly depth-dependent. Within one cell type at one
#      timepoint the fraction climbs monotonically with detected genes.
#
# The obvious worry is that (2) manufactures (1). It does not -- it works the other way,
# and that is the point of this script. ko_ml_depth_matching.csv records that at P7 the
# WILD TYPE is the deeper library (median 18.0k UMIs vs the KO's 12.7k, depth AUC 0.33).
# Deeper sequencing inflates the cycling call, so the deeper genotype is the one being
# flattered -- and it is the one with the LOWER number. The raw comparison is therefore
# biased AGAINST the observed difference, and the honest estimate should be larger, not
# smaller. That is the opposite of the usual confound story, so it needs measuring rather
# than asserting.
#
# WHY NOT JUST STRATIFY. Binning cells by depth and comparing within bins (which is what
# cellcycle_tricycle.R does for P0 vs P7) conditions on the confound but does not remove
# it: the two genotypes barely overlap in depth at P7, so the common bins are the tails of
# both distributions and the weights land on whichever bins happen to be shared. Worse, the
# cycling CALL is itself computed from depth-sensitive counts, so a cell keeps its inflated
# label inside its bin. The only way to answer the question is to make the two libraries
# the same depth and then ASK TRICYCLE AGAIN.
#
# HOW. Quantile-floor binomial thinning, per cell type x timepoint stratum -- the same
# procedure ko_signature_ml.R uses, deliberately, so the two analyses cannot disagree about
# what "depth-matched" means. For cell i in genotype g at ECDF position u_i, the target
# depth is the smaller of the two genotypes' u_i-quantiles, so both marginal depth
# distributions become their pointwise minimum; cells at or below target keep everything.
# Then the thinned counts are re-normalised and pushed through project_cycle_space() and
# estimate_cycle_position() from scratch.
#
# BOTH projections run here, raw and matched, rather than reusing the earlier run's theta.
# project_cycle_space() centres each gene across the cells it is given, so a projection of
# the thinned matrix is not on the same footing as one of the raw matrix; computing both in
# one pass is what makes the comparison a comparison.
#
#   Reads  processing/seurat.combined.annotated.rds
#   Writes results/tables/cellcycle_tricycle_matched_percell.csv.gz
#          results/tables/cellcycle_tricycle_matched_genotype.csv
#
#   Rscript cellcycle_tricycle_depthmatched.R [--max-cells=N]
#
# Runs in e2f-tricycle:latest (our_analysis/Dockerfile.tricycle).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({
  library(tricycle); library(SingleCellExperiment); library(Matrix)
}))

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
MAXCELLS <- as.integer(argval("--max-cells", "0"))
SEED     <- 42L
MIN_CELLS <- 50L                      # per genotype, to match a stratum at all
CYC_LO <- 0.5 * pi; CYC_HI <- 1.5 * pi   # the arc, as established by cellcycle_tricycle.R
NOTE <- "descriptive_n1_sex_confounded"
say <- function(...) cat(sprintf("[tri_depth] %s\n", paste0(...)))

## ---- 1. cells -------------------------------------------------------------
obj <- file.path(PROC, "seurat.combined.annotated.rds")
stopifnot(file.exists(obj))
say("reading ", basename(obj), " ...")
comb <- readRDS(obj)
DefaultAssay(comb) <- "RNA"
counts <- SeuratObject::GetAssayData(comb, assay = "RNA", layer = "counts")

md <- comb@meta.data[colnames(counts), , drop = FALSE]   # one authoritative order
stopifnot(identical(rownames(md), colnames(counts)))
md$genotype <- genotype_of(md$orig.ident)
md$nCount   <- Matrix::colSums(counts)
md$stratum  <- paste(md$celltype, md$timepoint, sep = "|")

if (MAXCELLS > 0 && MAXCELLS < ncol(counts)) {
  set.seed(SEED); keep <- sort(sample(ncol(counts), MAXCELLS))
  counts <- counts[, keep]; md <- md[keep, , drop = FALSE]
  say(sprintf("SMOKE TEST: %s cells", format(ncol(counts), big.mark = ",")))
}
say(sprintf("%s cells | %d strata", format(ncol(counts), big.mark = ","),
            length(unique(md$stratum))))

## ---- 2. quantile-floor binomial thinning ----------------------------------
strat_n <- table(md$stratum, md$genotype)
ok_strata <- rownames(strat_n)[strat_n[, "KO"] >= MIN_CELLS & strat_n[, "WT"] >= MIN_CELLS]
say("matching ", length(ok_strata), " strata with >= ", MIN_CELLS, " cells per genotype")

p_keep <- rep(1, nrow(md))
for (s in ok_strata) {
  iK <- which(md$stratum == s & md$genotype == "KO")
  iW <- which(md$stratum == s & md$genotype == "WT")
  nK <- md$nCount[iK]; nW <- md$nCount[iW]
  for (side in list(list(i = iK, n = nK), list(i = iW, n = nW))) {
    u   <- (rank(side$n, ties.method = "first") - 0.5) / length(side$n)
    tgt <- pmin(quantile(nK, u, type = 7, names = FALSE),
                quantile(nW, u, type = 7, names = FALSE))
    p_keep[side$i] <- pmin(1, tgt / side$n)
  }
}
set.seed(SEED)
thin <- counts
pc <- rep(p_keep, diff(thin@p))
thin@x <- as.numeric(rbinom(length(thin@x), size = as.integer(round(thin@x)), prob = pc))
thin <- Matrix::drop0(thin); rm(pc)
md$nCount_matched <- Matrix::colSums(thin)
say(sprintf("UMIs retained: %.1f%% overall", 100 * sum(md$nCount_matched) / sum(md$nCount)))

## ---- 3. tricycle, twice ----------------------------------------------------
# Mann-Whitney AUC of depth against genotype: 0.5 means the two are indistinguishable on
# depth, which is what the matching is for and what has to be checked rather than assumed.
auc_rank <- function(score, pos) {
  pos <- as.logical(pos); n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score); (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

theta_of <- function(cnt, label) {
  say("projecting (", label, ") ...")
  lognorm <- Seurat::NormalizeData(Seurat::CreateSeuratObject(cnt), verbose = FALSE)
  X <- SeuratObject::GetAssayData(lognorm, assay = "RNA", layer = "data")
  sce <- SingleCellExperiment(assays = list(logcounts = X))
  msg <- character(0)
  sce <- withCallingHandlers(
    project_cycle_space(sce, gname.type = "SYMBOL", species = "mouse"),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
  ng <- suppressWarnings(as.integer(sub(".*is (\\d+).*", "\\1",
          grep("number of projection genes", msg, value = TRUE)[1])))
  say(sprintf("  %s reference genes matched", if (is.na(ng)) "?" else format(ng)))
  sce <- estimate_cycle_position(sce)
  th <- sce$tricyclePosition
  stopifnot(length(th) == ncol(cnt))
  th
}

theta_raw     <- theta_of(counts, "raw")
theta_matched <- theta_of(thin,   "depth-matched")

cyc <- function(th) !is.na(th) & th >= CYC_LO & th <= CYC_HI
out <- data.frame(
  cell = colnames(counts), celltype = md$celltype, timepoint = md$timepoint,
  genotype = md$genotype, stratum = md$stratum,
  ngene = md$nFeature_RNA, numi = md$nCount, numi_matched = md$nCount_matched,
  theta_raw = round(theta_raw, 4), theta_matched = round(theta_matched, 4),
  cycling_raw = cyc(theta_raw), cycling_matched = cyc(theta_matched),
  matched = md$stratum %in% ok_strata,
  stringsAsFactors = FALSE)
# gzipped, like pcdims_cm_percell.csv.gz: 58,917 rows of mostly-numeric columns is ~7 MB
# flat and about a fifth of that compressed, and this file is tracked.
.pc <- gzfile(file.path(OUTTAB, "cellcycle_tricycle_matched_percell.csv.gz"), "wt")
write.csv(out, .pc, row.names = FALSE); close(.pc)

## ---- 4. the answer ---------------------------------------------------------
gt <- do.call(rbind, lapply(ok_strata, function(s) {
  z <- out[out$stratum == s, , drop = FALSE]
  iK <- z$genotype == "KO"; iW <- z$genotype == "WT"
  r_ko <- 100 * mean(z$cycling_raw[iK]);     r_wt <- 100 * mean(z$cycling_raw[iW])
  m_ko <- 100 * mean(z$cycling_matched[iK]); m_wt <- 100 * mean(z$cycling_matched[iW])
  data.frame(
    celltype = z$celltype[1], timepoint = z$timepoint[1],
    n_KO = sum(iK), n_WT = sum(iW),
    median_numi_KO_raw = median(z$numi[iK]), median_numi_WT_raw = median(z$numi[iW]),
    median_numi_KO_matched = median(z$numi_matched[iK]),
    median_numi_WT_matched = median(z$numi_matched[iW]),
    auc_depth_raw     = round(auc_rank(log(z$numi), iK), 3),
    auc_depth_matched = round(auc_rank(log(z$numi_matched), iK), 3),
    pct_cycling_KO_raw = round(r_ko, 1), pct_cycling_WT_raw = round(r_wt, 1),
    gap_raw = round(r_ko - r_wt, 1),
    pct_cycling_KO_matched = round(m_ko, 1), pct_cycling_WT_matched = round(m_wt, 1),
    gap_matched = round(m_ko - m_wt, 1),
    NOTE = NOTE, stringsAsFactors = FALSE)
}))
gt <- gt[order(gt$celltype, gt$timepoint), ]; rownames(gt) <- NULL
write.csv(gt, file.path(OUTTAB, "cellcycle_tricycle_matched_genotype.csv"), row.names = FALSE)

cat("\n=== did the matching actually remove depth? (AUC 0.5 = indistinguishable) ===\n")
print(gt[, c("celltype","timepoint","n_KO","n_WT","median_numi_KO_raw","median_numi_WT_raw",
             "auc_depth_raw","auc_depth_matched")], row.names = FALSE)

cat("\n=== KO - WT cycling fraction, raw vs at matched depth ===\n")
print(gt[, c("celltype","timepoint","pct_cycling_KO_raw","pct_cycling_WT_raw","gap_raw",
             "pct_cycling_KO_matched","pct_cycling_WT_matched","gap_matched")], row.names = FALSE)

cm <- gt[gt$celltype == "Cardiomyocyte", ]
for (i in seq_len(nrow(cm))) {
  r <- cm[i, ]
  dir <- if (is.na(r$gap_matched) || is.na(r$gap_raw)) "?" else
    if (abs(r$gap_matched) > abs(r$gap_raw)) "WIDENS -- depth was masking it" else
    if (sign(r$gap_matched) != sign(r$gap_raw)) "REVERSES" else "narrows"
  cat(sprintf("  cardiomyocytes %s: KO-WT %+.1f raw -> %+.1f matched  (%s)\n",
              r$timepoint, r$gap_raw, r$gap_matched, dir))
}
cat("\n  n = 1 animal per genotype x timepoint, and genotype is confounded with sex.\n")
cat("  Depth matching removes DEPTH. It cannot remove sex, animal or library.\n")

say("wrote tables")
cat("=== DONE cellcycle_tricycle_depthmatched ===\n")
