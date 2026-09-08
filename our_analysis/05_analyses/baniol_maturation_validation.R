#!/usr/bin/env Rscript
# EXTERNAL VALIDATION of the CM maturation model against Baniol et al. 2021
# (Exp Cell Res 408:112880; ENA PRJEB47622) -- 285 FACS-sorted cardiomyocytes,
# P0 and P7, Smart-seq2. This is the ONE external test this dataset supports.
#
# WHAT THIS CAN AND CANNOT TEST.
#   CAN: cm_maturation_glmnet.rds, which predicts TIMEPOINT (P0 vs P7). Baniol's
#        `stage` is a genuine external label -- different lab, protocol and mice --
#        so this is a real out-of-sample test, not a restatement of our clustering.
#   CANNOT: cell_type_glmnet.rds. Baniol is 100% cardiomyocytes (vCM/aCM), so a
#        7-class classifier has no negative class there. "285/285 called
#        Cardiomyocyte" would be a tautology, not evidence. It is run below purely
#        as a one-sided smoke test: it can falsify, it cannot confirm.
#
# WHICH STORE, AND WHY IT MATTERS. Two Baniol stores exist with identical genes,
# cells, indptr and indices but DIFFERENT data.bin -- they are two normalisations
# of one count matrix, not copies (they are not even monotone transforms of each
# other). Measured per-cell:
#     Baniol2021_FUCCI    (ingest.py)         sum(expm1(x)) = 10,000  <- natural-log CP10K
#     Baniol2021_FUCCI_R  (sce_final_umap.rds) sum(expm1(x)) = 1.4e7  <- some other scale
# Seurat's NormalizeData default is natural-log CP10K, so the ingest.py store is on
# OUR scale and the _R store is not. We use the former. BANIOL_DATASET can override.
#
# PRE-REGISTERED CAVEAT -- the sort is asymmetric and it cuts AGAINST us.
# model/cmcycle/baniol.py records that P7 is 4.5-5.2x enriched for cycling cells
# while P0 is essentially unenriched. Cycling is an IMMATURE phenotype, so Baniol's
# P7 population is deliberately skewed to look less mature than P7 really is.
# Therefore, decided before seeing the result:
#   high AUC  -> strong evidence: the model separated P0 from P7 despite a sort
#                that pushes P7 toward the P0 phenotype.
#   low AUC   -> AMBIGUOUS, not a refutation: could be the model, the platform gap,
#                or the sort. The P(P7)-vs-cycling_score correlation below is what
#                distinguishes them -- if P(P7) tracks cycling rather than stage,
#                the sort is the explanation.
#
# THREE ARMS, because the shipped applier has a bug worth quantifying.
#   naive     : the shipped model + the shipped zero-fill applier. Absent features
#               are set to 0, but in log-normalised space 0 means "not expressed",
#               NOT "unmeasured" -- so ~27% of the coefficient genes actively PUSH
#               the prediction instead of abstaining. This arm measures that damage.
#   shared_raw: refit on our CMs using only genes present in BOTH datasets. No
#               zero-fill, so no phantom absences.
#   shared_z  : as shared_raw, but each gene z-scored WITHIN each dataset first, so
#               per-gene scale differences between Smart-seq2 and 3' droplet data
#               cannot masquerade as biology. This is the primary arm.
# Report AUC, not accuracy: the two platforms are not on a common absolute scale,
# so a 0.5 probability threshold is not meaningful, but ranking is.
#
# Reads  processing/seurat.combined.annotated.rds, results/models/*.rds,
#        and the Baniol CSR store (BANIOL_ROOT, mounted read-only under Docker).
# Writes results/tables/baniol_*.csv, results/figures/baniol_maturation_validation.png
#
# Runs in e2f-ml:latest (our_analysis/Dockerfile.ml -- the only image with glmnet).
#   docker run --rm -u $(id -u):$(id -g) \
#     -v "$PWD/our_analysis":/work \
#     -v /home/justin/Projects/lab-server/apps/cardiac-rnaseq-db:/baniol:ro \
#     -e BANIOL_ROOT=/baniol -w /work e2f-ml:latest \
#     Rscript 05_analyses/baniol_maturation_validation.R [--smoke]

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(glmnet); library(Matrix); library(ggplot2) }))
set.seed(1)

SMOKE   <- any(grepl("^--smoke$", commandArgs(TRUE)))
BROOT   <- Sys.getenv("BANIOL_ROOT", "/home/justin/Projects/lab-server/apps/cardiac-rnaseq-db")
BDSET   <- Sys.getenv("BANIOL_DATASET", "Baniol2021_FUCCI")
MODELS  <- file.path(RESULTS, "models")
say <- function(...) { cat(sprintf(...), "\n"); flush(stdout()) }
say("=== baniol_maturation_validation%s ===", if (SMOKE) " [SMOKE]" else "")

## ---- 1. read the Baniol gene-major CSR store -------------------------------
# Layout (little-endian, headerless): indptr.bin int32[n_genes+1], indices.bin
# int32 cell indices, data.bin float32 log-normalised values. Assert the shapes --
# these are raw dumps, so a truncated copy would otherwise read as plausible.
read_baniol <- function(root, dset) {
  e <- file.path(root, "expr", dset)
  d <- file.path(root, "data", dset)
  stopifnot(dir.exists(e), dir.exists(d))
  genes <- readLines(file.path(e, "genes.txt"))
  cells <- readLines(file.path(e, "cells.txt"))
  ng <- length(genes); nc <- length(cells)
  indptr <- readBin(file.path(e, "indptr.bin"), "integer", n = ng + 1L, size = 4L)
  stopifnot(length(indptr) == ng + 1L, indptr[1] == 0L)
  nnz <- indptr[ng + 1L]
  idx <- readBin(file.path(e, "indices.bin"), "integer", n = nnz, size = 4L)
  dat <- readBin(file.path(e, "data.bin"),    "numeric", n = nnz, size = 4L)
  stopifnot(length(idx) == nnz, length(dat) == nnz,
            min(idx) >= 0L, max(idx) < nc)
  m <- sparseMatrix(i = rep.int(seq_len(ng), diff(indptr)), j = idx + 1L, x = dat,
                    dims = c(ng, nc), dimnames = list(genes, cells))
  md <- read.csv(gzfile(file.path(d, "meta.csv.gz")), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(md$ID), setequal(md$ID, cells))
  md <- md[match(cells, md$ID), ]
  # Confirm the normalisation is the one we think it is before trusting any transfer.
  cp <- median(Matrix::colSums(expm1(m)))
  say("Baniol store '%s': %d genes x %d cells, nnz=%d, median sum(expm1)=%.0f",
      dset, ng, nc, nnz, cp)
  if (abs(cp - 1e4) > 1) stop(sprintf(
    "store '%s' is not natural-log CP10K (median sum(expm1)=%.0f, expected 10000); ",
    dset, cp), "it is not on Seurat's NormalizeData scale -- refusing to transfer.")
  list(m = m, meta = md)
}
ban  <- read_baniol(BROOT, BDSET)
bmat <- ban$m; bmeta <- ban$meta
bstage <- factor(bmeta$stage, levels = c("P0", "P7"))
stopifnot(!anyNA(bstage))
say("Baniol stage: P0=%d P7=%d | CellType: %s", sum(bstage == "P0"), sum(bstage == "P7"),
    paste(names(table(bmeta$CellType)), table(bmeta$CellType), sep = "=", collapse = " "))

## ---- 2. our cardiomyocytes, same subset rule as cell_state_classifier.R ----
cm_bundle <- readRDS(file.path(MODELS, "cm_maturation_glmnet.rds"))
ct_bundle <- readRDS(file.path(MODELS, "cell_type_glmnet.rds"))
feats <- cm_bundle$features                # the shipped model's own 2000 HVGs
stopifnot(cm_bundle$family == "binomial", cm_bundle$positive == "P7")

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
comb <- comb[, comb$lane == "lane1" & comb$celltype == "Cardiomyocyte"]
DefaultAssay(comb) <- "RNA"
comb <- NormalizeData(comb, verbose = FALSE)
ycm <- factor(comb$timepoint, levels = c("P0", "P7"))
if (SMOKE) {                               # keep the code path, shrink the data
  k <- unlist(lapply(split(seq_along(ycm), ycm), function(i) sample(i, min(600, length(i)))))
  comb <- comb[, k]; ycm <- droplevels(ycm[k])
}
ourlog <- GetAssayData(comb, assay = "RNA", layer = "data")
say("our CMs: %d (P0=%d, P7=%d)", ncol(comb), sum(ycm == "P0"), sum(ycm == "P7"))

shared <- intersect(feats, rownames(bmat))
say("model features: %d | present in Baniol: %d (%.1f%%) | zero-filled by the shipped applier: %d",
    length(feats), length(shared), 100 * length(shared) / length(feats),
    length(feats) - length(shared))

## ---- 3. helpers ------------------------------------------------------------
# AUC as the Mann-Whitney U statistic; ties get mid-ranks, so this is exact.
auc <- function(score, pos) {
  if (!any(pos) || all(pos)) return(NA_real_)
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
zscore <- function(M) {                    # per gene (row), within one dataset
  mu <- Matrix::rowMeans(M)
  sdv <- sqrt(pmax(Matrix::rowMeans(M^2) - mu^2, 0))
  sdv[sdv < 1e-8] <- 1                     # constant gene -> contributes 0, not NaN
  (as.matrix(M) - mu) / sdv
}
fit_predict <- function(Xtr, ytr, Xte, nfolds = 5) {
  f <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0.9,
                 type.measure = "class", nfolds = nfolds)
  list(fit = f,
       p = as.vector(predict(f, Xte, s = "lambda.1se", type = "response")),
       nz = sum(as.matrix(coef(f, s = "lambda.1se"))[-1, 1] != 0))
}
# Internal held-out AUC on OUR data: the ceiling this arm could reach externally.
internal_auc <- function(X, y) {
  te_i <- unlist(lapply(split(seq_along(y), y), function(i) sample(i, ceiling(0.2 * length(i)))))
  te <- seq_along(y) %in% te_i
  p <- fit_predict(X[!te, , drop = FALSE], y[!te], X[te, , drop = FALSE])$p
  auc(p, y[te] == "P7")
}

## ---- 3b. how much of the model's ACTUAL signal is missing -------------------
# 867/2000 features are absent from Baniol, but most carry a zero coefficient and so
# cost nothing. The number that matters is how many NON-ZERO coefficient genes get
# zero-filled, because only those actively push the naive prediction.
cm_coef <- as.matrix(coef(cm_bundle$model, s = "lambda.1se"))
coef_genes   <- setdiff(rownames(cm_coef)[cm_coef[, 1] != 0], "(Intercept)")
coef_missing <- setdiff(coef_genes, rownames(bmat))
say("shipped model: %d non-zero coefficients | %d absent from Baniol (%.0f%% of the model's real signal is zero-filled by the naive arm)",
    length(coef_genes), length(coef_missing),
    100 * length(coef_missing) / length(coef_genes))

## ---- 4. arm 1: the shipped model + the shipped zero-fill applier -----------
# Replicates predict_cell_state() exactly, including the bug being measured.
Xn <- matrix(0, nrow = ncol(bmat), ncol = length(feats),
             dimnames = list(colnames(bmat), feats))
Xn[, shared] <- t(as.matrix(bmat[shared, , drop = FALSE]))
p_naive <- as.vector(predict(cm_bundle$model, Xn, s = "lambda.1se", type = "response"))

## ---- 5. arms 2 and 3: refit on the shared gene space -----------------------
Xour_s <- Matrix::t(ourlog[shared, , drop = FALSE])
Xban_s <- Matrix::t(bmat[shared, , drop = FALSE])
r_raw  <- fit_predict(Xour_s, ycm, Xban_s)

Xour_z <- t(zscore(ourlog[shared, , drop = FALSE]))
Xban_z <- t(zscore(bmat[shared, , drop = FALSE]))
r_z    <- fit_predict(Xour_z, ycm, Xban_z)

## ---- 6. score every arm ----------------------------------------------------
pos <- bstage == "P7"
arms <- list(naive = p_naive, shared_raw = r_raw$p, shared_z = r_z$p)
int_auc <- c(naive = NA_real_,                       # shipped model: no refit to split
             shared_raw = internal_auc(Xour_s, ycm),
             shared_z   = internal_auc(Xour_z, ycm))
# COLUMN NAMES IN THIS STORE ARE MISLEADING, verified against the file rather than
# assumed: `cycling_score` holds the CATEGORICAL call ("cycling"/"noncycling") and
# `cc_score` holds the numeric score. Reading them by name alone yields all-NA.
stopifnot(all(bmeta$cycling_score %in% c("cycling", "noncycling")))
cyc     <- suppressWarnings(as.numeric(bmeta$cc_score))
cycling <- bmeta$cycling_score == "cycling"
stopifnot(!all(is.na(cyc)))
say("Baniol cycling: %d cycling / %d noncycling | cc_score range %.2f-%.2f",
    sum(cycling), sum(!cycling), min(cyc, na.rm = TRUE), max(cyc, na.rm = TRUE))
# Never let a degenerate confound column abort the run that produced the result.
safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

res <- do.call(rbind, lapply(names(arms), function(a) {
  p <- arms[[a]]
  data.frame(
    arm            = a,
    n_genes            = if (a == "naive") length(feats) else length(shared),
    n_feat_zero_filled = if (a == "naive") length(feats) - length(shared) else 0L,
    n_coef_zero_filled = if (a == "naive") length(coef_missing) else 0L,
    auc_baniol     = round(auc(p, pos), 3),
    auc_vCM        = round(auc(p[bmeta$CellType == "vCM"], pos[bmeta$CellType == "vCM"]), 3),
    auc_aCM        = round(auc(p[bmeta$CellType == "aCM"], pos[bmeta$CellType == "aCM"]), 3),
    acc_at_0.5     = round(mean((p > 0.5) == pos), 3),
    median_P7_in_P0 = round(median(p[!pos]), 3),
    median_P7_in_P7 = round(median(p[pos]), 3),
    internal_heldout_auc = round(int_auc[[a]], 3),
    # The confound check the pre-registered caveat calls for: if P(P7) tracks the
    # cycling score rather than stage, the asymmetric FACS sort is the explanation.
    spearman_vs_cc_score = round(safe_cor(p, cyc), 3),
    auc_vs_cycling_call  = round(auc(p, cycling), 3),
    stringsAsFactors = FALSE)
}))
res$NOTE <- paste("Judge by auc_baniol, NOT acc_at_0.5: the two platforms are not on a",
                  "common absolute scale, so the 0.5 threshold is uncalibrated across",
                  "them while ranking is unaffected. See baniol_stratified_auc.csv.")
print(res[, setdiff(names(res), "NOTE")], row.names = FALSE)
write.csv(res, file.path(OUTTAB, "baniol_maturation_validation.csv"), row.names = FALSE)

## ---- 6b. stratified RANK AUCs -- the confound argument rests on these -------
# acc_at_0.5 above is a CALIBRATION statement, and cross-platform calibration is
# expected to be off (Smart-seq2 vs 3' droplet are not on a common absolute scale),
# so a cell can rank correctly yet sit on the wrong side of 0.5. These are RANK
# statements and are immune to that. The decisive row is within_P0_cyc_vs_noncyc:
# if the model were really a cycling detector wearing a maturation label, it would
# be high. The P0 group is the clean place to ask, because P0 is the arm Baniol did
# NOT enrich, so cycling and stage are not confounded there.
strat <- do.call(rbind, lapply(names(arms), function(a) {
  p <- arms[[a]]; s7 <- bstage == "P7"; s0 <- bstage == "P0"; cy <- cycling
  f <- function(pos, neg)
    round(auc(c(p[pos], p[neg]), c(rep(TRUE, sum(pos)), rep(FALSE, sum(neg)))), 3)
  data.frame(arm = a,
             P7_vs_P0               = f(s7, s0),
             P7noncyc_vs_all_P0     = f(s7 & !cy, s0),
             P7cyc_vs_all_P0        = f(s7 &  cy, s0),
             within_P7_cyc_vs_nonc  = f(s7 &  cy, s7 & !cy),
             within_P0_cyc_vs_nonc  = f(s0 &  cy, s0 & !cy),
             stringsAsFactors = FALSE)
}))
strat$NOTE <- paste("RANK AUCs (calibration-free). within_P0_cyc_vs_nonc ~0.5 means the",
                    "model is NOT a cycling detector: in the group Baniol did not sort-enrich,",
                    "cycling carries no signal for it.")
print(strat, row.names = FALSE)
write.csv(strat, file.path(OUTTAB, "baniol_stratified_auc.csv"), row.names = FALSE)

percell <- data.frame(cell = colnames(bmat), stage = bmeta$stage, celltype = bmeta$CellType,
                      phase = bmeta$phase, cycling = bmeta$cycling_score,
                      p7_naive = round(p_naive, 4), p7_shared_raw = round(r_raw$p, 4),
                      p7_shared_z = round(r_z$p, 4), stringsAsFactors = FALSE)
write.csv(percell, file.path(OUTTAB, "baniol_percell_predictions.csv"), row.names = FALSE)
say("refit sparsity: shared_raw nonzero=%d, shared_z nonzero=%d", r_raw$nz, r_z$nz)

## ---- 7. cell-type classifier: one-sided smoke test, NOT validation ---------
# Baniol is pure cardiomyocyte, so a correct answer is uninformative and only a
# wrong one carries signal. Recorded with that caveat attached to the file.
ct_feats <- ct_bundle$features
Xc <- matrix(0, nrow = ncol(bmat), ncol = length(ct_feats),
             dimnames = list(colnames(bmat), ct_feats))
ct_shared <- intersect(ct_feats, rownames(bmat))
Xc[, ct_shared] <- t(as.matrix(bmat[ct_shared, , drop = FALSE]))
ct_pred <- as.vector(predict(ct_bundle$model, Xc, s = "lambda.1se", type = "class"))
ct_tab <- as.data.frame(table(predicted = ct_pred), stringsAsFactors = FALSE)
ct_tab$frac <- round(ct_tab$Freq / sum(ct_tab$Freq), 3)
ct_tab$NOTE <- paste("ONE-SIDED SMOKE TEST, NOT VALIDATION: Baniol is 100% sorted",
                     "cardiomyocytes, so there is no negative class. A CM call is a",
                     "tautology; only a non-CM call is informative.")
say("cell-type smoke test (%d/%d features present): %s", length(ct_shared), length(ct_feats),
    paste(ct_tab$predicted, ct_tab$Freq, sep = "=", collapse = "  "))
write.csv(ct_tab, file.path(OUTTAB, "baniol_celltype_smoketest.csv"), row.names = FALSE)

## ---- 8. figure -------------------------------------------------------------
# Encode cycling, not chamber: the visible bimodality within P7 IS the cycling split,
# and showing it is the whole confound argument. Chamber separates nothing here
# (auc_vCM/auc_aCM in the table above are both high), so it would only add noise.
pl <- do.call(rbind, lapply(names(arms), function(a)
  data.frame(arm = a, p7 = arms[[a]], stage = bstage,
             cycling = ifelse(cycling, "cycling", "noncycling"))))
pl$arm <- factor(pl$arm, levels = c("naive", "shared_raw", "shared_z"))
lab <- setNames(sprintf("%s (AUC = %.2f)", res$arm, res$auc_baniol), res$arm)
ggsave(file.path(OUTFIG, "baniol_maturation_validation.png"),
  ggplot(pl, aes(stage, p7, fill = stage)) +
    geom_violin(alpha = .5, scale = "width") +
    geom_jitter(aes(shape = cycling), width = .15, size = 1, alpha = .6) +
    scale_shape_manual(values = c(cycling = 17, noncycling = 1)) +
    facet_wrap(~ arm, labeller = labeller(arm = lab)) +
    scale_fill_manual(values = c(P0 = "#1565c0", P7 = "#c62828")) +
    geom_hline(yintercept = .5, linetype = 2, colour = "grey40") +
    theme_minimal(base_size = 12) +
    labs(title = "CM maturation model on Baniol 2021 (external, 285 sorted CMs)",
         subtitle = paste("Baniol P7 is cycling-enriched by the sort; within P0 (unenriched)",
                          "cycling separates at AUC 0.50 -- not a cycling detector"),
         x = "Baniol stage (true)", y = "predicted P(P7)"),
  width = 10, height = 4.5, dpi = 120)

say("=== DONE baniol_maturation_validation%s ===", if (SMOKE) " [SMOKE]" else "")
