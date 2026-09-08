#!/usr/bin/env Rscript
# KO-signature machine learning -- what elastic-net classifiers can honestly say about the
# E2f7/E2f8 knockout when there is ONE animal per condition.
#
# cell_state_classifier.R refuses to predict genotype, and the reason stands: KO is one
# male animal and WT one female animal per timepoint, the knockout is not visible in the
# transcript (E2f7/E2f8 are flat-to-higher in KO; see e2f_ko_verification.csv), and at P7
# the KO library has about HALF the UMI depth of the WT library. A KO-vs-WT classifier on
# these cells is therefore trivially near-perfect and, on its own, tells you nothing: it
# has learned sex, animal, library and depth. This script does not report that number as a
# result. It uses the classifier as an INSTRUMENT, under three controls, to ask questions
# whose answers are not fixed in advance by the design:
#
#   1. WHERE does the KO act?  Sex, animal and depth confounds are shared by every cell
#      type. If cardiomyocytes separate KO from WT far better than fibroblasts do, after
#      the controls, that difference is not explained by the shared confounds. Every
#      celltype x timepoint stratum is fitted at equal cell numbers against a
#      label-permutation null. Red blood cells at P0 are carried as an AMBIENT CONTROL:
#      if the RBC compartment separates as well as the cardiomyocytes, the signal is
#      animal-level, not cell-type biology.
#   2. Does the signature TRANSFER across animals?  P0 and P7 are different animal pairs.
#      A model trained KO-vs-WT at P0 and tested at P7 (and the reverse) is the only
#      genuine held-out-animal test this dataset allows. Transfer well above the
#      permuted-training null means the signature is consistent across two independent
#      pairs, not one animal's noise. Both pairs are male-KO / female-WT, so transfer can
#      NOT break the sex confound -- only the gene-level controls address that.
#   3. Is it E2F-TARGET DE-REPRESSION?  E2F7/8 repress the activating-E2F programme, so
#      the pre-specified E2F_TARGETS set (and MSigDB Hallmark E2F targets) is fitted as
#      the ONLY feature set and compared with hundreds of expression-matched random gene
#      sets of the same size. The percentile, and the sign of the univariate KO-WT
#      differences, test the mechanism instead of assuming it.
#   4. Does KO SHIFT MATURATION, and where in the CM compartment does the KO signal sit?
#      The shipped P0-vs-P7 stage model is applied to KO vs WT cardiomyocytes within each
#      timepoint as a phenotype readout (its sex/ambient weights neutralised), and the
#      out-of-fold P(KO) of a full CM model is summarised by cell-cycle phase and by the
#      production CM subclusters.
#   5. A confounder-screened CM KO GENE PANEL, cross-referenced against the pseudobulk
#      DESeq2 lists, next to a NAIVE panel fitted with no controls so the reader can see
#      what the controls removed.
#
# THE THREE CONTROLS
#   lane1 only     lane1/lane6 are the same library sequenced twice (97-100 % barcode
#                  overlap). Keeping both puts a cell's twin on both sides of every split.
#   sex genes      the project's 7-gene blocklist (Y genes, Xist, Tsix, the ROSA26
#                  construct) is removed, THEN a data-driven ubiquity filter drops genes
#                  whose KO-vs-WT difference has the same sign in every cell type at both
#                  timepoints and clears a magnitude threshold in most of them. A signal
#                  identical in every compartment is animal or sex, not knockout biology.
#                  One exception is built in: a gene whose difference is several-fold larger
#                  in cardiomyocytes than anywhere else (Myh7, Tcf4 in the first run) is NOT
#                  removed -- in a tissue that is 70-80 % cardiomyocyte, that pattern is a
#                  cardiomyocyte effect whose ambient-RNA shadow reaches every compartment,
#                  and the filter is not entitled to call it animal-level. Such genes are kept
#                  and flagged (cm_dominant) wherever they appear.
#                  The filter is run over the full universe as a positive control (it must
#                  rediscover Xist) and over a random WT/WT split as a negative control
#                  (it must remove ~nothing). Every removed gene is written out.
#   UMI depth      within each celltype x timepoint stratum, the deeper genotype's cells
#                  are binomially thinned so both genotypes share one depth distribution
#                  before normalisation. Thinning only, never upsampling.
#
# WHAT THIS IS NOT. n = 1 animal per condition, cells are pseudoreplicates, and genotype ==
# animal == sex == library. Nothing here is a p-value about the knockout; the permutation
# nulls test label exchangeability within this fixed pair of animals. Every table carries
# NOTE = descriptive_n1_sex_confounded_depth_matched. Depth matching removes UMI depth,
# not other library-level differences; the ubiquity filter is the second line, and a
# signal that is animal-specific AND cell-type-specific in magnitude is not separable
# from KO at n = 1. No fitted model is saved: these models are animal-specific by
# construction and must not be reused as classifiers.
#
# FINDINGS (first full run, 2026-09-04, seed 42, lane1 = 29,083 cells after dropping 73
# immune-contaminated CMs; recorded here so the write-up can be checked against the code):
#
#   Controls.  Depth matching took every stratum's depth-AUC to 0.50 (P7 cardiomyocytes were
#   at 0.33 raw: KO median 12.7k vs WT 18.0k UMIs). The ubiquity filter rediscovered Xist,
#   Ddx3y, Eif2s3y and Uty from the full universe, removed 0 genes on a random WT/WT split,
#   and removed 31 genes beyond the blocklist: rRNA-repeat artefacts (Gm42418, Gm40271,
#   Rps27rt, Rpl37rt), iron/stress/translation genes (Ftl1, Fth1, Cirbp, Ubb, Uba52, Eef1a1),
#   6 on chrX/Y. Myh7 (3.7x) and Tcf4 (5.0x) matched the ubiquitous pattern but were CM-
#   dominant and were kept.
#
#   1. WHERE.  Nowhere distinguishable from the animal effect. Every stratum separates KO
#      from WT at AUC 0.95-1.00 after the controls (permutation nulls 0.44-0.48), and the
#      RBC ambient control is the HIGHEST (0.999), above the cardiomyocytes (0.994/0.998).
#      Only P0 immune cells are lower (0.79, n = 100). A pervasive library/animal signal
#      survives the sex blocklist, the ubiquity filter and depth matching; the per-cell-type
#      comparison does not single out cardiomyocytes.
#   2. TRANSFER.  Cardiomyocytes are the exception that matters. A CM model trained at P0
#      scores P7 CMs at AUC 0.98 (and P7 -> P0 at 0.98; permuted-training null 0.50), while
#      endothelial, fibroblast, mural and immune signatures transfer at 0.42-0.57, i.e. not
#      at all. Whatever separates KO from WT cardiomyocytes is CONSISTENT across the two
#      independent animal pairs; whatever separates the other cell types is pair-specific.
#      Both pairs are male-KO / female-WT, so a CM-restricted sex effect is not excluded.
#      The consistent genes (same sign at both timepoints, 21 of 353): Tcf4 (+2.1 / +4.5,
#      the strongest weight in both models, CM-dominant), Adamts9, Myh7, Actg1, Apoe, Sdc4,
#      Tmod3 up in KO; Cox8b, Cox7a1, Acer2, Ebf1, Actn1, Fcrls down in KO. Tcf4 and Adamts9
#      are two of the seven non-confounder genes in shared_KO_up_P0_and_P7.txt.
#   3. E2F TARGETS.  Weakly and consistently UP in KO cardiomyocytes, but nearly useless as a
#      discriminator. 79 % (P0) and 97 % (P7) of the 29 E2F_TARGETS are higher in KO at
#      matched depth, mean difference +0.03 / +0.06 log units. Within a timepoint the set
#      separates KO from WT at AUC 0.57-0.59, BELOW expression-matched random 29-gene sets
#      (percentile 0.10 / 0.02) -- any 29 expressed genes carry more animal signal. In
#      transfer the set reaches only 0.54-0.57 but ranks at the 0.90-0.98 percentile of
#      random sets, i.e. what little it carries is the consistent part. Hallmark E2F (199
#      genes) behaves the same way. Read: the de-repression direction is there, the effect
#      size is small, and it is not what the KO classifier is learning.
#   4. MATURATION / STATE.  The shipped stage model scores KO cardiomyocytes slightly LESS
#      mature at both timepoints (Cliff's delta -0.17 at P0, -0.09 at P7, matched depth; raw
#      depth gives the same). P(KO) is uniformly separated in every phase and every
#      production subcluster except CM8 (the endothelial-transcript CM state; Cldn5/Esam/
#      Ly6a markers), where both genotypes drift toward the middle.
#   5. PANEL.  353 genes (225 P0, 159 P7, 31 in both); sign agreement with the pseudobulk
#      DESeq2 log2FC 0.80. The naive, uncontrolled panel puts Xist first and Tcf4 second;
#      the controlled panel puts Tcf4 first. Tcf4 is therefore the one candidate that
#      survives every control AND transfers across animal pairs -- and it is the gene whose
#      KO-WT difference is 5x larger in cardiomyocytes than in any other compartment.
#
#   Recorded against overclaiming: none of this separates knockout from the one male animal
#   per timepoint. What it does establish is that a cardiomyocyte-specific, cross-pair-
#   consistent signature exists, that it is not the E2F-target programme, and that the
#   E2F-target shift is real in direction but small in size.
#
#   Reads  processing/seurat.combined.annotated.rds
#          results/models/cm_maturation_glmnet.rds          (shipped stage model, applied)
#          results/tables/cm_immune_contamination.csv       (optional; flagged CMs dropped)
#          results/tables/pcdims_cm_percell.csv.gz          (production CM subcluster labels)
#          results/tables/percelltype_{P0,P7}_Cardiomyocyte_KOvsWT.descriptive.DE.csv
#          results/tables/shared_KO_up_P0_and_P7.txt
#   Writes results/tables/ko_ml_*.csv     results/figures/ko_ml_*.png
#
#   Rscript ko_signature_ml.R [--cores=8] [--smoke] [--no-naive] [--featureset-celltypes=all]
#     --smoke   tiny caps/repeats, a few minutes: exercises every code path, numbers not usable
#
# Runs in e2f-ml:latest (our_analysis/Dockerfile.ml).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({
  library(glmnet); library(Matrix); library(ggplot2); library(patchwork); library(parallel)
}))

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
NCORES   <- max(1L, as.integer(argval("--cores", "8")))
SMOKE    <- "--smoke" %in% commandArgs(TRUE)
NAIVE    <- !("--no-naive" %in% commandArgs(TRUE))
FS_CT    <- argval("--featureset-celltypes", "Cardiomyocyte")
SEED     <- 42L
say <- function(...) cat(sprintf("[ko_ml] %s\n", paste0(...)))
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
t0 <- Sys.time(); elapsed <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

NOTE     <- "descriptive_n1_sex_confounded_depth_matched"
CONFOUND <- c("Eif2s3y", "Kdm5d", "Uty", "Ddx3y", "Xist", "Tsix", "Gt(ROSA)26Sor")
ALPHA    <- 0.9
MIN_CELLS_ML     <- 100L   # per genotype, for any fitted stratum
MIN_CELLS_FILTER <- 50L    # per genotype, for the depth match and the ubiquity filter
CAP_SEP   <- 500L;  R_REP    <- 10L; N_PERM    <- 20L    # analysis 1
CAP_TRAIN <- 1000L; R_REP_TR <- 5L;  N_PERM_TR <- 20L    # analysis 2
N_RAND_E2F <- 200L; N_RAND_HALLMARK <- 100L; CAP_FS <- 1000L   # analysis 3
DELTA <- 0.10; FRAC_K <- 0.8; M_MIN <- 6L; DETECT_MIN <- 0.02  # ubiquity filter
CM_RATIO <- 3        # |d| in CM strata / median |d| elsewhere at or above this = CM-dominant, kept
MIN_STATE_N <- 30L   # cells per genotype for a subcluster/phase state to be summarised
N_HVG <- 2000L
if (SMOKE) {
  CAP_SEP <- 150L; CAP_TRAIN <- 150L; CAP_FS <- 150L
  R_REP <- 2L; N_PERM <- 3L; R_REP_TR <- 2L; N_PERM_TR <- 3L
  N_RAND_E2F <- 10L; N_RAND_HALLMARK <- 5L
  say("SMOKE MODE: caps/repeats reduced; outputs exercise the code paths only")
}
say("cores=", NCORES, "  seed=", SEED, "  alpha=", ALPHA, "  naive control=", NAIVE)
TP_COLS <- c(P0 = "#1565c0", P7 = "#c62828")
GT_COLS <- c(WT = "#546e7a", KO = "#ef6c00")
theme_set(theme_minimal(base_size = 12))
with_note <- function(df) { df$NOTE <- NOTE; df }
wtab <- function(df, name) {
  write.csv(with_note(df), file.path(OUTTAB, name), row.names = FALSE)
  say("wrote tables/", name, "  (", nrow(df), " rows)")
}

## ---- 1. load, lane-dedup, extract, free --------------------------------------------
say("loading combined annotated object ...")
comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
need <- c("lane", "celltype", "timepoint", "orig.ident", "Phase")
stopifnot(all(need %in% colnames(comb@meta.data)))
comb <- comb[, comb$lane == "lane1"]
stopifnot(length(unique(comb$lane)) == 1L)
DefaultAssay(comb) <- "RNA"
if (sum(grepl("^counts", SeuratObject::Layers(comb[["RNA"]]))) > 1L)
  comb[["RNA"]] <- SeuratObject::JoinLayers(comb[["RNA"]])
counts <- SeuratObject::GetAssayData(comb, assay = "RNA", layer = "counts")
stopifnot(inherits(counts, "dgCMatrix"))
md <- comb@meta.data[, c("orig.ident", "timepoint", "celltype", "Phase"), drop = FALSE]
md$genotype <- genotype_of(md$orig.ident)
if ("genotype" %in% colnames(comb@meta.data))
  stopifnot(identical(as.character(comb$genotype), md$genotype))
rm(comb); invisible(gc())
md$nCount <- Matrix::colSums(counts)
say("lane1 only: ", ncol(counts), " cells x ", nrow(counts), " genes  [", elapsed(), "]")

contam_f <- file.path(OUTTAB, "cm_immune_contamination.csv")
if (file.exists(contam_f)) {
  ct <- read.csv(contam_f, stringsAsFactors = FALSE)
  bad <- rownames(md) %in% ct$cell[ct$immune_contam %in% c(TRUE, "TRUE")] & md$celltype == "Cardiomyocyte"
  say("dropping ", sum(bad), " lane1 cardiomyocytes flagged immune-contaminated (cm_immune_contamination.csv)")
  counts <- counts[, !bad]; md <- md[!bad, , drop = FALSE]
} else say("cm_immune_contamination.csv not found; no CM contamination filter applied")
md$stratum <- paste(md$celltype, md$timepoint, sep = "|")
strat_n <- table(md$stratum, md$genotype)
print(strat_n)

## ---- 2. depth matching: quantile-floor binomial thinning per stratum --------------
# For cell i in genotype g with ECDF position u_i, the target depth is the smaller of the
# two genotypes' u_i-quantiles, so both marginal depth distributions become their pointwise
# minimum. Cells already at or below target keep everything (p = 1).
filter_strata <- rownames(strat_n)[strat_n[, "KO"] >= MIN_CELLS_FILTER & strat_n[, "WT"] >= MIN_CELLS_FILTER]
say("depth-matching ", length(filter_strata), " strata with >= ", MIN_CELLS_FILTER, " cells per genotype")
p_keep <- rep(1, nrow(md))
for (s in filter_strata) {
  iK <- which(md$stratum == s & md$genotype == "KO"); iW <- which(md$stratum == s & md$genotype == "WT")
  nK <- md$nCount[iK]; nW <- md$nCount[iW]
  for (side in list(list(i = iK, n = nK), list(i = iW, n = nW))) {
    u <- (rank(side$n, ties.method = "first") - 0.5) / length(side$n)
    tgt <- pmin(quantile(nK, u, type = 7, names = FALSE), quantile(nW, u, type = 7, names = FALSE))
    p_keep[side$i] <- pmin(1, tgt / side$n)
  }
}
set.seed(SEED)
thin <- counts
pc <- rep(p_keep, diff(thin@p))
thin@x <- as.numeric(rbinom(length(thin@x), size = as.integer(round(thin@x)), prob = pc))
thin <- Matrix::drop0(thin); rm(pc)
md$nCount_matched <- Matrix::colSums(thin)

auc_rank <- function(score, pos) {              # Mann-Whitney AUC via mid-ranks; constant -> 0.5
  pos <- as.logical(pos); n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score); (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
ksD <- function(a, b) unname(suppressWarnings(ks.test(a, b))$statistic)
depth_tab <- do.call(rbind, lapply(filter_strata, function(s) {
  iK <- md$stratum == s & md$genotype == "KO"; iW <- md$stratum == s & md$genotype == "WT"
  isK <- md$genotype[iK | iW] == "KO"
  data.frame(celltype = sub("\\|.*", "", s), timepoint = sub(".*\\|", "", s),
             n_KO = sum(iK), n_WT = sum(iW),
             median_nCount_KO_raw = median(md$nCount[iK]), median_nCount_WT_raw = median(md$nCount[iW]),
             median_nCount_KO_matched = median(md$nCount_matched[iK]),
             median_nCount_WT_matched = median(md$nCount_matched[iW]),
             ks_D_raw = ksD(md$nCount[iK], md$nCount[iW]),
             ks_D_matched = ksD(md$nCount_matched[iK], md$nCount_matched[iW]),
             auc_depth_raw = auc_rank(log(md$nCount[iK | iW]), isK),
             auc_depth_matched = auc_rank(log(md$nCount_matched[iK | iW]), isK),
             frac_umis_retained = sum(md$nCount_matched[iK | iW]) / sum(md$nCount[iK | iW]))
}))
print(depth_tab[, c("celltype", "timepoint", "median_nCount_KO_raw", "median_nCount_WT_raw",
                    "auc_depth_raw", "auc_depth_matched")], digits = 3)
wtab(depth_tab, "ko_ml_depth_matching.csv")

say("normalising (LogNormalize, 1e4) matched and raw counts ...")
lognorm     <- Seurat::NormalizeData(thin,   normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
lognorm_raw <- Seurat::NormalizeData(counts, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
stopifnot(inherits(lognorm, "dgCMatrix"))
say("normalised  [", elapsed(), "]")

## ---- 3. confounder removal: blocklist, then the ubiquity filter -------------------
grp <- factor(paste(md$stratum, md$genotype, sep = "#"))
G <- Matrix::sparse.model.matrix(~ 0 + grp); colnames(G) <- levels(grp)
n_g <- Matrix::colSums(G)
group_means <- function(X) { m <- sweep(as.matrix(X %*% G), 2, n_g, "/"); colnames(m) <- levels(grp); m }
Mn <- group_means(lognorm)
Dt <- lognorm; Dt@x[] <- 1
Det <- group_means(Dt); rm(Dt)
D_mat <- sapply(filter_strata, function(s) Mn[, paste0(s, "#KO")] - Mn[, paste0(s, "#WT")])
T_mat <- sapply(filter_strata, function(s) Det[, paste0(s, "#KO")] >= DETECT_MIN | Det[, paste0(s, "#WT")] >= DETECT_MIN)
rownames(D_mat) <- rownames(T_mat) <- rownames(lognorm)

ubiquity_rule <- function(D, Tt) {
  m <- rowSums(Tt)
  pos <- rowSums(Tt & D > 0); neg <- rowSums(Tt & D < 0)
  same_sign <- m >= M_MIN & (pos == m | neg == m)
  n_pass <- rowSums(Tt & abs(D) >= DELTA)
  absD <- ifelse(Tt, abs(D), NA); is_cm <- grepl("^Cardiomyocyte\\|", colnames(D))
  cm_abs  <- rowMeans(absD[, is_cm, drop = FALSE], na.rm = TRUE)
  non_med <- apply(absD[, !is_cm, drop = FALSE], 1, median, na.rm = TRUE)
  cm_ratio <- cm_abs / pmax(non_med, 1e-3)
  cm_dominant <- is.finite(cm_ratio) & cm_ratio >= CM_RATIO
  pattern <- same_sign & n_pass >= ceiling(FRAC_K * m)          # the ubiquitous pattern itself
  data.frame(gene = rownames(D), m_testable = m, n_pass_delta = n_pass,
             sign = ifelse(pos == m & m > 0, "+", ifelse(neg == m & m > 0, "-", "mixed")),
             same_sign = same_sign, ubiquitous_pattern = pattern,
             cm_ratio = cm_ratio, cm_dominant = cm_dominant,
             lenient = pattern & !cm_dominant,
             strict  = same_sign & n_pass == m & !cm_dominant,
             median_d = apply(ifelse(Tt, D, NA), 1, median, na.rm = TRUE),
             min_abs_d = suppressWarnings(apply(absD, 1, min, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}
rule_all <- ubiquity_rule(D_mat, T_mat)
rule_all$min_abs_d[!is.finite(rule_all$min_abs_d)] <- NA
# positive control: over the FULL universe the rule must rediscover the blocklist
redisc <- intersect(CONFOUND, rule_all$gene[rule_all$lenient])
testable_conf <- intersect(CONFOUND, rule_all$gene[rule_all$m_testable >= M_MIN])
say("ubiquity filter, positive control: rediscovers blocklist genes {", paste(redisc, collapse = ", "),
    "} of the testable {", paste(testable_conf, collapse = ", "), "}")
if (!"Xist" %in% redisc) warning("[ko_ml] the ubiquity filter did NOT rediscover Xist -- check DELTA/FRAC_K")
# negative control: random WT/WT split with the identical rule
set.seed(SEED + 1L)
wt_i <- which(md$genotype == "WT" & md$stratum %in% filter_strata)
split_lab <- rep(NA_character_, nrow(md))
for (s in filter_strata) { i <- wt_i[md$stratum[wt_i] == s]; split_lab[i] <- sample(rep(c("A", "B"), length.out = length(i))) }
g2 <- factor(paste(md$stratum[wt_i], split_lab[wt_i], sep = "#"))
G2 <- Matrix::sparse.model.matrix(~ 0 + g2); colnames(G2) <- levels(g2)
gm2 <- function(X) { m <- sweep(as.matrix(X %*% G2), 2, Matrix::colSums(G2), "/"); colnames(m) <- levels(g2); m }
Ln_wt <- lognorm[, wt_i]; Mn2 <- gm2(Ln_wt); Dt2 <- Ln_wt; Dt2@x[] <- 1; Det2 <- gm2(Dt2); rm(Dt2, Ln_wt)
D2 <- sapply(filter_strata, function(s) Mn2[, paste0(s, "#A")] - Mn2[, paste0(s, "#B")])
T2 <- sapply(filter_strata, function(s) Det2[, paste0(s, "#A")] >= DETECT_MIN | Det2[, paste0(s, "#B")] >= DETECT_MIN)
rownames(D2) <- rownames(T2) <- rownames(lognorm)
rule_null <- ubiquity_rule(D2, T2)
say("ubiquity filter, negative control (random WT/WT split): removes ", sum(rule_null$lenient),
    " genes lenient, ", sum(rule_null$strict), " strict  (expect ~0)")
rm(Mn2, Det2, D2, T2)

removed <- setdiff(rule_all$gene[rule_all$lenient], CONFOUND)
say("ubiquity filter removes ", length(removed), " genes beyond the ", length(CONFOUND), "-gene blocklist (",
    sum(rule_all$strict[match(removed, rule_all$gene)]), " of them strict)")
cm_dom_kept <- setdiff(rule_all$gene[rule_all$ubiquitous_pattern & rule_all$cm_dominant], CONFOUND)
say("ubiquitous pattern but CM-dominant (>= ", CM_RATIO, "x), KEPT and flagged: ",
    if (length(cm_dom_kept)) paste(sprintf("%s (%.1fx)", cm_dom_kept, rule_all$cm_ratio[match(cm_dom_kept, rule_all$gene)]), collapse = ", ") else "none")
hallmark <- tryCatch({
  h <- msigdbr::msigdbr(species = "Mus musculus", collection = "H")
  unique(h$gene_symbol[h$gs_name == "HALLMARK_E2F_TARGETS"])
}, error = function(e) { say("msigdbr unavailable: ", conditionMessage(e)); character(0) })
shared_up <- tryCatch(readLines(file.path(OUTTAB, "shared_KO_up_P0_and_P7.txt")), error = function(e) character(0))
chr_of <- function(genes) tryCatch({
  x <- suppressWarnings(suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = genes, column = "CHR",
                                              keytype = "SYMBOL", multiVals = "first")))
  unname(x[genes])
}, error = function(e) rep(NA_character_, length(genes)))
gene_class <- function(g) ifelse(grepl("^mt-", g), "mito", ifelse(grepl("^Rp[sl]\\d", g), "ribo",
                          ifelse(grepl("^Hb[ab]-", g), "hemoglobin", ifelse(grepl("^Gm\\d+", g), "Gm", "other"))))
cand <- rule_all[rule_all$same_sign | rule_all$gene %in% CONFOUND, ]
cand$removed <- cand$gene %in% c(CONFOUND, removed)
cand$reason <- ifelse(cand$gene %in% CONFOUND, "blocklist",
                      ifelse(cand$strict, "ubiquitous_strict", ifelse(cand$lenient, "ubiquitous_lenient",
                      ifelse(cand$ubiquitous_pattern & cand$cm_dominant, "candidate_cm_dominant_kept", "candidate_not_removed"))))
cand$chr <- chr_of(cand$gene); cand$is_XY <- cand$chr %in% c("X", "Y")
cand$class <- gene_class(cand$gene)
cand$in_E2F_TARGETS <- cand$gene %in% E2F_TARGETS
cand$in_hallmark_E2F <- cand$gene %in% hallmark
cand$in_shared_KO_up <- cand$gene %in% shared_up
cand$is_E2f7_E2f8 <- cand$gene %in% c("E2f7", "E2f8")
dcols <- as.data.frame(D_mat[cand$gene, , drop = FALSE]); colnames(dcols) <- paste0("d_", gsub("\\|", "_", colnames(D_mat)))
cand$cm_ratio <- round(cand$cm_ratio, 2)
removed_tab <- cbind(cand[, c("gene", "removed", "reason", "m_testable", "n_pass_delta", "sign", "median_d",
                              "min_abs_d", "cm_ratio", "cm_dominant", "chr", "is_XY", "class", "in_E2F_TARGETS",
                              "in_hallmark_E2F", "in_shared_KO_up", "is_E2f7_E2f8")], round(dcols, 4))
removed_tab <- removed_tab[order(!removed_tab$removed, -abs(removed_tab$median_d)), ]
wtab(removed_tab, "ko_ml_removed_genes.csv")
say("removed genes by class: ", paste(names(table(removed_tab$class[removed_tab$removed])),
    table(removed_tab$class[removed_tab$removed]), sep = "=", collapse = "  "))
say("removed genes on chrX/Y: ", sum(removed_tab$is_XY[removed_tab$removed], na.rm = TRUE))
if (any(removed_tab$removed & (removed_tab$in_E2F_TARGETS | removed_tab$is_E2f7_E2f8)))
  say("NOTE: ubiquity filter removed E2F-related gene(s): ",
      paste(removed_tab$gene[removed_tab$removed & (removed_tab$in_E2F_TARGETS | removed_tab$is_E2f7_E2f8)], collapse = ", "))
say("top removed by |median d|: ", paste(head(removed_tab$gene[removed_tab$removed & removed_tab$reason != "blocklist"], 15), collapse = ", "))

EXCLUDE <- union(CONFOUND, removed)
U <- setdiff(rownames(lognorm), EXCLUDE)
say("feature universe: ", length(U), " genes after removing ", length(EXCLUDE))
hvg_of <- function(X) {
  so <- Seurat::CreateSeuratObject(X, min.cells = 0, min.features = 0)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
  v <- Seurat::VariableFeatures(so); rm(so); v
}
hvg <- hvg_of(thin[U, ])
stopifnot(!any(hvg %in% EXCLUDE), length(hvg) == N_HVG)
X_all <- Matrix::t(lognorm[hvg, ]); stopifnot(inherits(X_all, "dgCMatrix"))
if (NAIVE) {
  hvg_naive <- hvg_of(counts)
  X_naive <- Matrix::t(lognorm_raw[hvg_naive, ])
  say("naive control feature space: ", sum(hvg_naive %in% EXCLUDE), " of its ", N_HVG, " HVGs are excluded genes")
}
rm(thin); invisible(gc())
say("features ready  [", elapsed(), "]")

## ---- 4. model helpers --------------------------------------------------------------
yfac <- function(g) factor(g, levels = c("WT", "KO"))          # KO is the modelled (positive) class
fit_ko <- function(X, y, weights = NULL) {
  y <- yfac(y)
  fit <- if (is.null(weights))
    cv.glmnet(X, y, family = "binomial", alpha = ALPHA, nfolds = 5, type.measure = "deviance", keep = TRUE)
  else
    cv.glmnet(X, y, family = "binomial", alpha = ALPHA, nfolds = 5, type.measure = "deviance", keep = TRUE, weights = weights)
  j <- match(fit$lambda.1se, fit$lambda)
  oof <- fit$fit.preval[, j]                                     # link scale
  stopifnot(!anyNA(oof))
  list(fit = fit, j = j, oof_link = oof, cv_auc = auc_rank(oof, y == "KO"), nzero = unname(fit$nzero[j]))
}
coef_1se <- function(fit) {
  b <- as.matrix(coef(fit, s = "lambda.1se")); b <- b[rownames(b) != "(Intercept)", 1]; b[b != 0]
}
bal_sample <- function(iK, iW, n) c(sample(iK, n), sample(iW, n))
run_tasks <- function(tasks, fn) {
  go <- function(i) { set.seed(SEED + 1000L + i); fn(tasks[[i]], i) }
  res <- if (NCORES > 1L) mclapply(seq_along(tasks), go, mc.cores = NCORES, mc.set.seed = FALSE, mc.preschedule = FALSE)
         else lapply(seq_along(tasks), go)
  err <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(err)) stop("[ko_ml] ", sum(err), " task(s) failed; first: ",
                     paste(as.character(res[[which(err)[1]]]), collapse = " "))
  res
}
qlo <- function(x) unname(quantile(x, 0.025)); qhi <- function(x) unname(quantile(x, 0.975))
ml_strata <- rownames(strat_n)[strat_n[, "KO"] >= MIN_CELLS_ML & strat_n[, "WT"] >= MIN_CELLS_ML]
ct_of <- function(s) sub("\\|.*", "", s); tp_of <- function(s) sub(".*\\|", "", s)

## ---- 5. analysis 1: where does the KO act? -----------------------------------------
say("analysis 1: separability in ", length(ml_strata), " strata, ", R_REP, " real + ", N_PERM, " null fits each ...")
tasks <- do.call(c, lapply(ml_strata, function(s)
  c(lapply(seq_len(R_REP), function(r) list(s = s, kind = "real", rep = r)),
    lapply(seq_len(N_PERM), function(r) list(s = s, kind = "null", rep = r)))))
res <- run_tasks(tasks, function(t, i) {
  iK <- which(md$stratum == t$s & md$genotype == "KO"); iW <- which(md$stratum == t$s & md$genotype == "WT")
  n <- min(CAP_SEP, length(iK), length(iW)); sel <- bal_sample(iK, iW, n)
  y <- md$genotype[sel]; if (t$kind == "null") y <- sample(y)
  r <- fit_ko(X_all[sel, , drop = FALSE], y)
  data.frame(stratum = t$s, kind = t$kind, rep = t$rep, n_per_genotype = n, auc = r$cv_auc, nzero = r$nzero)
})
sep_raw <- do.call(rbind, res)
sep <- do.call(rbind, lapply(ml_strata, function(s) {
  a <- sep_raw$auc[sep_raw$stratum == s & sep_raw$kind == "real"]; z <- sep_raw$auc[sep_raw$stratum == s & sep_raw$kind == "null"]
  data.frame(celltype = ct_of(s), timepoint = tp_of(s),
             n_KO = strat_n[s, "KO"], n_WT = strat_n[s, "WT"], n_per_genotype = sep_raw$n_per_genotype[sep_raw$stratum == s][1],
             auc_mean = mean(a), auc_lo = qlo(a), auc_hi = qhi(a),
             null_mean = mean(z), null_lo = qlo(z), null_hi = qhi(z), null_sd = sd(z),
             delta = mean(a) - mean(z), z_score = (mean(a) - mean(z)) / max(sd(z), 1e-6),
             frac_null_ge = mean(z >= mean(a)), mean_nzero = mean(sep_raw$nzero[sep_raw$stratum == s & sep_raw$kind == "real"]),
             role = ifelse(ct_of(s) == "RBC", "ambient_control", "tissue"))
}))
sep <- sep[order(-sep$auc_mean), ]
print(sep[, c("celltype", "timepoint", "n_per_genotype", "auc_mean", "null_mean", "delta", "mean_nzero")], digits = 3)
wtab(sep, "ko_ml_separability.csv")
# Nulls a little BELOW 0.5 are expected: cross-validation on exchangeable labels is pessimistic
# (each training fold's class imbalance is anti-correlated with its held-out fold). A null
# ABOVE 0.5 is the leakage direction and is the one to worry about.
if (any(sep$null_mean > 0.53)) warning("[ko_ml] a permutation null sits ABOVE 0.53 -- check for leakage (a second lane?)")
sep$label <- ifelse(sep$role == "ambient_control", paste0(sep$celltype, "\n(ambient ctrl)"), sep$celltype)
ord <- sep$label[order(-sep$auc_mean)]; sep$label <- factor(sep$label, levels = unique(ord))
p1 <- ggplot(sep, aes(label, auc_mean, colour = timepoint)) +
  geom_hline(yintercept = 0.5, colour = "grey70") +
  geom_linerange(aes(ymin = null_lo, ymax = null_hi), colour = "grey55", position = position_dodge(width = .6), linewidth = 3, alpha = .4) +
  geom_linerange(aes(ymin = auc_lo, ymax = auc_hi), position = position_dodge(width = .6)) +
  geom_point(size = 2.6, position = position_dodge(width = .6)) +
  scale_colour_manual(values = TP_COLS) + coord_cartesian(ylim = c(0.4, 1)) +
  labs(title = "KO-vs-WT separability per cell type (elastic net, equal n, depth-matched, sex genes removed)",
       subtitle = sprintf("points: cross-validated AUC over %d subsamples; grey: label-permutation null (%d). Descriptive, n = 1 animal per condition.", R_REP, N_PERM),
       x = NULL, y = "out-of-fold AUC", colour = "timepoint") +
  theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9), axis.text.x = element_text(size = 9))
ggsave(file.path(OUTFIG, "ko_ml_separability.png"), p1, width = 9, height = 5, dpi = 120)
say("analysis 1 done  [", elapsed(), "]")

## ---- 6. analysis 2: does the signature transfer across animals (P0 <-> P7)? --------
tr_ct <- Reduce(intersect, list(ct_of(ml_strata[tp_of(ml_strata) == "P0"]), ct_of(ml_strata[tp_of(ml_strata) == "P7"])))
say("analysis 2: transfer for ", paste(tr_ct, collapse = ", "), " ...")
tasks <- do.call(c, lapply(tr_ct, function(ct) do.call(c, lapply(c("P0->P7", "P7->P0"), function(dir)
  c(lapply(seq_len(R_REP_TR), function(r) list(ct = ct, dir = dir, kind = "real", rep = r)),
    lapply(seq_len(N_PERM_TR), function(r) list(ct = ct, dir = dir, kind = "null", rep = r)))))))
res <- run_tasks(tasks, function(t, i) {
  src <- substr(t$dir, 1, 2); tgt <- substr(t$dir, 5, 6)
  iK <- which(md$celltype == t$ct & md$timepoint == src & md$genotype == "KO")
  iW <- which(md$celltype == t$ct & md$timepoint == src & md$genotype == "WT")
  n <- min(CAP_TRAIN, length(iK), length(iW)); sel <- bal_sample(iK, iW, n)
  y <- md$genotype[sel]; if (t$kind == "null") y <- sample(y)
  r <- fit_ko(X_all[sel, , drop = FALSE], y)
  it <- which(md$celltype == t$ct & md$timepoint == tgt)
  link <- as.vector(predict(r$fit, X_all[it, , drop = FALSE], s = "lambda.1se", type = "link"))
  data.frame(celltype = t$ct, direction = t$dir, kind = t$kind, rep = t$rep, n_train_per_genotype = n,
             n_test = length(it), within_cv_auc = r$cv_auc, transfer_auc = auc_rank(link, md$genotype[it] == "KO"), nzero = r$nzero)
})
tr_raw <- do.call(rbind, res)
tr <- do.call(rbind, lapply(split(tr_raw, list(tr_raw$celltype, tr_raw$direction), drop = TRUE), function(d) {
  a <- d[d$kind == "real", ]; z <- d[d$kind == "null", ]
  data.frame(celltype = d$celltype[1], direction = d$direction[1], n_train_per_genotype = a$n_train_per_genotype[1], n_test = a$n_test[1],
             within_cv_auc = mean(a$within_cv_auc), within_cv_lo = qlo(a$within_cv_auc), within_cv_hi = qhi(a$within_cv_auc),
             transfer_auc = mean(a$transfer_auc), transfer_lo = qlo(a$transfer_auc), transfer_hi = qhi(a$transfer_auc),
             null_transfer_mean = mean(z$transfer_auc), null_transfer_lo = qlo(z$transfer_auc), null_transfer_hi = qhi(z$transfer_auc),
             transfer_minus_null = mean(a$transfer_auc) - mean(z$transfer_auc), mean_nzero = mean(a$nzero))
}))
rownames(tr) <- NULL
print(tr[, c("celltype", "direction", "within_cv_auc", "transfer_auc", "null_transfer_mean", "mean_nzero")], digits = 3)
wtab(tr, "ko_ml_transfer.csv")
trl <- rbind(data.frame(tr[, c("celltype", "direction")], what = "within-timepoint CV", auc = tr$within_cv_auc, lo = tr$within_cv_lo, hi = tr$within_cv_hi),
             data.frame(tr[, c("celltype", "direction")], what = "transfer to other timepoint", auc = tr$transfer_auc, lo = tr$transfer_lo, hi = tr$transfer_hi),
             data.frame(tr[, c("celltype", "direction")], what = "permuted-training null", auc = tr$null_transfer_mean, lo = tr$null_transfer_lo, hi = tr$null_transfer_hi))
trl$what <- factor(trl$what, levels = c("within-timepoint CV", "transfer to other timepoint", "permuted-training null"))
p2 <- ggplot(trl, aes(celltype, auc, colour = what)) +
  geom_hline(yintercept = 0.5, colour = "grey70") +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .6)) +
  geom_point(size = 2.4, position = position_dodge(width = .6)) +
  facet_wrap(~ direction) + scale_colour_manual(values = c("#37474f", "#d84315", "grey60")) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(title = "Does a KO-vs-WT signature learned in one animal pair transfer to the other?",
       subtitle = "Train at one timepoint, score every cell of that type at the other. Both pairs are male-KO / female-WT: transfer cannot break the sex confound.",
       x = NULL, y = "AUC", colour = NULL) +
  theme(legend.position = "bottom", plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTFIG, "ko_ml_transfer.png"), p2, width = 9, height = 5.5, dpi = 120)
say("analysis 2 done  [", elapsed(), "]")

## ---- 7. analysis 3: E2F targets vs expression-matched random gene sets -------------
fs_ct <- if (identical(FS_CT, "all")) tr_ct else intersect(strsplit(FS_CT, ",")[[1]], tr_ct)
cc <- tryCatch(cc_lists(), error = function(e) { say("cc_lists() failed: ", conditionMessage(e)); list(S = character(0), G2M = character(0)) })
SETS <- list(E2F_TARGETS = E2F_TARGETS, HALLMARK_E2F_TARGETS = hallmark, S_phase = cc$S, G2M_phase = cc$G2M,
             CM_MATURE = CM_MATURE, CM_IMMATURE = CM_IMMATURE)
RAND_N <- c(E2F_TARGETS = N_RAND_E2F, HALLMARK_E2F_TARGETS = N_RAND_HALLMARK)
say("analysis 3: feature sets for ", paste(fs_ct, collapse = ", "), " -- sizes in universe: ",
    paste(names(SETS), vapply(SETS, function(g) length(intersect(g, U)), integer(1)), sep = "=", collapse = " "))
set.seed(SEED + 2L)
CTX <- list()
for (ct in fs_ct) for (mode in c("P0 CV", "P7 CV", "P0->P7", "P7->P0")) {
  src <- substr(mode, 1, 2); tgt <- if (grepl("CV", mode)) NA else substr(mode, 5, 6)
  iK <- which(md$celltype == ct & md$timepoint == src & md$genotype == "KO")
  iW <- which(md$celltype == ct & md$timepoint == src & md$genotype == "WT")
  n <- min(CAP_FS, length(iK), length(iW)); sel <- bal_sample(iK, iW, n)
  Xtr <- Matrix::t(lognorm[U, sel]); ytr <- md$genotype[sel]
  it <- if (is.na(tgt)) integer(0) else which(md$celltype == ct & md$timepoint == tgt)
  Xtg <- if (length(it)) Matrix::t(lognorm[U, it]) else NULL
  det <- Matrix::colMeans(Xtr > 0); expressed <- colnames(Xtr)[det >= DETECT_MIN]
  mu <- Matrix::colMeans(Xtr[, expressed, drop = FALSE])
  bins <- cut(mu, breaks = unique(quantile(mu, seq(0, 1, by = 0.05))), include.lowest = TRUE, labels = FALSE)
  names(bins) <- expressed
  CTX[[paste(ct, mode)]] <- list(ct = ct, mode = mode, src = src, tgt = tgt, Xtr = Xtr, ytr = ytr, Xtg = Xtg,
                                 ytg = if (length(it)) md$genotype[it] else NULL, n = n, bins = bins,
                                 d = D_mat[U, paste(ct, src, sep = "|")])
}
draw_matched <- function(target, bins) {
  tb <- bins[intersect(target, names(bins))]; pool_bins <- bins[setdiff(names(bins), target)]
  unlist(lapply(split(names(tb), tb), function(g) {
    pool <- names(pool_bins)[pool_bins == bins[g[1]]]
    if (length(pool) >= length(g)) sample(pool, length(g)) else c(pool, sample(names(pool_bins), length(g) - length(pool)))
  }), use.names = FALSE)
}
eval_set <- function(cx, genes) {
  genes <- intersect(genes, colnames(cx$Xtr)); if (length(genes) < 2L) return(NULL)
  r <- fit_ko(cx$Xtr[, genes, drop = FALSE], cx$ytr)
  auc <- if (is.null(cx$Xtg)) r$cv_auc else
    auc_rank(as.vector(predict(r$fit, cx$Xtg[, genes, drop = FALSE], s = "lambda.1se", type = "link")), cx$ytg == "KO")
  b <- coef_1se(r$fit)
  list(auc = auc, n_used = length(genes), nzero = length(b), n_pos = sum(b > 0), n_neg = sum(b < 0),
       mean_coef = if (length(b)) mean(b) else 0, frac_KO_up = mean(cx$d[genes] > 0), mean_d = mean(cx$d[genes]))
}
tasks <- list()
for (k in names(CTX)) {
  for (sn in names(SETS)) tasks[[length(tasks) + 1]] <- list(k = k, set = sn, kind = "named", rep = 0L)
  tasks[[length(tasks) + 1]] <- list(k = k, set = "HVG_2000", kind = "named", rep = 0L)
  for (sn in names(RAND_N)) for (r in seq_len(RAND_N[[sn]])) tasks[[length(tasks) + 1]] <- list(k = k, set = sn, kind = "random", rep = r)
}
res <- run_tasks(tasks, function(t, i) {
  cx <- CTX[[t$k]]
  genes <- if (t$set == "HVG_2000") hvg else SETS[[t$set]]
  if (t$kind == "random") genes <- draw_matched(intersect(genes, names(cx$bins)), cx$bins)
  e <- eval_set(cx, genes)
  if (is.null(e)) return(data.frame(k = t$k, set = t$set, kind = t$kind, rep = t$rep, auc = NA, n_used = 0L, nzero = NA, n_pos = NA, n_neg = NA, mean_coef = NA, frac_KO_up = NA, mean_d = NA))
  data.frame(k = t$k, set = t$set, kind = t$kind, rep = t$rep, auc = e$auc, n_used = e$n_used, nzero = e$nzero, n_pos = e$n_pos,
             n_neg = e$n_neg, mean_coef = e$mean_coef, frac_KO_up = e$frac_KO_up, mean_d = e$mean_d)
})
fs_raw <- do.call(rbind, res)
fs <- do.call(rbind, lapply(names(CTX), function(k) {
  cx <- CTX[[k]]; nm <- fs_raw[fs_raw$k == k & fs_raw$kind == "named", ]
  do.call(rbind, lapply(seq_len(nrow(nm)), function(i) {
    rn <- fs_raw$auc[fs_raw$k == k & fs_raw$kind == "random" & fs_raw$set == nm$set[i]]
    data.frame(celltype = cx$ct, context = cx$mode, train_timepoint = cx$src, test_timepoint = ifelse(is.na(cx$tgt), cx$src, cx$tgt),
               evaluation = ifelse(is.na(cx$tgt), "5-fold CV", "transfer"), n_train_per_genotype = cx$n,
               set_name = nm$set[i], n_in_set = if (nm$set[i] == "HVG_2000") length(hvg) else length(SETS[[nm$set[i]]]),
               n_used = nm$n_used[i], auc = nm$auc[i],
               null_n = length(rn), null_mean = if (length(rn)) mean(rn) else NA, null_p95 = if (length(rn)) unname(quantile(rn, .95)) else NA,
               percentile = if (length(rn)) mean(rn <= nm$auc[i]) else NA,
               nzero = nm$nzero[i], n_pos_coef = nm$n_pos[i], n_neg_coef = nm$n_neg[i], mean_coef = nm$mean_coef[i],
               frac_genes_KO_up = nm$frac_KO_up[i], mean_d_KO_minus_WT = nm$mean_d[i])
  }))
}))
print(fs[, c("context", "set_name", "n_used", "auc", "null_mean", "percentile", "frac_genes_KO_up")], digits = 3)
wtab(fs, "ko_ml_featuresets.csv")
rnd <- fs_raw[fs_raw$kind == "random" & !is.na(fs_raw$auc), ]
rnd$context <- vapply(rnd$k, function(k) CTX[[k]]$mode, character(1)); rnd$celltype <- vapply(rnd$k, function(k) CTX[[k]]$ct, character(1))
lines <- fs[fs$set_name %in% c("E2F_TARGETS", "HALLMARK_E2F_TARGETS", "HVG_2000") & !is.na(fs$auc), ]
p3 <- ggplot(rnd, aes(auc, fill = set)) + geom_histogram(bins = 30, alpha = .55, position = "identity") +
  geom_vline(data = lines, aes(xintercept = auc, colour = set_name, linetype = set_name), linewidth = .9) +
  facet_grid(celltype ~ context, scales = "free_y") +
  scale_fill_manual(values = c(E2F_TARGETS = "#90a4ae", HALLMARK_E2F_TARGETS = "#bcaaa4"), name = "random sets matched to") +
  scale_colour_manual(values = c(E2F_TARGETS = "#c62828", HALLMARK_E2F_TARGETS = "#6a1b9a", HVG_2000 = "black"), name = NULL) +
  scale_linetype_manual(values = c(E2F_TARGETS = "solid", HALLMARK_E2F_TARGETS = "solid", HVG_2000 = "dashed"), name = NULL) +
  labs(title = "Do E2F-target genes carry more KO-vs-WT information than expression-matched random gene sets?",
       subtitle = "Histograms: random sets of the same size and expression profile. Lines: the named set; dashed black = all 2,000 HVGs (ceiling).",
       x = "AUC", y = "random sets") +
  guides(fill = guide_legend(nrow = 1, order = 1), colour = guide_legend(nrow = 1, order = 2), linetype = guide_legend(nrow = 1, order = 2)) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5))
ggsave(file.path(OUTFIG, "ko_ml_featuresets.png"), p3, width = 10, height = 3 + 2.2 * length(fs_ct), dpi = 120)
rm(CTX); invisible(gc())
say("analysis 3 done  [", elapsed(), "]")

## ---- 8. analysis 4a: the shipped maturation model as a KO phenotype readout --------
cm_i <- which(md$celltype == "Cardiomyocyte")
mat_ok <- FALSE
bundle_f <- file.path(RESULTS, "models", "cm_maturation_glmnet.rds")
if (file.exists(bundle_f)) {
  bundle <- readRDS(bundle_f)
  stopifnot(identical(bundle$positive, "P7"), identical(bundle$family, "binomial"))
  b_st <- as.matrix(coef(bundle$model, s = "lambda.1se"))[, 1]
  neut <- b_st[intersect(names(b_st), EXCLUDE)]; neut <- neut[neut != 0]
  say("stage model: neutralising ", length(neut), " excluded-gene weights: ",
      paste(names(neut), round(neut, 3), sep = "=", collapse = "  "))
  stage_link <- function(L) {
    Xn <- matrix(0, length(cm_i), length(bundle$features), dimnames = list(rownames(md)[cm_i], bundle$features))
    common <- setdiff(intersect(rownames(L), bundle$features), EXCLUDE)
    Xn[, common] <- as.matrix(Matrix::t(L[common, cm_i]))
    as.vector(predict(bundle$model, Xn, s = "lambda.1se", type = "link"))
  }
  mat <- tryCatch({
    lk_m <- stage_link(lognorm); lk_r <- stage_link(lognorm_raw)
    do.call(rbind, lapply(c("P0", "P7"), function(tp) do.call(rbind, lapply(c("matched", "raw"), function(dep) {
      lk <- if (dep == "matched") lk_m else lk_r
      sel <- md$timepoint[cm_i] == tp; g <- md$genotype[cm_i][sel]; l <- lk[sel]
      q <- function(x, p) unname(quantile(x, p))
      data.frame(timepoint = tp, depth = dep, n_KO = sum(g == "KO"), n_WT = sum(g == "WT"),
                 median_logit_KO = median(l[g == "KO"]), median_logit_WT = median(l[g == "WT"]),
                 iqr_logit_KO = IQR(l[g == "KO"]), iqr_logit_WT = IQR(l[g == "WT"]),
                 KO_minus_WT_median_logit = median(l[g == "KO"]) - median(l[g == "WT"]),
                 auc_KO_vs_WT = auc_rank(l, g == "KO"), cliff_delta = 2 * auc_rank(l, g == "KO") - 1,
                 frac_P7_KO = mean(plogis(l[g == "KO"]) > 0.5), frac_P7_WT = mean(plogis(l[g == "WT"]) > 0.5),
                 median_pP7_KO = median(plogis(l[g == "KO"])), median_pP7_WT = median(plogis(l[g == "WT"])),
                 neutralised_weights = paste(names(neut), collapse = ";"))
    }))))
  }, error = function(e) { say("stage model could not be applied: ", conditionMessage(e)); NULL })
  if (!is.null(mat)) {
    mat_ok <- TRUE
    print(mat[, c("timepoint", "depth", "median_logit_KO", "median_logit_WT", "cliff_delta", "frac_P7_KO", "frac_P7_WT")], digits = 3)
    wtab(mat, "ko_ml_maturation_shift.csv")
    vd <- rbind(data.frame(logit = lk_m, depth = "matched", timepoint = md$timepoint[cm_i], genotype = md$genotype[cm_i]),
                data.frame(logit = lk_r, depth = "raw",     timepoint = md$timepoint[cm_i], genotype = md$genotype[cm_i]))
    vd$genotype <- yfac(vd$genotype)
    p4 <- ggplot(vd, aes(genotype, logit, fill = genotype)) + geom_violin(alpha = .7, colour = NA) +
      geom_boxplot(width = .15, outlier.shape = NA, fill = "white") +
      facet_grid(depth ~ timepoint) + scale_fill_manual(values = GT_COLS) +
      labs(title = "Shipped P0-vs-P7 stage model applied to KO vs WT cardiomyocytes",
           subtitle = sprintf("logit P(P7); sex/ambient weights neutralised (%s). 'matched' = depth-matched counts, 'raw' = as sequenced.", paste(names(neut), collapse = ", ")),
           x = NULL, y = "maturation logit (higher = more P7-like)") +
      theme(legend.position = "none", plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5))
    ggsave(file.path(OUTFIG, "ko_ml_maturation_shift.png"), p4, width = 8, height = 6, dpi = 120)
  }
} else say("no shipped stage model at ", bundle_f, "; skipping analysis 4a")
say("analysis 4a done  [", elapsed(), "]")

## ---- 9. analysis 4b + 5: full CM models, P(KO) by state, gene panels ---------------
say("full cardiomyocyte models per timepoint (class-weighted, all lane1 CMs) ...")
full_tasks <- list(list(tp = "P0", space = "controlled"), list(tp = "P7", space = "controlled"))
if (NAIVE) full_tasks <- c(full_tasks, list(list(tp = "P0", space = "naive"), list(tp = "P7", space = "naive")))
res <- run_tasks(full_tasks, function(t, i) {
  ii <- cm_i[md$timepoint[cm_i] == t$tp]; y <- yfac(md$genotype[ii])
  w <- ifelse(y == "KO", 0.5 / sum(y == "KO"), 0.5 / sum(y == "WT")) * length(y)
  X <- if (t$space == "controlled") X_all[ii, , drop = FALSE] else X_naive[ii, , drop = FALSE]
  r <- fit_ko(X, y, weights = w)
  list(tp = t$tp, space = t$space, idx = ii, oof_link = r$oof_link, cv_auc = r$cv_auc, coef = coef_1se(r$fit), nzero = r$nzero)
})
names(res) <- vapply(res, function(r) paste(r$space, r$tp), character(1))
for (r in res) say(sprintf("  %-10s %s: %d CMs, OOF AUC %.3f, %d non-zero genes", r$space, r$tp, length(r$idx), r$cv_auc, r$nzero))

# subcluster labels: production variant cm_dims30_res0.2 == dims 30 arm of the per-cell sweep table
sub_lab <- rep(NA_character_, nrow(md))
pc_f <- file.path(OUTTAB, "pcdims_cm_percell.csv.gz")
if (file.exists(pc_f)) {
  pc <- read.csv(gzfile(pc_f), stringsAsFactors = FALSE); pc <- pc[pc$dims == 30, ]
  lab <- paste0("CM", pc$SCT_snn_res.0.2); names(lab) <- pc$cell
  sub_lab <- unname(lab[rownames(md)])
  cov <- mean(!is.na(sub_lab[cm_i]))
  say(sprintf("production subcluster labels (dims30_res0.2) cover %.1f%% of lane1 CMs", 100 * cov))
  if (cov < 0.9) warning("[ko_ml] subcluster labels cover < 90% of CMs -- check pcdims_cm_percell.csv.gz")
  rec_f <- file.path(OUTTAB, "cm_subcluster_dims30_res0.2_cellcycle.csv")
  if (file.exists(rec_f)) {
    rec <- read.csv(rec_f); ref_n <- tapply(rec$n, rec$cm_subcluster, sum); ours <- table(lab)
    common <- intersect(names(ref_n), names(ours))
    say("  reconciliation with cm_subcluster_dims30_res0.2_cellcycle.csv (all lanes): ",
        sum(ref_n[common] == ours[common]), "/", length(common), " subcluster counts identical")
  }
} else say("pcdims_cm_percell.csv.gz not found; subcluster breakdown skipped")

pko_rows <- list(); pko_plot <- list()
for (tp in c("P0", "P7")) {
  r <- res[[paste("controlled", tp)]]; ii <- r$idx
  d <- data.frame(cell = rownames(md)[ii], timepoint = tp, genotype = yfac(md$genotype[ii]), Phase = md$Phase[ii],
                  cm_subcluster = sub_lab[ii], logit = r$oof_link, p_ko = plogis(r$oof_link))
  say(sprintf("  %s: median P(KO) KO cells %.3f, WT cells %.3f", tp, median(d$p_ko[d$genotype == "KO"]), median(d$p_ko[d$genotype == "WT"])))
  gmed <- tapply(d$logit, d$genotype, median)
  for (st in c("Phase", "cm_subcluster")) {
    dd <- d[!is.na(d[[st]]), ]
    for (s in sort(unique(dd[[st]]))) {
      ds <- dd[dd[[st]] == s, ]
      wa <- if (all(table(ds$genotype) >= 30)) auc_rank(ds$logit, ds$genotype == "KO") else NA_real_
      for (g in c("WT", "KO")) {
        x <- ds[ds$genotype == g, ]
        pko_rows[[length(pko_rows) + 1]] <- data.frame(timepoint = tp, state_type = st, state = s, genotype = g, n = nrow(x),
          median_logit = median(x$logit), median_p_ko = median(x$p_ko), frac_p_ko_gt_0.5 = mean(x$p_ko > 0.5),
          rel_median_logit = median(x$logit) - gmed[[g]], within_state_auc = wa)
      }
    }
  }
  pko_plot[[tp]] <- d
}
pko <- do.call(rbind, pko_rows); pko <- pko[!is.na(pko$median_logit), ]
wtab(pko, "ko_ml_pko_by_state.csv")
ext <- pko[pko$state_type == "cm_subcluster" & pko$n >= MIN_STATE_N, ]; ext <- ext[order(-abs(ext$rel_median_logit)), ]
say("most KO-shifted subcluster states (|rel median logit|): ",
    paste(head(sprintf("%s %s %s (%+.2f, n=%d)", ext$timepoint, ext$state, ext$genotype, ext$rel_median_logit, ext$n), 6), collapse = "; "))
pd <- do.call(rbind, pko_plot)
pA <- if (any(!is.na(pd$cm_subcluster))) {
  pd2 <- pd[!is.na(pd$cm_subcluster), ]
  keep <- names(which(table(pd2$cm_subcluster) >= 2 * MIN_STATE_N)); pd2 <- pd2[pd2$cm_subcluster %in% keep, ]
  pd2$cm_subcluster <- factor(pd2$cm_subcluster, levels = unique(pd2$cm_subcluster[order(as.integer(sub("CM", "", pd2$cm_subcluster)))]))
  ggplot(pd2, aes(cm_subcluster, logit, fill = genotype)) + geom_boxplot(outlier.size = .3, outlier.alpha = .2, linewidth = .3) +
    facet_wrap(~ timepoint, ncol = 1, scales = "free_x") + scale_fill_manual(values = GT_COLS) +
    labs(title = sprintf("Out-of-fold KO logit by production CM subcluster (states with >= %d cells)", 2 * MIN_STATE_N), x = NULL, y = "OOF logit P(KO)") +
    theme(legend.position = "bottom", plot.title = element_text(size = 11))
} else NULL
pB <- ggplot(pd, aes(Phase, logit, fill = genotype)) + geom_boxplot(outlier.size = .3, outlier.alpha = .2, linewidth = .3) +
  facet_wrap(~ timepoint, ncol = 1) + scale_fill_manual(values = GT_COLS) +
  labs(title = "... and by cell-cycle phase", x = NULL, y = NULL) + theme(legend.position = "none", plot.title = element_text(size = 11))
p5 <- if (is.null(pA)) pB else (pA + pB + plot_layout(widths = c(3, 1)))
ggsave(file.path(OUTFIG, "ko_ml_pko_by_state.png"), p5, width = 11, height = 6.5, dpi = 120)

# gene panels
read_de <- function(tp) {
  f <- file.path(OUTTAB, sprintf("percelltype_%s_Cardiomyocyte_KOvsWT.descriptive.DE.csv", tp))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE); d[!duplicated(d$gene), c("gene", "log2FoldChange", "padj")]
}
de0 <- read_de("P0"); de7 <- read_de("P7")
build_panel <- function(c0, c7, space) {
  genes <- union(names(c0), names(c7)); if (!length(genes)) return(NULL)
  pnl <- data.frame(gene = genes, weight_P0 = unname(c0[genes]), weight_P7 = unname(c7[genes]), stringsAsFactors = FALSE)
  pnl$weight_P0[is.na(pnl$weight_P0)] <- 0; pnl$weight_P7[is.na(pnl$weight_P7)] <- 0
  pnl$nonzero_P0 <- pnl$weight_P0 != 0; pnl$nonzero_P7 <- pnl$weight_P7 != 0
  pnl$in_both <- pnl$nonzero_P0 & pnl$nonzero_P7
  pnl$same_sign <- pnl$in_both & sign(pnl$weight_P0) == sign(pnl$weight_P7)
  pnl$direction <- ifelse(pnl$weight_P0 + pnl$weight_P7 > 0, "up_in_KO", "down_in_KO")
  if (!is.null(de0)) { m <- match(pnl$gene, de0$gene); pnl$DE_log2FC_P0 <- de0$log2FoldChange[m]; pnl$DE_padj_P0 <- de0$padj[m] }
  if (!is.null(de7)) { m <- match(pnl$gene, de7$gene); pnl$DE_log2FC_P7 <- de7$log2FoldChange[m]; pnl$DE_padj_P7 <- de7$padj[m] }
  pnl$in_shared_KO_up <- pnl$gene %in% shared_up
  pnl$in_CONFOUND <- pnl$gene %in% CONFOUND
  pnl$removed_by_ubiquity <- pnl$gene %in% removed
  pnl$ubiquitous_cm_dominant <- pnl$gene %in% cm_dom_kept
  pnl$in_E2F_TARGETS <- pnl$gene %in% E2F_TARGETS
  pnl$in_hallmark_E2F <- pnl$gene %in% hallmark
  pnl$class <- gene_class(pnl$gene)
  pnl$chr <- chr_of(pnl$gene)
  pnl$feature_space <- space
  pnl[order(-pmax(abs(pnl$weight_P0), abs(pnl$weight_P7))), ]
}
panel <- build_panel(res[["controlled P0"]]$coef, res[["controlled P7"]]$coef, "controlled")
stopifnot(!any(panel$in_CONFOUND), !any(panel$removed_by_ubiquity))
agree <- function(pnl) {
  a <- c(if ("DE_log2FC_P0" %in% names(pnl)) sign(pnl$weight_P0[pnl$nonzero_P0]) == sign(pnl$DE_log2FC_P0[pnl$nonzero_P0]),
         if ("DE_log2FC_P7" %in% names(pnl)) sign(pnl$weight_P7[pnl$nonzero_P7]) == sign(pnl$DE_log2FC_P7[pnl$nonzero_P7]))
  mean(a, na.rm = TRUE)
}
say(sprintf("controlled CM panel: %d genes (%d P0, %d P7, %d in both, %d same sign); sign agreement with pseudobulk DESeq2 log2FC = %.2f",
            nrow(panel), sum(panel$nonzero_P0), sum(panel$nonzero_P7), sum(panel$in_both), sum(panel$same_sign), agree(panel)))
say("  top by |weight|: ", paste(head(panel$gene, 15), collapse = ", "))
say("  E2F targets in panel: ", paste(panel$gene[panel$in_E2F_TARGETS | panel$in_hallmark_E2F], collapse = ", "))
say("  ubiquitous-but-CM-dominant genes in panel: ", paste(panel$gene[panel$ubiquitous_cm_dominant], collapse = ", "))
wtab(panel, "ko_ml_cm_panel.csv")
if (NAIVE) {
  naive <- build_panel(res[["naive P0"]]$coef, res[["naive P7"]]$coef, "naive_no_controls")
  say(sprintf("naive CM panel (no controls): %d genes, %d blocklist, %d ubiquity-removed; sign agreement with DESeq2 = %.2f",
              nrow(naive), sum(naive$in_CONFOUND), sum(naive$removed_by_ubiquity), agree(naive)))
  say("  naive top by |weight|: ", paste(head(naive$gene, 15), collapse = ", "))
  wtab(naive, "ko_ml_cm_panel_naive.csv")
}

## ---- 10. summary ---------------------------------------------------------------------
cat("\n=== SUMMARY (descriptive; n = 1 animal per condition; genotype == sex == library) ===\n")
cat(sprintf("  controls: lane1 only (%d cells); blocklist %d + ubiquity-removed %d genes; depth matched in %d strata\n",
            nrow(md), length(CONFOUND), length(removed), length(filter_strata)))
top <- sep[order(-sep$delta), ]
cat(sprintf("  separability: highest %s %s (AUC %.2f vs null %.2f); lowest %s %s (AUC %.2f); RBC ambient control %s\n",
            top$celltype[1], top$timepoint[1], top$auc_mean[1], top$null_mean[1],
            top$celltype[nrow(top)], top$timepoint[nrow(top)], top$auc_mean[nrow(top)],
            if (any(top$role == "ambient_control")) sprintf("AUC %.2f", top$auc_mean[top$role == "ambient_control"][1]) else "n/a"))
cm_tr <- tr[tr$celltype == "Cardiomyocyte", ]
if (nrow(cm_tr)) cat(sprintf("  CM transfer: %s AUC %.2f (within-CV %.2f, null %.2f); %s AUC %.2f (within-CV %.2f, null %.2f)\n",
                             cm_tr$direction[1], cm_tr$transfer_auc[1], cm_tr$within_cv_auc[1], cm_tr$null_transfer_mean[1],
                             cm_tr$direction[2], cm_tr$transfer_auc[2], cm_tr$within_cv_auc[2], cm_tr$null_transfer_mean[2]))
e2 <- fs[fs$set_name == "E2F_TARGETS" & fs$celltype == "Cardiomyocyte", ]
for (i in seq_len(nrow(e2))) cat(sprintf("  E2F_TARGETS %-6s AUC %.2f, percentile vs matched random %.2f, %.0f%% of genes higher in KO\n",
                                          e2$context[i], e2$auc[i], e2$percentile[i], 100 * e2$frac_genes_KO_up[i]))
if (mat_ok) for (i in which(mat$depth == "matched")) cat(sprintf("  maturation shift %s (matched depth): Cliff's delta %+.2f (KO vs WT logit)\n", mat$timepoint[i], mat$cliff_delta[i]))
cat(sprintf("  total runtime %s\n", elapsed()))
cat("=== DONE ko_signature_ml ===\n")
