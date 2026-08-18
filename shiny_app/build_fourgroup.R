# build_fourgroup.R
# ---------------------------------------------------------------------------
# Four-group (WT-P0 / WT-P7 / KO-P0 / KO-P7) analysis within CM subclusters,
# written back into app_data.rds as app$fourgroup. Drives the two new
# "Four-group (WT/KO x P0/P7)" and "Maturation n P7 KO" tabs.
#
#   1. counts      per-subcluster x four-group cell counts / percentages, with an
#                  explicit status flag on under-powered arms
#   2. phase       per-subcluster x four-group cell-cycle phase composition
#                  (the G1-proportion question)
#   3. scores      maturation / metabolic / proliferation score summaries per
#                  subcluster x group, computed overall AND within G1 only
#   4. de          descriptive Wilcoxon DE for four contrasts x two phase strata
#   5. maturation  gene-level maturation association (mature- vs immature-linked)
#   6. intersect   maturation axis x P7 KO-vs-WT, quadrant-assigned
#
# Also adds a `cm_subcluster` column to app$meta / app$deg_meta ("CM2", ... ;
# NA for non-CM cells) so CM subcluster becomes a first-class filter in the
# existing "Subset & DEGs" tab and a colour/split option on the UMAP.
#
# Run LOCALLY, in a real R session (NOT the bash sandbox — it segfaults on the
# ~100 MB gzip rds), from the repo root or from shiny_app/:
#     Rscript shiny_app/build_fourgroup.R
#
# Run AFTER build_signature_scores.R — this script consumes its sig_* columns.
#
# FEASIBILITY GATE: run once with --probe to print the per-subcluster x group
# cell counts, the phase composition, gene coverage, and an estimate of how much
# this will add to the bundle, then EXIT WITHOUT WRITING:
#     Rscript shiny_app/build_fourgroup.R --probe
#
# Matrix choice (--probe reports which was used):
#     --matrix=broad     app$deg_expr  — many genes, ~8k downsampled cells (default)
#     --matrix=curated   app$expr      — 2,181 curated genes, all ~30k cells
#     --seurat=<path>    a Seurat .rds — all genes x all CM cells (best; needs Seurat)
#
# TWO CONFOUNDS THIS SCRIPT HANDLES EXPLICITLY
#   * Sort artefact. P7 was FACS cycling-enriched 4.5-5.2x and P0 essentially
#     unenriched (model/README.md), so a raw P0-vs-P7 contrast reads out the sort,
#     not development. Every contrast is therefore computed in TWO strata: "G1"
#     (phase-matched, the default the app shows) and "all" (raw, labelled).
#   * Maturation circularity. sig_maturation's immature program contains Mki67 /
#     Top2a / Ccnd1, so it is partly a cell-cycle score. The maturation axis here
#     defaults to sig_maturation_nocc (those three dropped); see
#     build_signature_scores.R.
#
# DESCRIPTIVE ONLY. n = 1 animal per genotype x timepoint, so no contrast has
# biological replication. Wilcoxon is run cell-level; the p/padj columns are
# pseudoreplicated and are carried for ranking only. Tables sort by effect size.
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))

ARGS   <- commandArgs(trailingOnly = TRUE)
PROBE  <- "--probe" %in% ARGS
argval <- function(flag, default = "") {
  hit <- grep(paste0("^", flag, "="), ARGS, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
MATRIX_PICK <- argval("--matrix", "broad")
SEURAT_PATH <- argval("--seurat", "")

if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

RES        <- "0.2"                    # the only resolution the app exposes
MIN_CELLS  <- 10                       # hard floor, same as deg_compute() in app.R
THIN_CELLS <- 50                       # soft floor: computed, but flagged as thin
MIN_PCT    <- 5                        # row gate: expressed in >=5% of one side
MIN_LFC    <- 0.5                      # row gate: OR |log2FC| >= this
MAX_PADJ   <- 0.05                     # row gate: OR padj < this
TERTILE_Q  <- c(1/3, 2/3)              # maturation-axis split
MAT_SCORE  <- "sig_maturation_nocc"    # cycle-free maturation axis (see header)
GROUPS     <- c("WT-P0", "WT-P7", "KO-P0", "KO-P7")
# Genes that bypass the row gate, so a table never silently lacks a gene someone
# went looking for. The shortlist is from the collaborator's email; the G2/M
# members (Foxm1, Aurkb, Birc5, Prc1) fall below the gate in the G1 stratum by
# construction, and "absent" would read as "not measured" rather than "not
# expressed in G1". Their pct_A/pct_B columns carry the real answer.
KEEP_ALWAYS <- c("Birc5","Foxm1","Rrm2","Aurkb","Prc1","Gabbr2","Tcf4","Adamts9",
                 "E2f7","E2f8","Mki67","Top2a","Ect2","Myh6","Myh7","Nppa","Tnni3")

# The four contrasts. A/B are four-group labels; logFC > 0 means "up in A".
# The first two are the email's within-genotype temporal comparisons; the third
# is its headline and the only sort-clean contrast; the fourth is not requested
# but is required to answer "is this change P7-specific?".
CONTRASTS <- list(
  list(key = "WT_P0_vs_P7", label = "WT: P0 vs P7", A = "WT-P7", B = "WT-P0",
       pos = "up at P7", neg = "up at P0", xlab = "log2 fold change (P7 / P0)"),
  list(key = "KO_P0_vs_P7", label = "KO: P0 vs P7", A = "KO-P7", B = "KO-P0",
       pos = "up at P7", neg = "up at P0", xlab = "log2 fold change (P7 / P0)"),
  list(key = "P7_KO_vs_WT", label = "P7: KO vs WT", A = "KO-P7", B = "WT-P7",
       pos = "up in KO", neg = "up in WT", xlab = "log2 fold change (KO / WT)"),
  list(key = "P0_KO_vs_WT", label = "P0: KO vs WT", A = "KO-P0", B = "WT-P0",
       pos = "up in KO", neg = "up in WT", xlab = "log2 fold change (KO / WT)"))
STRATA <- list(G1 = "G1", all = NULL)   # NULL = no phase restriction

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
stopifnot(!is.null(app$meta), !is.null(app$cm$meta))
CONF <- app$confound

# ---- CM subcluster label, and the four-group variable ----------------------
subcol <- paste0("SCT_snn_res.", RES)
stopifnot(subcol %in% names(app$cm$meta))
cmm <- app$cm$meta
cmm$cm_subcluster <- paste0("CM", cmm[[subcol]])
SUBS <- { u <- unique(cmm$cm_subcluster); u[order(as.integer(sub("CM", "", u)))] }
CLUSTERS <- c("AllCM", SUBS)           # "AllCM" = every CM cell, pooled

fourgrp <- function(df) factor(paste(df$genotype, df$timepoint, sep = "-"), levels = GROUPS)
cmm$fourgrp <- fourgrp(cmm)
cm_lookup <- setNames(cmm$cm_subcluster, cmm$cell)

# ---- pick the expression matrix -------------------------------------------
# Returns list(M = genes x cells matrix, MD = metadata aligned to colnames(M), name).
pick_matrix <- function() {
  if (nzchar(SEURAT_PATH)) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("--seurat given but the Seurat package is not installed in this R.")
    cat("Reading Seurat object:", SEURAT_PATH, "...\n")
    so <- readRDS(SEURAT_PATH)
    keep <- colnames(so)[colnames(so) %in% names(cm_lookup)]
    if (!length(keep)) stop("No cell barcodes in the Seurat object match app$cm$meta$cell.")
    M  <- Seurat::GetAssayData(so, assay = "RNA", layer = "data")[, keep, drop = FALSE]
    MD <- cmm[match(keep, cmm$cell), , drop = FALSE]
    return(list(M = M, MD = MD, name = paste0("seurat:", basename(SEURAT_PATH))))
  }
  if (MATRIX_PICK == "broad" && !is.null(app$deg_expr)) {
    M  <- app$deg_expr
    MD <- if (!is.null(app$deg_meta)) app$deg_meta else app$meta
    MD <- MD[match(colnames(M), MD$cell), , drop = FALSE]
    return(list(M = M, MD = MD, name = "broad (app$deg_expr)"))
  }
  if (MATRIX_PICK == "broad")
    cat("!! --matrix=broad requested but app$deg_expr is absent — falling back to the curated panel.\n")
  M  <- app$expr
  MD <- app$meta[match(colnames(M), app$meta$cell), , drop = FALSE]
  list(M = M, MD = MD, name = "curated (app$expr)")
}

mx <- pick_matrix(); M <- mx$M; MD <- mx$MD
MD$cm_subcluster <- unname(cm_lookup[MD$cell])
MD$fourgrp <- fourgrp(MD)
MD <- MD[, intersect(c("cell","genotype","timepoint","Phase","cm_subcluster","fourgrp",
                       "sig_maturation", MAT_SCORE), names(MD)), drop = FALSE]
is_cm <- !is.na(MD$cm_subcluster)

cat(sprintf("\nMatrix: %s  (%d genes x %d cells; %d are cardiomyocytes)\n",
            mx$name, nrow(M), ncol(M), sum(is_cm)))
if (grepl("^curated", mx$name))
  cat("!! Curated panel only — DE is restricted to", nrow(M),
      "genes and is NOT a discovery-scale DEG table.\n")

# ---- 1. counts: cluster x four-group --------------------------------------
# Counted on ALL CM cells (app$cm$meta), not the DE matrix, so the percentages
# the collaborator asked for describe the dataset rather than the downsample.
count_rows <- list()
for (cl in CLUSTERS) {
  sel  <- if (cl == "AllCM") rep(TRUE, nrow(cmm)) else cmm$cm_subcluster == cl
  tot  <- sum(sel)
  n_de <- table(factor(MD$fourgrp[is_cm & (if (cl == "AllCM") TRUE else MD$cm_subcluster == cl)],
                       levels = GROUPS))
  for (g in GROUPS) {
    n    <- sum(sel & cmm$fourgrp == g)
    n_g1 <- sum(sel & cmm$fourgrp == g & as.character(cmm$Phase) == "G1")
    nd   <- as.integer(n_de[[g]])
    # Two separate ways an arm can be too thin, and the G1 one bites hardest:
    # the app defaults to the phase-matched stratum, so an arm can clear the
    # floor overall and still be a handful of cells once restricted to G1.
    # CM2's KO-P0 arm is the worked example — 31 cells, ~12 of them G1.
    st <- if (nd < MIN_CELLS) sprintf("too few cells (<%d in DE matrix)", MIN_CELLS)
          else if (n_g1 < THIN_CELLS) sprintf("thin in G1 (%d cells)", n_g1)
          else "ok"
    count_rows[[length(count_rows) + 1]] <- data.frame(
      cluster = cl, group = g, n = n, n_G1 = n_g1,
      pct_of_cluster = if (tot) round(100 * n / tot, 1) else NA_real_,
      pct_of_group   = round(100 * n / sum(cmm$fourgrp == g), 1),
      n_in_de_matrix = nd, status = st, stringsAsFactors = FALSE)
  }
}
counts <- do.call(rbind, count_rows)

# ---- 2. phase: cluster x four-group x Phase --------------------------------
phase_rows <- list()
if ("Phase" %in% names(cmm)) {
  plevs <- c("G1", "S", "G2M")
  for (cl in CLUSTERS) {
    sel <- if (cl == "AllCM") rep(TRUE, nrow(cmm)) else cmm$cm_subcluster == cl
    for (g in GROUPS) {
      s <- sel & cmm$fourgrp == g; tot <- sum(s)
      for (ph in plevs) {
        n <- sum(s & as.character(cmm$Phase) == ph)
        phase_rows[[length(phase_rows) + 1]] <- data.frame(
          cluster = cl, group = g, Phase = ph, n = n,
          pct = if (tot) round(100 * n / tot, 1) else NA_real_,
          n_group = tot, stringsAsFactors = FALSE)
      }
    }
  }
}
phase <- if (length(phase_rows)) do.call(rbind, phase_rows) else NULL

# ---- 3. scores: cluster x four-group x score, overall and within G1 --------
SCORE_COLS <- intersect(
  c("sig_maturation", "sig_maturation_nocc", "sig_mat_mature", "sig_mat_immature",
    "sig_mat_immature_nocc", "sig_metabolic", "sig_prolif", "sig_cytokinesis", "sig_ccexit"),
  names(cmm))
score_rows <- list()
for (sc in SCORE_COLS) for (cl in CLUSTERS) for (g in GROUPS) for (st in c("all", "G1")) {
  sel <- (if (cl == "AllCM") rep(TRUE, nrow(cmm)) else cmm$cm_subcluster == cl) & cmm$fourgrp == g
  if (st == "G1") sel <- sel & as.character(cmm$Phase) == "G1"
  v <- cmm[[sc]][sel]; v <- v[!is.na(v)]
  if (!length(v)) next
  score_rows[[length(score_rows) + 1]] <- data.frame(
    cluster = cl, group = g, score = sc, stratum = st, n = length(v),
    mean = round(mean(v), 4), median = round(stats::median(v), 4),
    sd = round(stats::sd(v), 4),
    se = round(stats::sd(v) / sqrt(length(v)), 4), stringsAsFactors = FALSE)
}
scores <- if (length(score_rows)) do.call(rbind, score_rows) else NULL

if (length(SCORE_COLS) == 0)
  cat("!! No sig_* score columns found — run build_signature_scores.R first.\n")
if (!MAT_SCORE %in% SCORE_COLS)
  cat("!! ", MAT_SCORE, " missing — re-run build_signature_scores.R (it now emits the\n",
      "   cycle-free maturation score). Falling back to sig_maturation for the axis.\n", sep = "")
MAT_USE <- if (MAT_SCORE %in% names(MD)) MAT_SCORE else "sig_maturation"

# ---- shared DE core --------------------------------------------------------
# Descriptive cell-level Wilcoxon (presto), gated to keep the bundle small.
de_one <- function(idxA, idxB) {
  nA <- length(idxA); nB <- length(idxB)
  if (nA < MIN_CELLS || nB < MIN_CELLS) return(NULL)
  cols <- c(idxA, idxB)
  grp  <- rep(c("A", "B"), c(nA, nB))
  res  <- presto::wilcoxauc(M[, cols, drop = FALSE], grp)
  rA <- res[res$group == "A", ]; rB <- res[res$group == "B", ]
  m  <- match(rA$feature, rB$feature)
  d <- data.frame(
    gene = rA$feature,
    log2FoldChange = round(rA$logFC, 4),
    auc = round(rA$auc, 4),
    pvalue = signif(rA$pval, 3), padj = signif(rA$padj, 3),
    pct_A = round(rA$pct_in, 1), pct_B = round(rA$pct_out, 1),
    mean_A = round(rA$avgExpr, 4), mean_B = round(rB$avgExpr[m], 4),
    confounder = rA$feature %in% CONF,
    n_A = nA, n_B = nB, stringsAsFactors = FALSE)
  keep <- (pmax(d$pct_A, d$pct_B) >= MIN_PCT &
           (d$padj < MAX_PADJ | abs(d$log2FoldChange) >= MIN_LFC)) |
          d$gene %in% KEEP_ALWAYS
  d <- d[which(keep), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d[order(-abs(d$log2FoldChange)), , drop = FALSE]
}

# ---- 4. the DE grid: cluster x contrast x stratum --------------------------
de <- list(); skipped <- list()
if (!PROBE) {
  if (!requireNamespace("presto", quietly = TRUE))
    stop("presto is required: remotes::install_github('immunogenomics/presto')")
  cat("\n== four-group DE (", length(CLUSTERS), " clusters x ", length(CONTRASTS),
      " contrasts x ", length(STRATA), " strata) ==\n", sep = "")
  for (cl in CLUSTERS) {
    inCl <- is_cm & (if (cl == "AllCM") TRUE else MD$cm_subcluster == cl)
    de[[cl]] <- list()
    for (ct in CONTRASTS) for (sn in names(STRATA)) {
      ph  <- STRATA[[sn]]
      base <- inCl & (if (is.null(ph)) TRUE else as.character(MD$Phase) == ph)
      idxA <- which(base & MD$fourgrp == ct$A)
      idxB <- which(base & MD$fourgrp == ct$B)
      key  <- paste0(ct$key, "__", sn)
      d <- de_one(idxA, idxB)
      if (is.null(d)) {
        skipped[[length(skipped) + 1]] <- data.frame(
          cluster = cl, contrast = ct$key, stratum = sn,
          n_A = length(idxA), n_B = length(idxB),
          reason = if (length(idxA) < MIN_CELLS || length(idxB) < MIN_CELLS)
                     "too few cells" else "no genes passed the row gate",
          stringsAsFactors = FALSE)
      } else de[[cl]][[key]] <- d
    }
    cat(sprintf("  %-6s %d tables\n", cl, length(de[[cl]])))
  }
}
skipped <- if (length(skipped)) do.call(rbind, skipped) else NULL

# ---- 5. maturation association --------------------------------------------
# Rank every gene along the maturation axis by comparing top- vs bottom-tertile
# cells on MAT_USE. Done WITHIN each timepoint and averaged, so the axis is a
# maturation axis and not a restatement of P0-vs-P7.
maturation <- NULL
if (!PROBE && MAT_USE %in% names(MD)) {
  cat("\n== maturation association (", MAT_USE, ", tertile split within timepoint) ==\n", sep = "")
  per_tp <- list()
  for (tp in unique(as.character(MD$timepoint[is_cm]))) {
    sel <- is_cm & as.character(MD$timepoint) == tp & !is.na(MD[[MAT_USE]])
    v   <- MD[[MAT_USE]][sel]
    if (sum(sel) < 3 * MIN_CELLS) next
    qs  <- stats::quantile(v, TERTILE_Q, na.rm = TRUE)
    idxHi <- which(sel)[v >= qs[2]]; idxLo <- which(sel)[v <= qs[1]]
    d <- de_one(idxHi, idxLo)          # logFC > 0 => higher in MATURE cells
    if (is.null(d)) next
    per_tp[[tp]] <- data.frame(gene = d$gene, lfc = d$log2FoldChange,
                               auc = d$auc, padj = d$padj, stringsAsFactors = FALSE)
    cat(sprintf("  %-3s  %d mature-high vs %d mature-low cells, %d genes\n",
                tp, length(idxHi), length(idxLo), nrow(d)))
  }
  if (length(per_tp)) {
    genes_all <- Reduce(union, lapply(per_tp, function(x) x$gene))
    grab <- function(x, col) x[[col]][match(genes_all, x$gene)]
    lfcm <- do.call(cbind, lapply(per_tp, grab, col = "lfc"))
    aucm <- do.call(cbind, lapply(per_tp, grab, col = "auc"))
    padm <- do.call(cbind, lapply(per_tp, grab, col = "padj"))
    # curated program membership, for cross-checking the data-driven axis
    prog_mature   <- c("Myh6","Tnni3","Pln","Atp2a2","Ckm","Myl2","Cox6a2","Ckmt2","Actn2","Csrp3")
    prog_immature <- c("Myh7","Tnni1","Nppa","Nppb","Myl7","Actc1")
    maturation <- data.frame(
      gene = genes_all,
      mat_log2FC = round(rowMeans(lfcm, na.rm = TRUE), 4),
      mat_auc    = round(rowMeans(aucm, na.rm = TRUE), 4),
      mat_padj   = apply(padm, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)),
      n_timepoints = rowSums(!is.na(lfcm)),
      stringsAsFactors = FALSE)
    maturation$mat_class <- with(maturation, ifelse(
      is.na(mat_padj) | mat_padj >= MAX_PADJ | abs(mat_log2FC) < MIN_LFC, "ns",
      ifelse(mat_log2FC > 0, "mature-associated", "immature-associated")))
    maturation$in_curated_program <- ifelse(
      maturation$gene %in% prog_mature, "mature",
      ifelse(maturation$gene %in% prog_immature, "immature", NA_character_))
    maturation <- maturation[order(-abs(maturation$mat_log2FC)), ]
    cat(sprintf("  %d genes ranked; %d mature-associated, %d immature-associated\n",
                nrow(maturation), sum(maturation$mat_class == "mature-associated"),
                sum(maturation$mat_class == "immature-associated")))
  }
}

# ---- 6. intersection: maturation axis x P7 KO-vs-WT ------------------------
intersect_tab <- NULL
if (!PROBE && !is.null(maturation)) {
  rows <- list()
  for (cl in CLUSTERS) {
    d <- de[[cl]][["P7_KO_vs_WT__G1"]]
    stratum <- "G1"
    if (is.null(d)) { d <- de[[cl]][["P7_KO_vs_WT__all"]]; stratum <- "all" }
    if (is.null(d)) next
    m <- match(d$gene, maturation$gene)
    ok <- !is.na(m)
    if (!any(ok)) next
    r <- data.frame(
      cluster = cl, gene = d$gene[ok], stratum = stratum,
      mat_log2FC = maturation$mat_log2FC[m[ok]],
      mat_class  = maturation$mat_class[m[ok]],
      p7ko_log2FC = d$log2FoldChange[ok], p7ko_padj = d$padj[ok],
      p7ko_pct_KO = d$pct_A[ok], p7ko_pct_WT = d$pct_B[ok],
      confounder = d$confounder[ok], stringsAsFactors = FALSE)
    r$quadrant <- with(r, ifelse(
      mat_class == "immature-associated" & p7ko_log2FC > 0, "immature_up_in_KO",
      ifelse(mat_class == "mature-associated" & p7ko_log2FC < 0, "mature_down_in_KO",
      ifelse(mat_class == "immature-associated" & p7ko_log2FC < 0, "immature_down_in_KO",
      ifelse(mat_class == "mature-associated"   & p7ko_log2FC > 0, "mature_up_in_KO", "ns")))))
    rows[[length(rows) + 1]] <- r
  }
  if (length(rows)) {
    intersect_tab <- do.call(rbind, rows)
    intersect_tab <- intersect_tab[order(intersect_tab$cluster,
                                         -abs(intersect_tab$p7ko_log2FC)), ]
    cat("\n== intersection ==\n")
    print(table(intersect_tab$quadrant))
  }
}

# ---- probe report ----------------------------------------------------------
if (PROBE) {
  cat("\n== cell counts: subcluster x four-group (all CM cells) ==\n")
  wide <- reshape(counts[, c("cluster","group","n")], idvar = "cluster",
                  timevar = "group", direction = "wide")
  names(wide) <- sub("^n\\.", "", names(wide))
  print(wide, row.names = FALSE)
  cat("\n== cell counts IN THE DE MATRIX (", mx$name, ") ==\n", sep = "")
  wide2 <- reshape(counts[, c("cluster","group","n_in_de_matrix")], idvar = "cluster",
                   timevar = "group", direction = "wide")
  names(wide2) <- sub("^n_in_de_matrix\\.", "", names(wide2))
  print(wide2, row.names = FALSE)
  bad <- counts[counts$status != "ok", c("cluster","group","n","n_in_de_matrix")]
  if (nrow(bad)) {
    cat("\n!! arms below the ", MIN_CELLS, "-cell DE floor (these contrasts will be skipped):\n", sep = "")
    print(bad, row.names = FALSE)
  }
  if (!is.null(phase)) {
    cat("\n== G1 percentage: subcluster x four-group ==\n")
    g1 <- phase[phase$Phase == "G1", c("cluster","group","pct")]
    wg <- reshape(g1, idvar = "cluster", timevar = "group", direction = "wide")
    names(wg) <- sub("^pct\\.", "", names(wg))
    print(wg, row.names = FALSE)
  }
  cat("\n== score columns available ==\n")
  cat(if (length(SCORE_COLS)) paste(SCORE_COLS, collapse = ", ") else
      "NONE — run build_signature_scores.R first", "\n")
  cat("\n== estimated DE tables ==\n")
  est <- 0
  for (cl in CLUSTERS) for (ct in CONTRASTS) for (sn in names(STRATA)) {
    inCl <- is_cm & (if (cl == "AllCM") TRUE else MD$cm_subcluster == cl)
    base <- inCl & (if (is.null(STRATA[[sn]])) TRUE else as.character(MD$Phase) == STRATA[[sn]])
    if (sum(base & MD$fourgrp == ct$A) >= MIN_CELLS &&
        sum(base & MD$fourgrp == ct$B) >= MIN_CELLS) est <- est + 1
  }
  cat(sprintf("  %d of %d possible tables will be computed (rest skipped: too few cells)\n",
              est, length(CLUSTERS) * length(CONTRASTS) * length(STRATA)))
  # Measure rather than guess: build one real table and scale it. AllCM is the
  # widest contrast, so this over-estimates the per-cluster tables.
  if (requireNamespace("presto", quietly = TRUE)) {
    b  <- is_cm & as.character(MD$Phase) == "G1"
    ex <- de_one(which(b & MD$fourgrp == "KO-P7"), which(b & MD$fourgrp == "WT-P7"))
    if (!is.null(ex)) {
      tf <- tempfile(); saveRDS(ex, tf, compress = "gzip")
      cat(sprintf("  sample table (AllCM P7 KO-vs-WT, G1): %d rows, %.2f MB gzip\n",
                  nrow(ex), file.size(tf) / 1e6))
      cat(sprintf("  => app$fourgroup adds roughly %.1f MB gzip to the bundle\n",
                  est * file.size(tf) / 1e6))
      unlink(tf)
    }
  }
  cat("\n  NOTE: the bundle is re-saved with compress=\"gzip\" (house convention, see\n")
  cat("  README). If the current file is xz-compressed, it will grow ~1.8x on its own,\n")
  cat("  independently of anything this script adds. Check with:  xxd -l 6 app_data.rds\n")
  cat("  (1f8b = gzip, fd377a58 = xz)\n")
  cat("\n--probe: nothing written.\n")
  quit(save = "no")
}

# ---- assemble + write ------------------------------------------------------
app$fourgroup <- list(
  built = list(
    when = as.character(Sys.time()), matrix = mx$name,
    n_genes = nrow(M), n_cells_de = sum(is_cm), n_cells_total = nrow(cmm),
    res = RES, groups = GROUPS,
    contrasts = do.call(rbind, lapply(CONTRASTS, function(x)
      data.frame(key = x$key, label = x$label, A = x$A, B = x$B,
                 pos = x$pos, neg = x$neg, xlab = x$xlab, stringsAsFactors = FALSE))),
    strata = names(STRATA), maturation_score = MAT_USE,
    gate = sprintf("max(pct) >= %g%% AND (padj < %g OR |log2FC| >= %g)",
                   MIN_PCT, MAX_PADJ, MIN_LFC),
    keep_always = KEEP_ALWAYS,
    min_cells = MIN_CELLS, thin_cells = THIN_CELLS),
  counts = counts, phase = phase, scores = scores,
  de = de, skipped = skipped,
  maturation = maturation, intersect = intersect_tab)

# CM subcluster as a first-class metadata column, so the existing "Subset & DEGs"
# tab can filter to CM2 and the UMAP can colour/split by subcluster.
app$meta$cm_subcluster <- unname(cm_lookup[app$meta$cell])
if (!is.null(app$deg_meta)) app$deg_meta$cm_subcluster <- unname(cm_lookup[app$deg_meta$cell])
app$cm$meta$cm_subcluster <- cmm$cm_subcluster

cat("\n== summary ==\n")
cat(sprintf("  DE tables      : %d\n", sum(vapply(de, length, 0L))))
cat(sprintf("  skipped        : %d\n", if (is.null(skipped)) 0L else nrow(skipped)))
cat(sprintf("  maturation     : %s genes\n", if (is.null(maturation)) "0" else nrow(maturation)))
cat(sprintf("  intersection   : %s rows\n", if (is.null(intersect_tab)) "0" else nrow(intersect_tab)))
cat(sprintf("  cm_subcluster  : added to meta (%d CM cells labelled)\n",
            sum(!is.na(app$meta$cm_subcluster))))

cat("\nBacking up -> app_data.pre_fourgroup.bak.rds\n")
file.copy("app_data.rds", "app_data.pre_fourgroup.bak.rds", overwrite = TRUE)
cat("Saving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
