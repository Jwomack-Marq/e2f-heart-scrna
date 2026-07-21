# build_signature_scores.R
# ---------------------------------------------------------------------------
# Precompute per-cell module scores for two analyses and write them back into
# app_data.rds as new columns on app$meta (all cells) and app$cm$meta (CM cells),
# plus a definitions/coverage table at app$score_meta.
#
#   1. Cardiomyocyte cycle-exit / polyploidization signatures
#        sig_prolif       proliferation (G2/M program)
#        sig_cytokinesis  cytokinesis machinery
#        sig_ccexit       cell-cycle exit / maturation arrest
#        sig_ploidy       polyploidization proxy = prolif - cytokinesis
#                         (cycling WITHOUT cytokinesis -> binucleation/endoreduplication)
#   2. Maturation + metabolic-switch scoring
#        sig_maturation   CM maturation = mature - immature program
#        sig_metabolic    metabolic maturation = FAO/OXPHOS - glycolysis
#      (component scores sig_glycolysis / sig_faox / sig_mat_mature / sig_mat_immature
#       are also stored for the dedicated tab, but not surfaced in the UMAP dropdown.)
#
# Run LOCALLY, in a real R session (NOT the bash sandbox — it segfaults on the
# ~100 MB gzip rds):
#     & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" shiny_app/build_signature_scores.R
#   (run from the repo root, or cd into shiny_app first — paths below are relative
#    to shiny_app/).
#
# FEASIBILITY GATE: run once with --probe to print per-signature gene coverage in
# the curated panel (app$expr, full ~30k cells) vs the broad matrix (app$deg_expr,
# ~8k cells) and EXIT WITHOUT WRITING. Decide whether any set is too sparse, then
# re-run without --probe to compute + save.
#     & "...\Rscript.exe" shiny_app/build_signature_scores.R --probe
#
# Method: an AddModuleScore-equivalent (mean set expression minus a matched-
# expression-bin control set) in base R + Matrix — no Seurat. Sex/construct
# confounder genes (app$confound) are dropped from every set first; none of the
# curated sets below contain confounders, so the scores are confounder-free by
# construction (no runtime toggle needed).
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))

PROBE <- "--probe" %in% commandArgs(trailingOnly = TRUE)

if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

# ---- gene sets (mouse symbols) --------------------------------------------
# Reuse the app's curated cell-cycle programs where possible; extend with the
# cytokinesis / cell-cycle-exit / metabolic panels the app does not yet carry.
SETS <- list(
  prolif        = c("Mki67","Top2a","Ccnb1","Ccnb2","Cdk1","Cdc20","Aurka","Aurkb",
                    "Bub1","Birc5","Cenpa","Cenpe","Cenpf","Ube2c","Cks2","Nusap1","Tpx2"),
  cytokinesis   = c("Anln","Ect2","Racgap1","Kif23","Cit","Aurkb","Kif20b","Prc1",
                    "Cep55","Mklp1","Sept7","Sept9","Cdca8","Incenp"),
  ccexit        = c("Cdkn1a","Cdkn1c","Cdkn2a","Cdkn1b","Meis1","Rb1","Btg2","Gadd45a"),
  mat_mature    = c("Myh6","Tnni3","Pln","Atp2a2","Ckm","Myl2","Cox6a2","Ckmt2","Actn2","Csrp3"),
  mat_immature  = c("Myh7","Tnni1","Nppa","Nppb","Ccnd1","Mki67","Top2a","Myl7","Actc1"),
  glycolysis    = c("Slc2a1","Hk1","Hk2","Pfkm","Pfkl","Pkm","Ldha","Gapdh","Eno1","Aldoa","Pgk1"),
  faox          = c("Cpt1a","Cpt1b","Cpt2","Acadm","Acadvl","Acadl","Hadha","Hadhb",
                    "Ppargc1a","Cox6a2","Ndufa4","Sdha","Acaa2","Etfa"))

# difference-score groups: two component sets scored on the SAME matrix so the net
# (pos - neg) is well-defined; the components are stored too (some are headline scores).
GROUPS <- list(
  list(net = "sig_ploidy",     pos = "prolif",     neg = "cytokinesis",
       comp = c(sig_prolif = "prolif", sig_cytokinesis = "cytokinesis")),
  list(net = "sig_maturation", pos = "mat_mature", neg = "mat_immature",
       comp = c(sig_mat_mature = "mat_mature", sig_mat_immature = "mat_immature")),
  list(net = "sig_metabolic",  pos = "faox",       neg = "glycolysis",
       comp = c(sig_faox = "faox", sig_glycolysis = "glycolysis")))
# stand-alone single-set scores
SINGLES <- c(sig_ccexit = "ccexit")

cat("Loading app_data.rds ...\n")
app  <- readRDS("app_data.rds")
EXPRc <- app$expr        # curated panel (genes x full cells)
EXPRb <- app$deg_expr    # broad matrix (genes x ~8k cells)
CONF  <- app$confound
stopifnot(!is.null(EXPRc), !is.null(app$meta))
# drop confounders from every set up front
SETS <- lapply(SETS, function(g) setdiff(g, CONF))

# ---- coverage report ------------------------------------------------------
cov1 <- function(g, m) if (is.null(m)) 0L else length(intersect(g, rownames(m)))
cat("\n== gene coverage per set (curated app$expr | broad app$deg_expr) ==\n")
for (nm in names(SETS)) {
  g <- SETS[[nm]]
  ne <- cov1(g, EXPRc); nd <- cov1(g, EXPRb)
  miss <- setdiff(g, union(rownames(EXPRc), if (is.null(EXPRb)) character(0) else rownames(EXPRb)))
  cat(sprintf("  %-13s n=%2d  curated=%2d  broad=%2d  missing(both)=%s\n",
              nm, length(g), ne, nd, if (length(miss)) paste(miss, collapse = ",") else "-"))
}
if (PROBE) { cat("\n--probe: coverage only, nothing written.\n"); quit(save = "no") }

# ---- AddModuleScore-equivalent -------------------------------------------
module_score <- function(mat, features, nbin = 24, ctrl = 100, seed = 1) {
  features <- intersect(features, rownames(mat))
  if (!length(features)) return(NULL)
  avg  <- Matrix::rowMeans(mat)
  brks <- quantile(avg, probs = seq(0, 1, length.out = nbin + 1), na.rm = TRUE)
  brks <- unique(brks); brks[1] <- -Inf; brks[length(brks)] <- Inf
  bins <- cut(avg, breaks = brks, labels = FALSE, include.lowest = TRUE)
  names(bins) <- rownames(mat)
  set.seed(seed)
  ctrl_genes <- character(0)
  for (f in features) {
    pool <- names(bins)[bins == bins[[f]]]
    n <- min(ctrl, length(pool))
    if (n > 0) ctrl_genes <- c(ctrl_genes, sample(pool, n))
  }
  feat_mean <- Matrix::colMeans(mat[features, , drop = FALSE])
  sc <- if (length(ctrl_genes)) feat_mean - Matrix::colMeans(mat[ctrl_genes, , drop = FALSE]) else feat_mean
  setNames(as.numeric(sc), colnames(mat))
}
# pick the matrix that covers `features` best; ties -> curated (full cells)
pick_mat <- function(features) {
  ne <- cov1(features, EXPRc); nd <- cov1(features, EXPRb)
  if (nd > ne) list(mat = EXPRb, name = "broad") else list(mat = EXPRc, name = "curated")
}

meta_cells <- app$meta$cell
cm_cells   <- if (!is.null(app$cm$meta)) app$cm$meta$cell else NULL
score_meta <- list()
put <- function(colname, vec, setname, matname, n_used, n_set) {
  app$meta[[colname]]  <<- vec[meta_cells]
  if (!is.null(cm_cells)) app$cm$meta[[colname]] <<- vec[cm_cells]
  score_meta[[length(score_meta) + 1]] <<- data.frame(
    score = colname, sets = setname, matrix = matname,
    n_genes_used = n_used, n_genes_set = n_set,
    n_cells_scored = sum(!is.na(vec)), stringsAsFactors = FALSE)
}

cat("\n== computing scores ==\n")
# stand-alone single-set scores
for (col in names(SINGLES)) {
  g <- SETS[[SINGLES[[col]]]]; pm <- pick_mat(g)
  v <- module_score(pm$mat, g)
  if (is.null(v)) { cat(sprintf("  %-16s SKIP (no genes in either matrix)\n", col)); next }
  put(col, v, SINGLES[[col]], pm$name, cov1(g, pm$mat), length(g))
  cat(sprintf("  %-16s %-8s genes=%d/%d\n", col, pm$name, cov1(g, pm$mat), length(g)))
}
# difference-score groups — both components on the SAME matrix, each computed once
for (gr in GROUPS) {
  pos <- SETS[[gr$pos]]; neg <- SETS[[gr$neg]]
  pm  <- pick_mat(c(pos, neg))
  comps <- setNames(list(module_score(pm$mat, pos), module_score(pm$mat, neg)), c(gr$pos, gr$neg))
  if (is.null(comps[[gr$pos]]) || is.null(comps[[gr$neg]])) {
    cat(sprintf("  %-16s SKIP (a component has no genes)\n", gr$net)); next }
  # store each component under its exposed column name
  for (cn in names(gr$comp)) {
    setname <- gr$comp[[cn]]; v <- comps[[setname]]
    put(cn, v, setname, pm$name, cov1(SETS[[setname]], pm$mat), length(SETS[[setname]]))
  }
  # store the net (pos - neg)
  net <- setNames(comps[[gr$pos]] - comps[[gr$neg]], names(comps[[gr$pos]]))
  put(gr$net, net, paste(gr$pos, "-", gr$neg), pm$name,
      cov1(c(pos, neg), pm$mat), length(c(pos, neg)))
  cat(sprintf("  %-16s %-8s (%s − %s)\n", gr$net, pm$name, gr$pos, gr$neg))
}

app$score_meta <- do.call(rbind, score_meta)
# columns to expose in the UMAP "colour by" dropdown (headline scores only)
app$score_cols <- c("sig_prolif","sig_cytokinesis","sig_ccexit","sig_ploidy",
                    "sig_maturation","sig_metabolic")
app$score_cols <- intersect(app$score_cols, names(app$meta))

cat("\n== summary ==\n"); print(app$score_meta, row.names = FALSE)
cat(sprintf("\nExposed in UMAP dropdown: %s\n", paste(app$score_cols, collapse = ", ")))

cat("\nBacking up -> app_data.pre_scores.bak.rds\n")
file.copy("app_data.rds", "app_data.pre_scores.bak.rds", overwrite = TRUE)
cat("Saving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
