# build_refmap.R
# ---------------------------------------------------------------------------
# Reference-marker annotation check, written back into app_data.rds as app$refmap.
#
# For each cell, score it against a set of published-style developmental mouse-heart
# cell-type marker panels (module scores), assign the argmax lineage as a
# "predicted" label, and cross-tabulate against the existing app$meta$celltype call.
# The diagonal of that confusion table should dominate; off-diagonal mass flags
# populations that may be under/mis-annotated in this pilot.
#
# This is a marker-signature CONCORDANCE check, NOT probabilistic Seurat anchor-based
# label transfer (that would need the source query object + a full reference object).
#
# Run LOCALLY, in a real R session (NOT the bash sandbox):
#     & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" shiny_app/build_refmap.R
#
# Output app$refmap = list(
#   pred      : data.frame(cell, celltype, predicted, margin)   # margin = top - 2nd score
#   confusion : data.frame(celltype, predicted, n, prop)        # prop row-normalised
#   panels    : the marker panels used (after intersecting with available genes)
#   coverage  : per-panel gene coverage)
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))

if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

# ---- reference cell-type marker panels (mouse, developmental heart) --------
REF <- list(
  Cardiomyocyte = c("Tnnt2","Myh6","Actc1","Ttn","Tnni3","Nppa","Myl7","Myl2"),
  Fibroblast    = c("Col1a1","Col1a2","Dcn","Pdgfra","Gsn","Postn","Col3a1"),
  Endothelial   = c("Pecam1","Cdh5","Kdr","Fabp4","Egfl7","Emcn","Cldn5"),
  Immune        = c("Ptprc","Cd68","Lyz2","C1qa","Csf1r","Fcgr1","Cd52"),
  Mural         = c("Rgs5","Pdgfrb","Myh11","Acta2","Tagln","Notch3","Kcnj8"),
  Epicardial    = c("Wt1","Tbx18","Upk3b","Msln","Krt19"),
  Neural        = c("Plp1","Kcna1","Sox10","Cdh19","Kcna6"),
  Erythroid     = c("Hba-a1","Hbb-bs","Alas2","Gypa"),
  LymphaticEC   = c("Lyve1","Prox1","Mmrn1","Flt4"))

cat("Loading app_data.rds ...\n")
app  <- readRDS("app_data.rds")
EXPRc <- app$expr
CONF  <- app$confound
stopifnot(!is.null(EXPRc), !is.null(app$meta), "celltype" %in% names(app$meta))
REF <- lapply(REF, function(g) setdiff(g, CONF))

# ---- AddModuleScore-equivalent (mirrors build_signature_scores.R) ----------
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
    n <- min(ctrl, length(pool)); if (n > 0) ctrl_genes <- c(ctrl_genes, sample(pool, n))
  }
  feat_mean <- Matrix::colMeans(mat[features, , drop = FALSE])
  sc <- if (length(ctrl_genes)) feat_mean - Matrix::colMeans(mat[ctrl_genes, , drop = FALSE]) else feat_mean
  setNames(as.numeric(sc), colnames(mat))
}

cat("\n== panel coverage (curated app$expr) ==\n")
cov <- data.frame(panel = names(REF),
                  n_set = vapply(REF, length, 0L),
                  n_used = vapply(REF, function(g) length(intersect(g, rownames(EXPRc))), 0L),
                  stringsAsFactors = FALSE)
print(cov, row.names = FALSE)
REF <- REF[cov$n_used >= 2]                       # drop panels too sparse to score
cat(sprintf("keeping %d panels with >= 2 genes\n", length(REF)))
stopifnot("need >= 2 usable reference panels (curated panel lacks these markers?)" = length(REF) >= 2)

cat("\n== scoring cells against panels ==\n")
S <- vapply(REF, function(g) { v <- module_score(EXPRc, g); if (is.null(v)) rep(NA_real_, ncol(EXPRc)) else v },
            numeric(ncol(EXPRc)))
rownames(S) <- colnames(EXPRc)                    # cells x panels
pred    <- colnames(S)[max.col(S, ties.method = "first")]
sorted  <- t(apply(S, 1, sort, decreasing = TRUE))
margin  <- sorted[, 1] - sorted[, 2]              # confidence: top - runner-up

pdf <- data.frame(cell = colnames(EXPRc), predicted = pred, margin = margin, stringsAsFactors = FALSE)
pdf$celltype <- as.character(app$meta$celltype[match(pdf$cell, app$meta$cell)])
pdf <- pdf[!is.na(pdf$celltype), ]

tab <- as.data.frame(table(celltype = pdf$celltype, predicted = pdf$predicted), stringsAsFactors = FALSE)
names(tab)[names(tab) == "Freq"] <- "n"
tot <- tapply(tab$n, tab$celltype, sum)
tab$prop <- tab$n / as.numeric(tot[tab$celltype])

app$refmap <- list(pred = pdf, confusion = tab, panels = REF, coverage = cov)

cat("\n== confusion (row-normalised, predicted argmax vs existing celltype) ==\n")
wide <- reshape(tab[c("celltype","predicted","prop")], idvar = "celltype",
                timevar = "predicted", direction = "wide")
print(wide, row.names = FALSE, digits = 2)

cat("\nBacking up -> app_data.pre_refmap.bak.rds\n")
file.copy("app_data.rds", "app_data.pre_refmap.bak.rds", overwrite = TRUE)
cat("Saving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
