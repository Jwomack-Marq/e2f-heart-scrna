# build_subcluster_enrichment.R
# ---------------------------------------------------------------------------
# Precompute per-subcluster enrichment for the cardiomyocyte deep-dive (res 0.2)
# and write it back into app_data.rds as app$enrich$sub.
#
#   Run LOCALLY, in a real R session (NOT the bash sandbox — it segfaults on the
#   103 MB gzip rds):
#     & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" shiny_app/build_subcluster_enrichment.R
#   (run from the repo root, or cd into shiny_app first — paths below are relative
#    to shiny_app/).
#
# Produces three flat data.frames, each keyed by a `subcluster` column so the app
# can display them with the same helpers used for the cell-type enrichment tab:
#   app$enrich$sub$go          KO-vs-WT GO BP (direction = KO_up / KO_down)
#   app$enrich$sub$gsea        KO-vs-WT GSEA (Hallmark + KEGG_LEGACY)
#   app$enrich$sub$identity_go GO BP of each subcluster's up-markers (identity)
#
# Sources (all already inside app_data.rds — no offsite Seurat object needed):
#   KO-vs-WT   <- app$tables$sub_DE[["res0.2"]]  (per-subcluster DE tables)
#   identity   <- app$deg_expr (broad log-norm matrix) subset to CM cells,
#                 labelled by app$cm$meta$SCT_snn_res.0.2
# ---------------------------------------------------------------------------

suppressMessages({
  library(Matrix)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(fgsea)
  library(msigdbr)
  library(presto)
})

setwd_if <- function() if (file.exists("app_data.rds")) invisible(TRUE)
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

RES        <- "0.2"                 # user decision: new features work only at res 0.2
RES_KEY    <- paste0("res", RES)    # "res0.2"
RES_COL    <- paste0("SCT_snn_res.", RES)
GO_PCUT    <- 0.2                   # permissive so small clusters still surface terms
GO_QCUT    <- 0.2
KO_LFC     <- 1                     # |log2FC| threshold for KO-up / KO-down gene lists
KO_PADJ    <- 0.05
ID_MIN_CELLS <- 25                  # skip identity for tiny clusters (drops CM12 = 9 cells)
ID_TOPN    <- 150                   # top markers fed to enrichGO
ID_AUC     <- 0.55
ID_PADJ    <- 0.05

cat("Loading app_data.rds ...\n")
app  <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_subenrich.bak.rds\n")
saveRDS(app, "app_data.pre_subenrich.bak.rds", compress = "gzip")

cmm   <- app$cm$meta
subDE <- app$tables$sub_DE[[RES_KEY]]
EXPR  <- app$deg_expr
CONF  <- app$confound                       # sex/construct confounder genes
stopifnot(!is.null(subDE), !is.null(EXPR), RES_COL %in% names(cmm))

# ---- gene sets for GSEA (match the cell-type build: Hallmark + KEGG_LEGACY, keep prefixes)
cat("Fetching MSigDB gene sets (Hallmark + KEGG_LEGACY, mouse) ...\n")
H  <- msigdbr(species = "Mus musculus", collection = "H")
KG <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
pathways <- split(c(H$gene_symbol, KG$gene_symbol), c(H$gs_name, KG$gs_name))
pathways <- lapply(pathways, unique)
cat(sprintf("  %d gene sets\n", length(pathways)))

# safe enrichGO -> data.frame or NULL
run_go <- function(genes, universe) {
  genes <- unique(genes[nzchar(genes)])
  if (length(genes) < 10) return(NULL)
  eg <- tryCatch(
    suppressWarnings(suppressMessages(
      enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
               universe = unique(universe), pvalueCutoff = GO_PCUT, qvalueCutoff = GO_QCUT,
               minGSSize = 10, maxGSSize = 500, readable = FALSE))),
    error = function(e) NULL)
  if (is.null(eg)) return(NULL)
  df <- as.data.frame(eg)
  if (!nrow(df)) return(NULL)
  df
}

# ---------------------------------------------------------------- KO-vs-WT ----
go_rows <- list(); gsea_rows <- list()
cls <- names(subDE)                                     # CM0 .. CM11
cls <- cls[order(as.integer(sub("CM", "", cls)))]
cat("\n== KO-vs-WT GO + GSEA per subcluster ==\n")
for (cl in cls) {
  d <- subDE[[cl]]
  if (is.null(d) || !nrow(d)) { cat(sprintf("  %-5s (no DE)\n", cl)); next }
  universe <- unique(d$gene)
  nonconf  <- !d$confounder

  up   <- d$gene[d$log2FoldChange >=  KO_LFC & d$padj < KO_PADJ & nonconf]
  down <- d$gene[d$log2FoldChange <= -KO_LFC & d$padj < KO_PADJ & nonconf]
  g_up   <- run_go(up,   universe)
  g_down <- run_go(down, universe)
  n_up_go <- 0L; n_down_go <- 0L
  if (!is.null(g_up))   { g_up$subcluster   <- cl; g_up$direction   <- "KO_up";   go_rows[[paste0(cl,"_up")]]   <- g_up;   n_up_go   <- nrow(g_up) }
  if (!is.null(g_down)) { g_down$subcluster <- cl; g_down$direction <- "KO_down"; go_rows[[paste0(cl,"_down")]] <- g_down; n_down_go <- nrow(g_down) }

  # GSEA: rank all non-confounder tested genes by signed -log10 p
  dd <- d[nonconf & is.finite(d$pvalue) & is.finite(d$log2FoldChange), ]
  stat <- sign(dd$log2FoldChange) * -log10(pmax(dd$pvalue, 1e-300))
  names(stat) <- dd$gene
  stat <- stat[!duplicated(names(stat))]
  stat <- sort(stat[is.finite(stat)], decreasing = TRUE)
  n_gsea <- 0L
  if (length(stat) >= 15) {
    set.seed(1)
    fg <- tryCatch(suppressWarnings(fgsea(pathways, stat, minSize = 10, maxSize = 500)),
                   error = function(e) NULL)
    if (!is.null(fg) && nrow(fg)) {
      fg$leadingEdge <- vapply(fg$leadingEdge, paste, "", collapse = ", ")
      fg$subcluster  <- cl
      gsea_rows[[cl]] <- as.data.frame(fg)
      n_gsea <- nrow(fg)
    }
  }
  cat(sprintf("  %-5s  KO-up genes=%3d GO=%3d | KO-down genes=%3d GO=%3d | GSEA sets=%3d\n",
              cl, length(up), n_up_go, length(down), n_down_go, n_gsea))
}

# ---------------------------------------------------------------- identity ----
cat("\n== Identity GO per subcluster (markers from broad deg_expr) ==\n")
id_rows <- list()
cells <- intersect(cmm$cell, colnames(EXPR))
X   <- EXPR[, cells, drop = FALSE]
lab <- paste0("CM", cmm[[RES_COL]][match(cells, cmm$cell)])
uni <- rownames(X)[Matrix::rowSums(X) > 0]
tabl <- table(lab)
keep_cls <- names(tabl)[tabl >= ID_MIN_CELLS]
keep_cls <- keep_cls[order(as.integer(sub("CM", "", keep_cls)))]
cat(sprintf("  %d CM cells in broad matrix; testing %d clusters (>= %d cells)\n",
            length(cells), length(keep_cls), ID_MIN_CELLS))
wa <- suppressWarnings(presto::wilcoxauc(X, lab))       # one-vs-rest per group, single pass
for (cl in keep_cls) {
  w <- wa[wa$group == cl, ]
  mk <- w[w$auc > ID_AUC & w$padj < ID_PADJ & w$logFC > 0, ]
  mk <- head(mk[order(-mk$auc), ], ID_TOPN)
  g  <- run_go(mk$feature, uni)
  n_go <- 0L
  if (!is.null(g)) { g$subcluster <- cl; id_rows[[cl]] <- g; n_go <- nrow(g) }
  cat(sprintf("  %-5s  markers=%3d GO=%3d\n", cl, nrow(mk), n_go))
}

# ---------------------------------------------------------------- assemble ----
bind <- function(L) { L <- L[!vapply(L, is.null, logical(1))]
                      if (length(L)) do.call(rbind, L) else NULL }
app$enrich$sub <- list(res = RES,
                       go          = bind(go_rows),
                       gsea        = bind(gsea_rows),
                       identity_go = bind(id_rows))

cat("\n== Summary ==\n")
cat(sprintf("  go rows          : %s\n", if (is.null(app$enrich$sub$go)) 0 else nrow(app$enrich$sub$go)))
cat(sprintf("  gsea rows        : %s\n", if (is.null(app$enrich$sub$gsea)) 0 else nrow(app$enrich$sub$gsea)))
cat(sprintf("  identity_go rows : %s\n", if (is.null(app$enrich$sub$identity_go)) 0 else nrow(app$enrich$sub$identity_go)))

cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
