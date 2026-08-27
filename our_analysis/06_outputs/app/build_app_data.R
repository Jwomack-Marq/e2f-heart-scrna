#!/usr/bin/env Rscript
# Build the SLIM data object (app_data.rds) that powers the browser cell-explorer.
# The full Seurat objects are multi-GB and cannot load in a browser (shinylive/webR),
# so we extract only what the app needs:
#   * meta   : per-cell data frame (UMAP coords + metadata)         [downsampled]
#   * expr   : curated sparse log-norm gene x cell matrix (RNA)      [downsampled]
#   * cm     : CM-subcluster UMAP coords + per-resolution labels     [downsampled]
#   * tables : result CSVs + per-cell-type DE + per-subcluster DE + summaries
#   * heat   : precomputed subcluster marker-expression heatmaps (z-scored)
#   * genes  : searchable gene vector (rownames of expr)
#
# DE tables / heatmaps are precomputed from the FULL data, so downsampling the live
# cells only thins the UMAP/violin views (still ~30k cells) — DE is unaffected.
# webR reads this with readRDS() and plots with ggplot2 + Matrix (no Seurat).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))   # OUR_ROOT, PROC, OUTTAB, OUT, gene sets, Seurat
suppressWarnings(suppressMessages({ library(Matrix); library(base64enc) }))

APP_DIR <- file.path(OUT, "app"); if (!dir.exists(APP_DIR)) dir.create(APP_DIR, recursive = TRUE)

# ---- tunables --------------------------------------------------------------
TARGET_GENES      <- 2600    # curated panel + HVG fill, capped here
DOWNSAMPLE_N      <- 30000   # live cells kept (stratified); NULL = all
DE_MIN_BASEMEAN   <- 3       # drop ultra-low-expression genes from stored DE tables
MARKERS_PER_SUB   <- 6       # top markers per subcluster for the identity heatmap
RES_KEEP          <- c("0.1", "0.2", "0.3")   # CM subcluster resolutions to surface (those present + precomputed)
ROUND_EXPR_DIGITS <- 2

CM_SUBTYPES <- list(
  Ventricular = c("Myl2","Myh7"), Atrial = c("Myl7","Sln","Nppa"),
  Trabecular  = c("Bmp10","Nppa","Hey2"), Compact = c("Hey2","Irx3","Tbx20"),
  Cycling     = c("Mki67","Top2a","Ccnb1","Aurkb","Cdca8"))
want_cols <- c("orig.ident","genotype","timepoint","celltype","seurat_clusters",
               "Phase","cycling","cm_subtype","pseudotime","S.Score","G2M.Score")

# ---- helpers ---------------------------------------------------------------
get_norm <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]), error = function(e) obj[["RNA"]])
  m <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"), error = function(e) NULL)
  if (is.null(m) || !nrow(m) || max(m@x, na.rm = TRUE) > 100) {
    obj <- NormalizeData(obj, verbose = FALSE); m <- GetAssayData(obj, assay = "RNA", layer = "data")
  }
  m
}
umap_xy <- function(obj) {
  rn <- Reductions(obj); pick <- rn[match(TRUE, rn %in% c("umap","umap.harmony","UMAP"))]
  if (is.na(pick)) pick <- grep("umap", rn, ignore.case = TRUE, value = TRUE)[1]
  e <- Embeddings(obj, pick)[, 1:2]; colnames(e) <- c("UMAP1","UMAP2"); as.data.frame(e)
}
split_genes <- function(x) unique(trimws(unlist(strsplit(x[!is.na(x) & nzchar(x)], ","))))
read_tab <- function(name) { p <- file.path(OUTTAB, name); if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL }
trim_de <- function(d) {                         # slim a full DE CSV for storage/volcano
  if (is.null(d)) return(NULL)
  keep <- c("gene","log2FoldChange","baseMean","pvalue","padj","confounder")
  d <- d[, intersect(keep, names(d)), drop = FALSE]
  d <- d[is.finite(d$log2FoldChange) & is.finite(d$baseMean) & d$baseMean >= DE_MIN_BASEMEAN, ]
  d$log2FoldChange <- round(d$log2FoldChange, 3); d$baseMean <- round(d$baseMean, 1)
  if ("pvalue" %in% names(d)) d$pvalue <- signif(d$pvalue, 3)
  if ("padj"   %in% names(d)) d$padj   <- signif(d$padj, 3)
  if (!"confounder" %in% names(d)) d$confounder <- d$gene %in%
      c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")
  rownames(d) <- NULL; d
}

# ---- 1. annotated 4-group object (primary) ---------------------------------
cat("Loading annotated object ...\n")
obj <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
cat(sprintf("  %d cells\n", ncol(obj)))
meta_full <- obj@meta.data[, intersect(want_cols, colnames(obj@meta.data)), drop = FALSE]
meta_full <- cbind(umap_xy(obj), meta_full); meta_full$cell <- rownames(meta_full)
expr_all <- get_norm(obj)

# curated gene panel (for the live gene-colouring views)
de_top <- character(0)
sm <- read_tab("percelltype_KOvsWT_summary.csv")
if (!is.null(sm)) de_top <- c(split_genes(sm$top_KO_up), split_genes(sm$top_KO_down))
for (f in c("P0.cardiac.descriptive.DE.csv","P7.cardiac.descriptive.DE.csv")) {
  d <- read_tab(f); if (!is.null(d)) { d <- d[is.finite(d$log2FoldChange), ]
    de_top <- c(de_top, head(d$gene[order(-abs(d$log2FoldChange))], 60)) }
}
hvg <- tryCatch(head(VariableFeatures(obj, assay = "SCT"), 2000), error = function(e) character(0))
if (!length(hvg)) hvg <- tryCatch(head(VariableFeatures(obj), 2000), error = function(e) character(0))
curated <- intersect(unique(c(unlist(CELLTYPE_MARKERS), E2F_TARGETS, unlist(CM_SUBTYPES),
                    unlist(cc_lists()), CM_MATURE, CM_IMMATURE, de_top,
                    "E2f7","E2f8","Gabbr2","Tcf4","Adamts9","Ralyl")), rownames(expr_all))
panel <- head(c(curated, setdiff(intersect(hvg, rownames(expr_all)), curated)), TARGET_GENES)
cat(sprintf("  gene panel: %d genes\n", length(panel)))

# ---- 2. CM-subclustered object (labels + 2nd embedding) --------------------
cat("Loading CM-subclustered object ...\n")
cmobj  <- readRDS(file.path(PROC, "seurat.cm.subclustered.rds"))
res_cols <- paste0("SCT_snn_res.", RES_KEEP)
cm_cols  <- intersect(c("genotype","timepoint", res_cols, "Phase","cycling"), colnames(cmobj@meta.data))
cm_full  <- cbind(umap_xy(cmobj), cmobj@meta.data[, cm_cols, drop = FALSE]); cm_full$cell <- rownames(cm_full)
rm(cmobj); invisible(gc())
RES_KEEP <- RES_KEEP[paste0("SCT_snn_res.", RES_KEEP) %in% names(cm_full)]   # keep only resolutions present on the object
cat(sprintf("  resolutions present: %s\n", paste(RES_KEEP, collapse = ", ")))

# ---- 3. precompute subcluster marker-expression heatmaps (FULL data) -------
cat("Precomputing subcluster marker heatmaps ...\n")
heat <- list()
for (r in RES_KEEP) {
  col <- paste0("SCT_snn_res.", r); if (!col %in% names(cm_full)) next
  tm  <- read_tab(sprintf("cm_subcluster_top_markers_res%s.csv", r))
  if (is.null(tm)) tm <- read_tab(sprintf("cm_subcluster_markers_res%s.csv", r))
  if (is.null(tm)) next
  tm <- tm[order(tm$cluster, -tm$avg_log2FC), ]
  topg <- do.call(rbind, lapply(split(tm, tm$cluster), head, MARKERS_PER_SUB))
  genes_h <- intersect(unique(topg$gene), rownames(expr_all))
  cells_h <- intersect(cm_full$cell, colnames(expr_all))
  lab <- paste0("CM", as.character(cm_full[[col]])); names(lab) <- cm_full$cell; lab <- lab[cells_h]
  m   <- as.matrix(expr_all[genes_h, cells_h, drop = FALSE])
  subs <- unique(lab); subs <- subs[order(as.integer(sub("CM", "", subs)))]
  mat <- vapply(subs, function(s) rowMeans(m[, lab == s, drop = FALSE]), numeric(length(genes_h)))
  z <- t(scale(t(mat))); z[!is.finite(z)] <- 0; z[z > 2.5] <- 2.5; z[z < -2.5] <- -2.5
  gene_ord <- genes_h[order(max.col(z, "first"), -apply(z, 1, max))]    # group genes by their peak subcluster
  long <- expand.grid(gene = gene_ord, cluster = colnames(z), stringsAsFactors = FALSE)
  long$z <- z[cbind(match(long$gene, rownames(z)), match(long$cluster, colnames(z)))]
  heat[[paste0("res", r)]] <- list(long = long, genes = gene_ord, clusters = colnames(z))
}

# ---- 4. downsample live cells ----------------------------------------------
keep_cells <- meta_full$cell
if (!is.null(DOWNSAMPLE_N) && DOWNSAMPLE_N < nrow(meta_full)) {
  set.seed(1); frac <- DOWNSAMPLE_N / nrow(meta_full)
  grp <- interaction(meta_full$genotype, meta_full$timepoint, meta_full$celltype, drop = TRUE)
  keep_cells <- unlist(lapply(split(meta_full$cell, grp), function(ix) {
    k <- max(min(length(ix), 40L), round(length(ix) * frac)); if (length(ix) > k) sample(ix, k) else ix }))
  cat(sprintf("  downsampled live cells: %d -> %d\n", nrow(meta_full), length(keep_cells)))
}
meta <- meta_full[match(keep_cells, meta_full$cell), ]
expr <- expr_all[panel, keep_cells, drop = FALSE]
expr <- as(expr, "CsparseMatrix"); if (!is.null(ROUND_EXPR_DIGITS)) expr@x <- round(expr@x, ROUND_EXPR_DIGITS)
rm(expr_all, obj); invisible(gc())
cm_meta <- cm_full[cm_full$cell %in% keep_cells, ]

# ---- 5. DE tables: per cell type + per subcluster + summaries ---------------
cat("Reading DE tables ...\n")
celltypes <- c("Cardiomyocyte","Endothelial","Fibroblast","Immune_Myeloid","Mural_Pericyte","RBC")
ct_DE <- list()
for (tp in c("P0","P7")) for (ct in celltypes) {
  d <- trim_de(read_tab(sprintf("percelltype_%s_%s_KOvsWT.descriptive.DE.csv", tp, ct)))
  if (!is.null(d) && nrow(d)) ct_DE[[paste(tp, ct, sep = "_")]] <- d
}
sub_DE <- list(); sub_summary <- list(); sub_subtype <- list()
for (r in RES_KEEP) {
  files <- list.files(OUTTAB, pattern = sprintf("^cm_subcluster_res%s_KOvsWT_CM[0-9]+\\.descriptive\\.DE\\.csv$", r), full.names = FALSE)
  lst <- list()
  for (f in files) { d <- trim_de(read_tab(f)); sub <- sub("\\.descriptive.*", "", sub(".*KOvsWT_", "", f))
    if (!is.null(d) && nrow(d)) lst[[sub]] <- d }
  sub_DE[[paste0("res", r)]]      <- lst
  sub_summary[[paste0("res", r)]] <- read_tab(sprintf("cm_subcluster_res%s_KOvsWT_summary.csv", r))
  sub_subtype[[paste0("res", r)]] <- read_tab(sprintf("cm_subcluster_subtype_map_res%s.csv", r))
}

tables <- Filter(Negate(is.null), list(
  P0_cardiac_DE = trim_de(read_tab("P0.cardiac.descriptive.DE.csv")),
  P7_cardiac_DE = trim_de(read_tab("P7.cardiac.descriptive.DE.csv")),
  composition   = read_tab("composition.csv"),
  cellcycle_fraction = read_tab("cellcycle_fraction_by_celltype.csv"),
  pseudotime    = read_tab("pseudotime_KOvsWT.csv"),
  e2f_regulon   = read_tab("e2f_regulon_activity.csv"),
  doublet       = read_tab("doublet_comparison.csv"),
  ct_DE = ct_DE, sub_DE = sub_DE, sub_summary = sub_summary, sub_subtype = sub_subtype,
  res_keep = RES_KEEP))

# QC / normalization figures embedded as base64 data-URIs (so they ship in-browser)
fig_uri <- function(name) { p <- file.path(OUTFIG, name)
  if (file.exists(p)) base64enc::dataURI(file = p, mime = "image/png") else NA_character_ }
figs <- list(filtering = fig_uri("norm_filtering_summary.png"), qc_violins = fig_uri("norm_qc_violins.png"),
             doublet = fig_uri("QC_doublet_comparison.png"), hvg = fig_uri("norm_hvg.png"),
             harmony = fig_uri("combined_harmony_before_after.png"))

# ---- 6. assemble + save ----------------------------------------------------
app <- list(meta = meta, expr = expr, genes = sort(rownames(expr)),
            cm = list(meta = cm_meta, res = RES_KEEP), heat = heat, tables = tables, figs = figs,
            confound = c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor"),
            built = as.character(Sys.Date()))
out <- file.path(APP_DIR, "app_data.rds"); saveRDS(app, out, compress = "gzip")
sz <- file.info(out)$size / 1024^2
cat(sprintf("Wrote %s  (%.1f MB)\n", out, sz))
cat(sprintf("  live: %d cells x %d genes | CM: %d cells | ct_DE: %d | sub_DE res: %s | heat res: %s\n",
            nrow(meta), length(app$genes), nrow(cm_meta), length(ct_DE),
            paste(names(sub_DE), collapse = ","), paste(names(heat), collapse = ",")))
cat("=== DONE build_app_data ===\n")
