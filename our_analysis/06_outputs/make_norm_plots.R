#!/usr/bin/env Rscript
# Generate normalization / preprocessing figures for the interactive report's
# "Normalization & preprocessing" section. Faithful to 00_DOCS/NORMALIZATION.md:
#   QC filter (mito + nFeature)  ->  doublet removal  ->  SCTransform  ->  Harmony
# Writes PNGs into results/figures/ (consumed by interactive/report.Rmd via fig()).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))     # paths + Seurat + ggplot2
suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork) }))
ggsave2 <- function(name, p, w = 9, h = 4.6) ggplot2::ggsave(file.path(OUTFIG, name), p, width = w, height = h, dpi = 120, bg = "white")

# ---- 1. QC filtering summary (from qc_comparison.csv; no object needed) -----
qc <- read.csv(file.path(OUT, "qc_comparison.csv"), check.names = FALSE)
rm_long <- do.call(rbind, lapply(c("removed_doublet","removed_highmito","removed_uppercap"), function(col)
  data.frame(sample = qc$sample, reason = sub("removed_", "", col), n = qc[[col]])))
rm_long$reason <- factor(rm_long$reason, levels = c("doublet","highmito","uppercap"),
                         labels = c("doublet (scDblFinder)", "high mito (>20%)", "upper nFeature cap"))
p_filt <- ggplot(rm_long, aes(sample, n, fill = reason)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("doublet (scDblFinder)" = "#c62828", "high mito (>20%)" = "#ef6c00", "upper nFeature cap" = "#1565c0")) +
  labs(title = "Cells removed at QC, by reason and sample/lane", x = NULL, y = "cells removed", fill = NULL) +
  theme_minimal(base_size = 13)
ggsave2("norm_filtering_summary.png", p_filt, h = 5)

# ---- load object for distribution + feature-selection plots -----------------
cat("Loading annotated object ...\n")
obj <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]), error = function(e) obj[["RNA"]])
DefaultAssay(obj) <- "RNA"
if (!"percent.mt" %in% colnames(obj@meta.data))
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# ---- 2. QC distributions (post-filter) with threshold guides ----------------
qcdf <- obj@meta.data[, c("orig.ident","nFeature_RNA","nCount_RNA","percent.mt")]
qcdf$nCount_RNA <- qcdf$nCount_RNA
vln <- function(y, ylab, hline = NULL, log = FALSE) {
  p <- ggplot(qcdf, aes(orig.ident, .data[[y]], fill = orig.ident)) +
    geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2) +
    guides(fill = "none") + theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = ylab)
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "grey30")
  if (log) p <- p + scale_y_log10()
  p
}
p_qc <- vln("nFeature_RNA", "genes / cell", 1500) + vln("nCount_RNA", "UMIs / cell (log10)", log = TRUE) +
        vln("percent.mt", "% mitochondrial", 20) +
        patchwork::plot_annotation(title = "Per-cell QC distributions after filtering (dashed = thresholds: genes ≥1500, mito ≤20%)")
ggsave2("norm_qc_violins.png", p_qc, w = 11, h = 4.4)

# ---- 3. SCTransform / feature selection: mean-variance + HVGs ----------------
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top <- head(VariableFeatures(obj), 10)
p_hvg <- tryCatch({
  vp <- VariableFeaturePlot(obj) + theme_minimal(base_size = 12) +
    labs(title = "Feature selection: variance vs mean expression (red = 2,000 highly-variable genes)")
  LabelPoints(plot = vp, points = top, repel = TRUE, xnudge = 0, ynudge = 0)
}, error = function(e) {
  hv <- HVFInfo(obj); hv$variable <- rownames(hv) %in% VariableFeatures(obj)
  ggplot(hv, aes(mean, variance.standardized, color = variable)) + geom_point(size = .5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#c62828")) + scale_x_log10() +
    theme_minimal(base_size = 12) + labs(title = "Feature selection: standardized variance vs mean", color = "HVG")
})
ggsave2("norm_hvg.png", p_hvg, w = 8.5, h = 5)

cat("=== DONE make_norm_plots ===\n")
for (f in c("norm_filtering_summary.png","norm_qc_violins.png","norm_hvg.png"))
  cat(" wrote", file.path(OUTFIG, f), file.exists(file.path(OUTFIG, f)), "\n")
