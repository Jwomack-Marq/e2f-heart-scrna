#!/usr/bin/env Rscript
# One-time: bundle a broad-gene log-norm matrix for the interactive "Subset & DEGs"
# explorer. To fit the shinyapps.io free tier (~1 GB RAM) we keep broad gene
# coverage (genes detected in >= MIN_CELLS) but DOWNSAMPLE cells (stratified) to
# TARGET_CELLS — nonzero count (hence memory) scales with cells, not genes. The
# explorer gets its own self-consistent cells + metadata (app$deg_expr/meta/genes).
#
#   Rscript build_expr_full.R [path/to/app_data.rds]

suppressWarnings(suppressMessages({ library(Seurat); library(SeuratObject); library(Matrix) }))

args <- commandArgs(trailingOnly = TRUE)
RDS  <- if (length(args)) args[1] else
  "C:/Users/Justi/OneDrive/Documents/GitHub/e2f-heart-scrna/shiny_app/app_data.rds"
PROC <- "C:/Users/Justi/OneDrive - Marquette University/Personal/E2F 7_8/Han_scRNA_2025/our_analysis/processing"
MIN_CELLS    <- 5
TARGET_CELLS <- 8000
stopifnot(file.exists(RDS))
say <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

app  <- readRDS(RDS)
meta <- app$meta

# stratified downsample of the 30k live cells (preserve celltype x genotype x timepoint x Phase)
set.seed(1)
strat <- interaction(meta$celltype, meta$genotype, meta$timepoint,
                     if ("Phase" %in% names(meta)) meta$Phase else "x", drop = TRUE)
frac  <- TARGET_CELLS / nrow(meta)
idx <- unlist(lapply(split(seq_len(nrow(meta)), strat), function(ii)
  if (length(ii) <= 4) ii else sample(ii, max(4, round(length(ii) * frac)))))
idx <- sort(idx)
sub_cells <- meta$cell[idx]
say(sprintf("downsampled cells: %d -> %d", nrow(meta), length(sub_cells)))

say("loading annotated Seurat object ...")
obj <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]), error = function(e) obj[["RNA"]])
m <- GetAssayData(obj, assay = "RNA", layer = "data")
if (is.null(m) || !length(m@x) || max(m@x, na.rm = TRUE) > 100) {
  say("normalizing"); obj <- NormalizeData(obj, verbose = FALSE)
  m <- GetAssayData(obj, assay = "RNA", layer = "data")
}
sub_cells <- sub_cells[sub_cells %in% colnames(m)]
m <- m[, sub_cells, drop = FALSE]
det <- Matrix::rowSums(m > 0)
m <- as(m[det >= MIN_CELLS, , drop = FALSE], "CsparseMatrix")
say(sprintf("deg_expr: %d genes x %d cells | object.size = %.0f MB | nnz = %d",
            nrow(m), ncol(m), as.numeric(object.size(m)) / 1e6, length(m@x)))

# metadata aligned to the matrix columns
dm <- meta[match(colnames(m), meta$cell), , drop = FALSE]
stopifnot(identical(dm$cell, colnames(m)))

app$expr_full <- NULL                       # ensure no stale large object lingers
app$deg_expr  <- m
app$deg_meta  <- dm
app$deg_genes <- rownames(m)

bak <- sub("\\.rds$", ".pre_exprfull.bak.rds", RDS)
if (!file.exists(bak)) file.copy(RDS, bak)
say("saving (xz) ...")
saveRDS(app, RDS, compress = "xz")
cat(sprintf("saved: %s (%.0f MB on disk)\n", RDS, file.info(RDS)$size / 1e6))
cat("=== DONE build_expr_full ===\n")
