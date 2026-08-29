#!/usr/bin/env Rscript
# How many principal components should go into the neighbour graph and UMAP?
#
# Both production embeddings hard-code dims = 1:30 (combined.R:36-38 for the whole heart,
# cm_subcluster_build.R:39-40 for the cardiomyocytes) with no recorded justification.
# docs/01-upstream.qmd lists "the number of principal components carried into neighbours
# and UMAP" as one of the parameters the book cannot cite. This measures it: same object,
# same SCT/PCA/Harmony, only the dims cut varies -- 10, 30, 50.
#
# WHAT IS HELD FIXED, and why it is the right cut. SCTransform, RunPCA (npcs = 50) and
# RunHarmony run ONCE and are shared by all three arms. Harmony corrects every dimension
# of the PCA reduction, and production then takes the first 30 of the CORRECTED embedding
# -- so "carrying N PCs" means slicing the harmony reduction, not re-running Harmony on a
# shorter PCA. Holding the expensive stages fixed also means the three arms differ by
# exactly one argument, which is the only way the diff is attributable.
#
#   Rscript pcdims_sweep.R --object=cm         # cardiomyocyte compartment (42,416 cells)
#   Rscript pcdims_sweep.R --object=combined   # whole heart (58,917 cells)
#
# DIAGNOSTIC ONLY. Writes new tables/figures; never touches the production objects.
#   Reads  processing/seurat.combined.annotated.rds
#   Writes results/tables/pcdims_<obj>_comparison.csv
#          results/tables/pcdims_<obj>_percell.csv.gz
#          results/figures/pcdims_<obj>_umap.png
#          processing/pcdims_<obj>_embeddings.rds
#
# Runs in e2f-seurat-full:latest (our_analysis/Dockerfile.seurat).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(harmony); library(ggplot2); library(patchwork) }))

# _common.R's 3 GiB brake is below what SCTransform needs on either object here; see the
# note in cm_objectmode_check.R. plan() is sequential, so nothing is ever shipped to a
# worker and the limit only refuses to run. Raised, still bounded, never Inf.
options(future.globals.maxSize = 16 * 1024^3)

OBJ <- local({
  a <- grep("^--object=", commandArgs(TRUE), value = TRUE)
  o <- if (length(a)) sub("^--object=", "", a[1]) else "cm"
  if (!o %in% c("cm", "combined")) stop("--object must be 'cm' or 'combined'")
  o
})
DIMS_SWEEP <- c(10, 30, 50)          # 30 is what both production scripts use
NPCS       <- 50                     # RunPCA default; 50 is the ceiling of the sweep
SEED       <- 42L
# Cluster at each object's own production resolution, so the ARI between dims arms is
# answering "would the clusters people actually look at change?" rather than an arbitrary
# resolution's question. combined.R uses 0.8; cm_subcluster_build.R sweeps and the app
# shows 0.2.
RES_SWEEP  <- if (OBJ == "cm") c(0.1, 0.2, 0.3) else 0.8
RES_COLS   <- paste0("SCT_snn_res.", RES_SWEEP)
KEEP_META  <- if (OBJ == "cm") {
  c("genotype","timepoint","Phase","orig.ident")
} else {
  c("celltype","genotype","timepoint","orig.ident")
}
TAG        <- function(d) paste0("dims", d)

## --- comparison maths (shared with cm_objectmode_check.R) -------------------
# Procrustes: rotate/reflect Y onto X after centring and scaling both to unit Frobenius
# norm. Reflection is allowed deliberately -- a mirrored UMAP is the same embedding.
procrustes_rmsd <- function(X, Y) {
  X <- scale(X, TRUE, FALSE); Y <- scale(Y, TRUE, FALSE)
  X <- X / sqrt(sum(X^2));   Y <- Y / sqrt(sum(Y^2))
  s <- svd(crossprod(Y, X)); R <- s$u %*% t(s$v)
  sqrt(mean(rowSums((X - Y %*% R)^2)))
}
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  cc <- function(x) sum(choose(x, 2))
  idx <- cc(tab); ea <- cc(rowSums(tab)); eb <- cc(colSums(tab))
  expct <- ea * eb / choose(n, 2); mx <- (ea + eb) / 2
  if (isTRUE(all.equal(mx, expct))) return(1)
  (idx - expct) / (mx - expct)
}
greedy_agree <- function(a, b) {
  tab <- table(a, b); tot <- sum(tab); hit <- 0
  while (nrow(tab) > 0 && ncol(tab) > 0 && max(tab) > 0) {
    w <- which(tab == max(tab), arr.ind = TRUE)[1, ]
    hit <- hit + tab[w[1], w[2]]; tab <- tab[-w[1], -w[2], drop = FALSE]
  }
  hit / tot
}

## --- load, and build the shared SCT/PCA/Harmony once ------------------------
inp <- file.path(PROC, "seurat.combined.annotated.rds")
cat(sprintf("[%s] reading %s (mtime %s)\n", OBJ, basename(inp),
            format(file.mtime(inp), "%Y-%m-%d %H:%M:%S")))
obj <- readRDS(inp)
stopifnot("celltype" %in% colnames(obj@meta.data))
obj$genotype <- genotype_of(obj$orig.ident)
if (OBJ == "cm") {
  obj <- obj[, obj$celltype == "Cardiomyocyte"]
  cat(sprintf("[%s] cardiomyocytes: %d cells x %d genes\n", OBJ, ncol(obj), nrow(obj)))
} else {
  cat(sprintf("[%s] whole heart: %d cells x %d genes\n", OBJ, ncol(obj), nrow(obj)))
}
prov <- list(object = OBJ, input_mtime = format(file.mtime(inp)),
             cells = ncol(obj), genes = nrow(obj))
gc(verbose = FALSE)

DefaultAssay(obj) <- "RNA"                     # parent's DefaultAssay is "SCT"; see
                                               # cm_objectmode_check.R for why this matters
cat(sprintf("[%s] SCTransform ...\n", OBJ)); set.seed(SEED)
obj <- SCTransform(obj, method = "glmGamPoi", conserve.memory = TRUE,
                   verbose = FALSE, seed.use = SEED)
cat(sprintf("[%s] RunPCA (npcs = %d) ...\n", OBJ, NPCS)); set.seed(SEED)
obj <- RunPCA(obj, npcs = NPCS, verbose = FALSE, seed.use = SEED)
cat(sprintf("[%s] RunHarmony ...\n", OBJ)); set.seed(SEED)
obj <- RunHarmony(obj, group.by.vars = "orig.ident")
stopifnot(ncol(Embeddings(obj, "harmony")) >= max(DIMS_SWEEP))

# How much variance each cut carries, as a share of what the COMPUTED PCs hold.
# The denominator is deliberately named: Stdev() returns only the npcs we asked for, so
# this is a fraction of the top-50 variance, NOT of the total variance in the data --
# which is why the 50-PC row is 100% by construction and carries no information on its
# own. What is informative is the gap between cuts: whatever 30 leaves on the table is
# the most PCs 31..50 could possibly contribute.
sdev <- Stdev(obj, reduction = "pca")
varpct <- cumsum(sdev^2) / sum(sdev^2) * 100
cat(sprintf("[%s] %% of the variance in the top %d PCs: 10 PCs %.1f%%, 30 PCs %.1f%%, 50 PCs %.1f%%\n",
            OBJ, NPCS, varpct[10], varpct[30], varpct[NPCS]))

## --- the sweep: only `dims` varies ------------------------------------------
res <- list()
for (d in DIMS_SWEEP) {
  cat(sprintf("[%s] --- dims = 1:%d ---\n", OBJ, d))
  set.seed(SEED)
  o <- FindNeighbors(obj, reduction = "harmony", dims = 1:d, verbose = FALSE)
  set.seed(SEED)
  o <- FindClusters(o, resolution = RES_SWEEP, verbose = FALSE, random.seed = SEED)
  set.seed(SEED)
  o <- RunUMAP(o, reduction = "harmony", dims = 1:d, reduction.name = "umap",
               verbose = FALSE, seed.use = SEED)
  res[[TAG(d)]] <- list(
    dims     = d,
    cells    = colnames(o),
    umap     = Embeddings(o, "umap")[, 1:2, drop = FALSE],
    clusters = o@meta.data[, RES_COLS, drop = FALSE],
    nclust   = vapply(RES_COLS, function(rc) length(unique(o[[rc]][, 1])), integer(1)))
  for (rc in RES_COLS)
    cat(sprintf("     %-18s %d clusters\n", rc, res[[TAG(d)]]$nclust[[rc]]))
  # Carry this arm's labels and embedding back onto the SHARED object, under names that
  # say which cut produced them. One object ends up holding every arm, which is the point:
  # SCT/PCA/Harmony are identical across arms by construction, so three saved objects
  # would be three copies of the same 2.3 GB matrices differing in a metadata column.
  for (rc in RES_COLS) {
    vcol <- sprintf("dims%d_res%s", d, sub("^SCT_snn_res\\.", "", rc))
    obj[[vcol]] <- as.character(o@meta.data[[rc]])
  }
  obj[[sprintf("umap.dims%d", d)]] <- SeuratObject::CreateDimReducObject(
    embeddings = Embeddings(o, "umap")[, 1:2, drop = FALSE],
    key = sprintf("umapd%d_", d), assay = DefaultAssay(o))
  rm(o); gc(verbose = FALSE)
}
## --- the variant object + the registry that names its labellings -------------
# Only for the CM compartment: the per-subcluster downstream (markers, pseudobulk DE,
# enrichment) is CM-specific, and nothing consumes a whole-heart variant object, so
# writing 3 GB for it would be dead weight.
VARIANTS <- do.call(rbind, lapply(DIMS_SWEEP, function(d) do.call(rbind, lapply(RES_SWEEP,
  function(r) data.frame(
    variant_id    = sprintf("%s_dims%d_res%s", OBJ, d, r),
    object        = OBJ, dims = d, resolution = r,
    cluster_col   = sprintf("dims%d_res%s", d, r),
    umap_reduction= sprintf("umap.dims%d", d),
    n_clusters    = res[[TAG(d)]]$nclust[[paste0("SCT_snn_res.", r)]],
    harmony_var   = "orig.ident",
    source_script = "our_analysis/04_integrate_annotate/pcdims_sweep.R",
    created       = format(Sys.Date()),
    # dims 30 / res 0.2 is what cm_subcluster_build.R and the app ship today. Flagged so
    # the app can default to it and mark everything else as not-the-published-labelling.
    is_production = (OBJ == "cm" && d == 30 && r == 0.2),
    stringsAsFactors = FALSE)))))

if (OBJ == "cm") {
  vpath <- file.path(PROC, "cm_variants.rds")
  saveRDS(obj, vpath)
  cat(sprintf("\n[%s] variant object -> %s (%d cells, %d variant columns)\n",
              OBJ, basename(vpath), ncol(obj), nrow(VARIANTS)))
}
# Upsert rather than overwrite: the registry is meant to accumulate as new labellings are
# added (a resolution sweep, a lane-harmonised run), and running one object must not erase
# another's rows.
regf <- file.path(OUTTAB, "clustering_registry.csv")
if (file.exists(regf)) {
  old_reg <- read.csv(regf, stringsAsFactors = FALSE)
  old_reg <- old_reg[!old_reg$variant_id %in% VARIANTS$variant_id, , drop = FALSE]
  VARIANTS <- rbind(old_reg[, intersect(names(old_reg), names(VARIANTS)), drop = FALSE], VARIANTS)
}
VARIANTS <- VARIANTS[order(VARIANTS$object, VARIANTS$dims, VARIANTS$resolution), ]
write.csv(VARIANTS, regf, row.names = FALSE)
cat(sprintf("[%s] registry -> %s (%d variants)\n", OBJ, basename(regf), nrow(VARIANTS)))

meta_keep <- obj@meta.data[, intersect(KEEP_META, colnames(obj@meta.data)), drop = FALSE]
rm(obj); gc(verbose = FALSE)

cells <- res[[1]]$cells
stopifnot(all(vapply(res, function(r) identical(r$cells, cells), logical(1))))

## --- comparison table --------------------------------------------------------
row <- function(pair, metric, value, note = "") data.frame(
  pair = pair, metric = metric, value = value, note = note, stringsAsFactors = FALSE)
out <- list(
  row("(all)", "cells",              length(cells)),
  row("(all)", "npcs_computed",      NPCS),
  row("(all)", "varpct_10",   varpct[10],   "% of the variance in the top 50 PCs, not of total variance"),
  row("(all)", "varpct_30",   varpct[30],   "% of the variance in the top 50 PCs, not of total variance"),
  row("(all)", "varpct_50",   varpct[NPCS], "100% by construction: the denominator IS the top 50 PCs"))
for (d in DIMS_SWEEP) for (rc in RES_COLS)
  out <- c(out, list(row(TAG(d), paste0("nclust_", rc), res[[TAG(d)]]$nclust[[rc]])))
prs <- list(c(10, 30), c(30, 50), c(10, 50))
for (p in prs) {
  a <- res[[TAG(p[1])]]; b <- res[[TAG(p[2])]]
  tag <- sprintf("dims %d vs %d", p[1], p[2])
  out <- c(out, list(row(tag, "umap_procrustes_rmsd", procrustes_rmsd(a$umap, b$umap),
                         "0 = same shape; reflection allowed")))
  for (rc in RES_COLS) {
    ca <- as.character(a$clusters[[rc]]); cb <- as.character(b$clusters[[rc]])
    out <- c(out, list(
      row(tag, paste0("ARI_", rc),          ari(ca, cb)),
      row(tag, paste0("agree_", rc),        greedy_agree(ca, cb), "greedy max-overlap match"),
      row(tag, paste0("ncells_moved_", rc), round((1 - greedy_agree(ca, cb)) * length(ca)),
          "greedy max-overlap match")))
  }
}
cmp <- do.call(rbind, out)
write.csv(cmp, file.path(OUTTAB, sprintf("pcdims_%s_comparison.csv", OBJ)), row.names = FALSE)

## --- per-cell coordinates (feeds the app tab) --------------------------------
percell <- do.call(rbind, lapply(DIMS_SWEEP, function(d) {
  r <- res[[TAG(d)]]
  cbind(data.frame(cell = r$cells, dims = d, UMAP1 = r$umap[, 1], UMAP2 = r$umap[, 2],
                   stringsAsFactors = FALSE),
        r$clusters, meta_keep[r$cells, , drop = FALSE]) }))
gz <- gzfile(file.path(OUTTAB, sprintf("pcdims_%s_percell.csv.gz", OBJ)), "w")
write.csv(percell, gz, row.names = FALSE); close(gz)
saveRDS(c(res, list(provenance = prov, comparison = cmp, varpct = varpct)),
        file.path(PROC, sprintf("pcdims_%s_embeddings.rds", OBJ)))

## --- figure ------------------------------------------------------------------
colvar <- if (OBJ == "cm") RES_COLS[[if (length(RES_COLS) >= 2) 2 else 1]] else "celltype"
ps <- lapply(DIMS_SWEEP, function(d) {
  r <- res[[TAG(d)]]
  dd <- data.frame(r$umap); names(dd) <- c("UMAP1", "UMAP2")
  dd$grp <- if (OBJ == "cm") factor(paste0("CM", r$clusters[[colvar]]))
            else             factor(meta_keep[r$cells, "celltype"])
  ggplot(dd, aes(UMAP1, UMAP2, colour = grp)) +
    geom_point(size = .2, alpha = .6, shape = 16) + theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    ggtitle(sprintf("dims 1:%d  (%.1f%% of the top-%d PC variance)", d, varpct[d], NPCS))
})
png(file.path(OUTFIG, sprintf("pcdims_%s_umap.png", OBJ)), 1650, 560)
print(wrap_plots(ps, nrow = 1)); dev.off()

## --- console verdict ---------------------------------------------------------
cat(sprintf("\n=== VERDICT (%s) ===\n", OBJ))
for (p in prs) {
  tag <- sprintf("dims %d vs %d", p[1], p[2])
  g <- function(m) cmp$value[cmp$pair == tag & cmp$metric == m][1]
  cat(sprintf("  %-16s procrustes %-10s ARI(%s) %-9s cells moved %s\n", tag,
              signif(g("umap_procrustes_rmsd"), 4), sub("SCT_snn_res.", "res ", RES_COLS[1]),
              signif(g(paste0("ARI_", RES_COLS[1])), 4),
              g(paste0("ncells_moved_", RES_COLS[1]))))
}
cat(sprintf("\nwrote: pcdims_%s_{comparison.csv,percell.csv.gz}, pcdims_%s_umap.png\n", OBJ, OBJ))
cat(sprintf("=== DONE pcdims_sweep (%s) ===\n", OBJ))
