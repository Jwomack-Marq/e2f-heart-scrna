#!/usr/bin/env Rscript
# Does the CM subset need to be its OWN object? -- a measurement, not an opinion.
#
# The concern raised in review: `cm <- comb[, comb$celltype == "Cardiomyocyte"]` never
# makes the subset a standalone object, so Seurat might "still be using" the unfiltered
# data. cm_subcluster_build.R already re-runs SCTransform -> PCA -> Harmony -> UMAP on the
# subset, so no whole-heart HVG set, scaling or PCA can reach the coordinates. This script
# tests what is left, by building the SAME cells three ways and diffing the results:
#
#   A filter    comb[, celltype == "Cardiomyocyte"]          -- what production does
#   B roundtrip A, saveRDS -> rm/gc -> readRDS                -- "make it its own object"
#   C rebuilt   CreateSeuratObject on the CM counts, genes    -- the strongest form of it:
#               with zero CM counts dropped, nothing else     nothing inherited at all
#               carried over
#
# EXPECTATION, stated up front so the result is interpretable either way: A and B should be
# IDENTICAL. saveRDS serialises by value; a Seurat subset holds no live reference to its
# parent, so a round-trip cannot drop a hidden link because there is no hidden link. If B
# differs from A, something non-deterministic is in the pipeline and that is the finding.
#
# C is where a real difference CAN appear, and it has nothing to do with hidden references:
# dropping CM-zero genes changes the gene universe SCTransform sees, which moves HVG
# ranking -> PCs -> Harmony -> UMAP. That is a genuine sensitivity worth knowing.
#
# DIAGNOSTIC ONLY. Writes new tables/figures; never touches seurat.cm.subclustered.rds.
#   Reads  processing/seurat.combined.annotated.rds
#   Writes results/tables/cm_objectmode_comparison.csv
#          results/tables/cm_objectmode_percell.csv.gz
#          results/figures/cm_objectmode_umap.png
#          processing/cm_objectmode_embeddings.rds   (slim: embeddings + labels, not objects)
#
# Runs in e2f-seurat-full:latest (our_analysis/Dockerfile.seurat). NOTE the package
# versions there are not the ones the shipped object was built with -- see SHIPPED below.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(harmony); library(ggplot2); library(patchwork) }))

RES_SWEEP <- c(0.1, 0.2, 0.3, 0.4, 0.6)   # identical to cm_subcluster_build.R
PC_DIMS   <- 1:30                         # identical to cm_subcluster_build.R
SEED      <- 42L
VARIANTS  <- c("filter", "roundtrip", "rebuilt")
RES_COLS  <- paste0("SCT_snn_res.", RES_SWEEP)
TMP_RDS   <- file.path(PROC, "cm_objectmode_tmp_subset.rds")   # deleted at the end
# Compare against the shipped production object too. Reported SEPARATELY and labelled,
# because it was built on Windows / R 4.5.2 / unpinned Seurat and this container is not
# that -- a difference there is confounded with package version and proves nothing about
# object construction. The A/B/C comparisons are all same-session and are not confounded.
SHIPPED   <- file.path(PROC, "seurat.cm.subclustered.rds")

# cross-version counts accessor (slot= in Seurat v4, layer= in v5) -- same shim as
# cm_subcluster_analyze.R:30-33 and cm_cycling_investigate.R:35-37
getdata <- function(o, assay, what) tryCatch(
  SeuratObject::GetAssayData(o, assay = assay, slot  = what),
  error = function(e) SeuratObject::GetAssayData(o, assay = assay, layer = what))

## --- helpers: comparison maths ---------------------------------------------
# PC/Harmony component signs are arbitrary (an eigenvector and its negation are the same
# axis). Comparing without aligning sign reports large fake differences, so align first.
align_sign <- function(X, Y) {
  s <- sign(colSums(X * Y)); s[s == 0] <- 1
  sweep(Y, 2, s, `*`)
}
cmp_embed <- function(X, Y) {
  n <- min(ncol(X), ncol(Y)); X <- X[, seq_len(n), drop = FALSE]; Y <- Y[, seq_len(n), drop = FALSE]
  Ya <- align_sign(X, Y)
  list(max_abs_diff = max(abs(X - Ya)),
       min_abs_cor  = min(abs(vapply(seq_len(n), function(j) suppressWarnings(
                        stats::cor(X[, j], Ya[, j])), numeric(1))), na.rm = TRUE))
}
# Procrustes: rotate/reflect Y onto X after centring and scaling both to unit Frobenius
# norm. Reflection is allowed deliberately -- a mirrored UMAP is the same embedding.
# Implemented here rather than via vegan so the image needs no extra dependency.
procrustes_rmsd <- function(X, Y) {
  X <- scale(X, TRUE, FALSE); Y <- scale(Y, TRUE, FALSE)
  X <- X / sqrt(sum(X^2));   Y <- Y / sqrt(sum(Y^2))
  s <- svd(crossprod(Y, X)); R <- s$u %*% t(s$v)
  sqrt(mean(rowSums((X - Y %*% R)^2)))
}
ari <- function(a, b) {                       # adjusted Rand index; label-invariant
  tab <- table(a, b); n <- sum(tab)
  cc <- function(x) sum(choose(x, 2))
  idx <- cc(tab); ea <- cc(rowSums(tab)); eb <- cc(colSums(tab))
  expct <- ea * eb / choose(n, 2); mx <- (ea + eb) / 2
  if (isTRUE(all.equal(mx, expct))) return(1)
  (idx - expct) / (mx - expct)
}
# Fraction of cells that keep their partner cluster under a GREEDY max-overlap matching.
# Greedy, not Hungarian -- it can be very slightly pessimistic, which is the safe direction
# for a check that is trying to detect disagreement.
greedy_agree <- function(a, b) {
  tab <- table(a, b); tot <- sum(tab); hit <- 0
  while (nrow(tab) > 0 && ncol(tab) > 0 && max(tab) > 0) {
    w <- which(tab == max(tab), arr.ind = TRUE)[1, ]
    hit <- hit + tab[w[1], w[2]]
    tab <- tab[-w[1], -w[2], drop = FALSE]
  }
  hit / tot
}
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

## --- the pipeline under test (byte-identical to cm_subcluster_build.R:35-41) ------
# set.seed before EVERY stage, so each variant enters each stage with the same RNG state
# no matter what the preceding stages consumed. RunHarmony seeds its k-means from the R
# RNG and takes no seed argument of its own, so this is the only way to pin it.
embed_cm <- function(obj, tag) {
  # Production sets this at cm_subcluster_build.R:35 and it is load-bearing: after
  # combined.R the parent's DefaultAssay is "SCT", so variants A and B inherit "SCT"
  # and would transform the whole-heart residuals instead of the RNA counts. C is
  # already "RNA" from CreateSeuratObject; setting it for all three keeps the arms
  # comparable.
  DefaultAssay(obj) <- "RNA"
  cat(sprintf("  [%s] SCTransform ...\n", tag)); set.seed(SEED)
  obj <- SCTransform(obj, method = "glmGamPoi", conserve.memory = TRUE,
                     verbose = FALSE, seed.use = SEED)
  cat(sprintf("  [%s] RunPCA ...\n", tag)); set.seed(SEED)
  obj <- RunPCA(obj, verbose = FALSE, seed.use = SEED)
  cat(sprintf("  [%s] RunHarmony ...\n", tag)); set.seed(SEED)
  obj <- RunHarmony(obj, group.by.vars = "orig.ident")
  cat(sprintf("  [%s] RunUMAP ...\n", tag)); set.seed(SEED)
  obj <- RunUMAP(obj, reduction = "harmony", dims = PC_DIMS,
                 reduction.name = "umap", verbose = FALSE, seed.use = SEED)
  cat(sprintf("  [%s] FindNeighbors/FindClusters ...\n", tag)); set.seed(SEED)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = PC_DIMS, verbose = FALSE)
  set.seed(SEED)
  FindClusters(obj, resolution = RES_SWEEP, verbose = FALSE, random.seed = SEED)
}
# Keep only what the comparison needs. Three full objects would be ~7 GB of RDS for no
# reason; the embeddings and labels are all any downstream check or the app tab reads.
slim_of <- function(obj) list(
  cells    = colnames(obj),
  genes    = rownames(obj[["RNA"]]),
  hvg      = VariableFeatures(obj),
  pca      = Embeddings(obj, "pca")[, PC_DIMS, drop = FALSE],
  harmony  = Embeddings(obj, "harmony")[, PC_DIMS, drop = FALSE],
  umap     = Embeddings(obj, "umap")[, 1:2, drop = FALSE],
  clusters = obj@meta.data[, RES_COLS, drop = FALSE],
  nCount   = obj$nCount_RNA, nFeature = obj$nFeature_RNA)

## --- load parent and pin its state -----------------------------------------
inp <- file.path(PROC, "seurat.combined.annotated.rds")
cat(sprintf("reading %s\n  (mtime %s, %s)\n", basename(inp),
            format(file.mtime(inp), "%Y-%m-%d %H:%M:%S"),
            format(structure(file.size(inp), class = "object_size"), units = "auto")))
comb <- readRDS(inp)
stopifnot("celltype" %in% colnames(comb@meta.data))
comb$genotype <- genotype_of(comb$orig.ident)
cat(sprintf("combined cells: %d\n", ncol(comb))); print(table(comb$celltype))

## --- A: filter (exactly what production does) -------------------------------
cat("\n=== A: filter ===\n")
cmA <- comb[, comb$celltype == "Cardiomyocyte"]
cat(sprintf("cardiomyocytes: %d cells x %d genes\n", ncol(cmA), nrow(cmA)))
print(table(cmA$orig.ident, cmA$lane))
prov <- list(input = basename(inp), input_mtime = format(file.mtime(inp)),
             combined_cells = ncol(comb), cm_cells = ncol(cmA), cm_genes = nrow(cmA))
# The input cell set, captured BEFORE any embedding. Variant C must be built from
# these, not from resA$cells: SCTransform can drop a cell, and seeding C from A's
# surviving cells would quietly give the arms different inputs.
cells_in <- colnames(cmA)
rm(comb); gc(verbose = FALSE)

saveRDS(cmA, TMP_RDS)                       # variant B reads this back
resA <- slim_of(embed_cm(cmA, "filter"))
rm(cmA); gc(verbose = FALSE)

## --- B: round-trip through disk ---------------------------------------------
cat("\n=== B: roundtrip (saveRDS -> rm/gc -> readRDS) ===\n")
cmB <- readRDS(TMP_RDS)
resB <- slim_of(embed_cm(cmB, "roundtrip"))
rm(cmB); gc(verbose = FALSE)

## --- C: rebuilt from counts, nothing inherited ------------------------------
cat("\n=== C: rebuilt (CreateSeuratObject on CM counts) ===\n")
cmS <- readRDS(TMP_RDS)
cts  <- getdata(cmS, "RNA", "counts")
meta <- cmS@meta.data
rm(cmS); gc(verbose = FALSE)
# Cell ORDER must match A exactly. SCTransform subsamples cells (ncells = 5000) to fit its
# parameters; a different column order changes that subsample and would surface as a
# difference that is really just RNG.
cts  <- cts[, cells_in, drop = FALSE]
meta <- meta[cells_in, , drop = FALSE]
keep <- Matrix::rowSums(cts) > 0
cat(sprintf("genes with zero counts in CMs: %d of %d dropped (%d kept)\n",
            sum(!keep), length(keep), sum(keep)))
# Drop nCount/nFeature and let CreateSeuratObject recompute them -- then assert they match
# A. They must: every dropped gene is all-zero, so it contributed nothing to either.
meta <- meta[, setdiff(colnames(meta), c("nCount_RNA", "nFeature_RNA")), drop = FALSE]
cmC  <- CreateSeuratObject(counts = cts[keep, , drop = FALSE], meta.data = meta,
                           project = "CM_rebuilt", min.cells = 0, min.features = 0)
rm(cts, meta); gc(verbose = FALSE)
stopifnot(identical(colnames(cmC), cells_in))
resC <- slim_of(embed_cm(cmC, "rebuilt"))
rm(cmC); gc(verbose = FALSE)
unlink(TMP_RDS)

res <- list(filter = resA, roundtrip = resB, rebuilt = resC)
# If SCTransform retained different cells in different arms, every comparison below is
# meaningless -- fail here rather than report a misaligned diff as a finding.
stopifnot(identical(res$filter$cells, res$roundtrip$cells),
          identical(res$filter$cells, res$rebuilt$cells))
cat(sprintf("\nnCount_RNA identical A vs C: %s\nnFeature_RNA identical A vs C: %s\n",
            isTRUE(all.equal(resA$nCount,   resC$nCount,   tolerance = 0)),
            isTRUE(all.equal(resA$nFeature, resC$nFeature, tolerance = 0))))

## --- comparison table --------------------------------------------------------
row <- function(pair, metric, value, note = "") data.frame(
  pair = pair, metric = metric, value = value, note = note, stringsAsFactors = FALSE)
out <- list()
for (p in list(c("filter","roundtrip"), c("filter","rebuilt"), c("roundtrip","rebuilt"))) {
  a <- res[[p[1]]]; b <- res[[p[2]]]; tag <- paste(p, collapse = " vs ")
  out <- c(out, list(
    row(tag, "cells_identical",   as.numeric(identical(a$cells, b$cells))),
    row(tag, "genes_a",           length(a$genes)),
    row(tag, "genes_b",           length(b$genes)),
    row(tag, "hvg_jaccard",       jaccard(a$hvg, b$hvg)),
    row(tag, "hvg_n_shared",      length(intersect(a$hvg, b$hvg))),
    row(tag, "pca_max_abs_diff",     cmp_embed(a$pca, b$pca)$max_abs_diff),
    row(tag, "pca_min_abs_cor",      cmp_embed(a$pca, b$pca)$min_abs_cor,
        "sign-aligned; min over the 30 PCs"),
    row(tag, "harmony_max_abs_diff", cmp_embed(a$harmony, b$harmony)$max_abs_diff),
    row(tag, "harmony_min_abs_cor",  cmp_embed(a$harmony, b$harmony)$min_abs_cor,
        "sign-aligned; min over the 30 dims"),
    row(tag, "umap_max_abs_diff",    max(abs(a$umap - b$umap)), "raw, no alignment"),
    row(tag, "umap_procrustes_rmsd", procrustes_rmsd(a$umap, b$umap),
        "0 = same shape; reflection allowed")))
  for (rc in RES_COLS) {
    ca <- as.character(a$clusters[[rc]]); cb <- as.character(b$clusters[[rc]])
    out <- c(out, list(
      row(tag, paste0("ARI_", rc),        ari(ca, cb)),
      row(tag, paste0("nclust_a_", rc),   length(unique(ca))),
      row(tag, paste0("nclust_b_", rc),   length(unique(cb))),
      row(tag, paste0("agree_", rc),      greedy_agree(ca, cb), "greedy max-overlap match"),
      row(tag, paste0("ncells_moved_", rc), round((1 - greedy_agree(ca, cb)) * length(ca)),
          "greedy max-overlap match")))
  }
}

# Shipped object -- separate, and labelled as confounded.
if (file.exists(SHIPPED)) {
  cat("\n=== shipped production object (version-confounded comparison) ===\n")
  sh <- readRDS(SHIPPED)
  shu <- Embeddings(sh, "umap")[, 1:2, drop = FALSE]
  shc <- sh@meta.data[, RES_COLS, drop = FALSE]
  common <- intersect(rownames(shu), resA$cells)
  cat(sprintf("shipped cells: %d; shared with this run: %d\n", nrow(shu), length(common)))
  out <- c(out, list(
    row("filter vs SHIPPED", "cells_shared", length(common),
        "CONFOUNDED: shipped was built on R 4.5.2/Windows with different package versions"),
    row("filter vs SHIPPED", "cells_identical",
        as.numeric(setequal(rownames(shu), resA$cells)), "CONFOUNDED"),
    row("filter vs SHIPPED", "umap_procrustes_rmsd",
        procrustes_rmsd(resA$umap[common, ], shu[common, ]), "CONFOUNDED")))
  for (rc in RES_COLS) out <- c(out, list(
    row("filter vs SHIPPED", paste0("ARI_", rc),
        ari(as.character(resA$clusters[common, rc]), as.character(shc[common, rc])),
        "CONFOUNDED")))
  rm(sh); gc(verbose = FALSE)
} else cat("\n(shipped object absent -- skipping that comparison)\n")

cmp <- do.call(rbind, out)
cmp$value <- ifelse(abs(cmp$value) >= 1e4 | cmp$value == 0, cmp$value, signif(cmp$value, 6))
write.csv(cmp, file.path(OUTTAB, "cm_objectmode_comparison.csv"), row.names = FALSE)

## --- per-cell coordinates (feeds the app tab) --------------------------------
percell <- do.call(rbind, lapply(VARIANTS, function(v) {
  d <- data.frame(cell = res[[v]]$cells, variant = v,
                  UMAP1 = res[[v]]$umap[, 1], UMAP2 = res[[v]]$umap[, 2],
                  stringsAsFactors = FALSE)
  cbind(d, res[[v]]$clusters) }))
gz <- gzfile(file.path(OUTTAB, "cm_objectmode_percell.csv.gz"), "w")
write.csv(percell, gz, row.names = FALSE); close(gz)

saveRDS(c(res, list(provenance = prov, comparison = cmp)),
        file.path(PROC, "cm_objectmode_embeddings.rds"))

## --- figure ------------------------------------------------------------------
ps <- lapply(VARIANTS, function(v) {
  d <- data.frame(res[[v]]$umap); names(d) <- c("UMAP1", "UMAP2")
  d$sub <- factor(paste0("CM", res[[v]]$clusters[["SCT_snn_res.0.2"]]))
  ggplot(d, aes(UMAP1, UMAP2, colour = sub)) +
    geom_point(size = .2, alpha = .6, shape = 16) + theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    ggtitle(sprintf("%s (%d subclusters, res 0.2)", v, nlevels(d$sub)))
})
png(file.path(OUTFIG, "cm_objectmode_umap.png"), 1650, 560)
print(wrap_plots(ps, nrow = 1)); dev.off()

## --- console verdict ---------------------------------------------------------
say <- function(pair) {
  g <- function(m) cmp$value[cmp$pair == pair & cmp$metric == m][1]
  cat(sprintf("  %-22s UMAP max|diff| %-11s procrustes %-11s ARI(res0.2) %s\n",
              pair, signif(g("umap_max_abs_diff"), 4), signif(g("umap_procrustes_rmsd"), 4),
              signif(g("ARI_SCT_snn_res.0.2"), 6)))
}
cat("\n=== VERDICT ===\n")
for (p in c("filter vs roundtrip", "filter vs rebuilt", "roundtrip vs rebuilt")) say(p)
ident <- cmp$value[cmp$pair == "filter vs roundtrip" & cmp$metric == "umap_max_abs_diff"][1]
cat(sprintf("\nA vs B identical: %s\n", if (identical(ident, 0) || ident == 0) "YES" else
            paste0("NO (max|diff| = ", signif(ident, 6), ") -- investigate before trusting A vs C")))
cat("\nwrote:\n  ", file.path(OUTTAB, "cm_objectmode_comparison.csv"),
    "\n  ", file.path(OUTTAB, "cm_objectmode_percell.csv.gz"),
    "\n  ", file.path(OUTFIG, "cm_objectmode_umap.png"),
    "\n  ", file.path(PROC, "cm_objectmode_embeddings.rds"), "\n", sep = "")
cat("=== DONE cm_objectmode_check ===\n")
