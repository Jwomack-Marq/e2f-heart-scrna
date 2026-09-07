#!/usr/bin/env Rscript
# Cardiomyocyte SUBCLUSTERING -- build step (run ONCE; expensive).
# True re-clustering of the CM compartment (the existing cm_subtypes.R only assigns
# marker-module subtype LABELS, it does not re-cluster). We subset cardiomyocytes from
# the Harmony-integrated 4-group object, re-embed on the CM subset so structure reflects
# CM-INTERNAL variation, sweep a few resolutions, and save the object for the (cheap,
# re-runnable) analyze step. Embedding/clustering only -- all DE stays on RNA/pseudobulk.
#   Reads  processing/seurat.combined.annotated.rds
#   Writes processing/seurat.cm.subclustered.rds + a resolution-sweep UMAP/clustree.
# DESCRIPTIVE pilot (n=1, sex-confounded) -- see REPLICATES.md.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(harmony); library(ggplot2); library(patchwork) }))

RES_SWEEP <- c(0.1, 0.2, 0.3, 0.4, 0.6)   # subcluster resolutions to store (pick one in analyze step)
PC_DIMS   <- 1:30

## --- load + subset to cardiomyocytes ---------------------------------------
comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
comb$genotype <- genotype_of(comb$orig.ident)
stopifnot("celltype" %in% colnames(comb@meta.data))
cat(sprintf("combined cells: %d\n", ncol(comb))); print(table(comb$celltype))

cm <- comb[, comb$celltype == "Cardiomyocyte"]
cat(sprintf("\ncardiomyocytes retained: %d\n", ncol(cm)))
print(table(cm$orig.ident, cm$lane))
keep_meta <- intersect(c("Phase","S.Score","G2M.Score"), colnames(cm@meta.data))  # carry over if present
rm(comb); gc(verbose = FALSE)

## --- re-embed on the CM subset ---------------------------------------------
# Recompute SCT/PCA on JUST the CMs so variable genes & PCs capture CM-internal
# heterogeneity rather than whole-heart lineage differences. Harmony again removes
# the library/batch axis so subclusters co-embed across P0/P7 x KO/WT.
DefaultAssay(cm) <- "RNA"
cm <- SCTransform(cm, method = "glmGamPoi", conserve.memory = TRUE, verbose = FALSE)
cm <- RunPCA(cm, verbose = FALSE)
cm <- RunHarmony(cm, group.by.vars = "orig.ident")
cm <- RunUMAP(cm, reduction = "harmony", dims = PC_DIMS, reduction.name = "umap", verbose = FALSE)
cm <- FindNeighbors(cm, reduction = "harmony", dims = PC_DIMS, verbose = FALSE)
cm <- FindClusters(cm, resolution = RES_SWEEP, verbose = FALSE)   # stores SCT_snn_res.<r>

res_cols <- paste0("SCT_snn_res.", RES_SWEEP)
cat("\n--- clusters per resolution ---\n")
for (rc in res_cols) cat(sprintf("  %-18s %d clusters\n", rc, length(unique(cm[[rc]][,1]))))

## --- cell-cycle scoring (carry over, else compute on the CM subset) ---------
if (all(c("Phase","S.Score","G2M.Score") %in% keep_meta)) {
  cat("\nPhase/S.Score/G2M.Score carried over from combined object.\n")
} else {
  cc <- cc_lists()
  cm <- CellCycleScoring(cm, s.features = intersect(cc$S, rownames(cm)),
                         g2m.features = intersect(cc$G2M, rownames(cm)), set.ident = FALSE)
  cat("\nCellCycleScoring run on CM subset (mouse-mapped cc.genes).\n")
}

saveRDS(cm, file.path(PROC, "seurat.cm.subclustered.rds"))

## --- resolution-sweep figures (decision aids for 'how many subclusters') ----
ps <- lapply(seq_along(RES_SWEEP), function(i)
  DimPlot(cm, group.by = res_cols[i], reduction = "umap", label = TRUE) + NoLegend() +
    ggtitle(sprintf("res %.1f (%d)", RES_SWEEP[i], length(unique(cm[[res_cols[i]]][,1])))))
png(file.path(OUTFIG, "cm_subcluster_resolution_sweep_umap.png"), 1500, 950)
print(wrap_plots(ps, ncol = 3)); dev.off()

if (requireNamespace("clustree", quietly = TRUE)) {
  png(file.path(OUTFIG, "cm_subcluster_clustree.png"), 950, 1000)
  print(clustree::clustree(cm, prefix = "SCT_snn_res.")); dev.off()
  cat("clustree written.\n")
} else message("clustree not installed -- skipped (resolution-sweep UMAP still written).")

cat("\nsaved: ", file.path(PROC, "seurat.cm.subclustered.rds"), "\n")
cat("=== DONE cm_subcluster_build ===\n")
