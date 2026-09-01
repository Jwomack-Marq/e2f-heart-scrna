#!/usr/bin/env Rscript
# Phase A -- build ONE 4-group object (P0/P7 x KO/WT) and Harmony-integrate it for
# joint annotation and cross-group views. Integration is for the EMBEDDING ONLY;
# all DE stays on RNA/pseudobulk (see de_descriptive.R / cross_timepoint.R).
# Reads the 4 processing/merge.lanes.*.rds; writes seurat.combined.rds + UMAPs.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(harmony); library(ggplot2) }))

groups <- c("P0WT","P0KO","P7WT","P7KO")
objs <- lapply(groups, function(g) {
  o <- readRDS(file.path(PROC, paste0("merge.lanes.", g, ".rds")))
  o$orig.ident <- g
  o$timepoint  <- substr(g, 1, 2)
  o$genotype   <- genotype_of(g)
  o
})
names(objs) <- groups
comb <- merge(objs[[1]], objs[-1], add.cell.ids = groups, project = "Han_combined")
if (utils::packageVersion("SeuratObject") >= "5.0.0") comb <- JoinLayers(comb)
rm(objs); gc(verbose = FALSE)
cat(sprintf("combined cells: %d\n", ncol(comb))); print(table(comb$orig.ident))

DefaultAssay(comb) <- "RNA"
comb <- SCTransform(comb, method = "glmGamPoi", conserve.memory = TRUE, verbose = FALSE)
comb <- RunPCA(comb, verbose = FALSE)

# Pre-integration UMAP (to show the group/batch structure) ...
comb <- RunUMAP(comb, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
saveRDS(comb, file.path(PROC, "seurat.combined.prePCAcheckpoint.rds"))  # checkpoint after SCT/PCA
# ... then Harmony across the 4 libraries so cell types co-embed for annotation.
# (harmony >=2.0 dropped `assay.use`; it uses the active reduction "pca" by default.)
comb <- RunHarmony(comb, group.by.vars = "orig.ident")
comb <- comb %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap")

saveRDS(comb, file.path(PROC, "seurat.combined.rds"))

p_un <- DimPlot(comb, reduction = "umap.unintegrated", group.by = "orig.ident") + ggtitle("Pre-Harmony (by library)")
p_in <- DimPlot(comb, reduction = "umap",               group.by = "orig.ident") + ggtitle("Post-Harmony (by library)")
png(file.path(OUTFIG, "combined_harmony_before_after.png"), 1200, 520); print(p_un | p_in); dev.off()
png(file.path(OUTFIG, "combined_umap_facets.png"), 1300, 900)
print((DimPlot(comb, group.by = "seurat_clusters", label = TRUE) + NoLegend() + ggtitle("clusters")) |
      (DimPlot(comb, group.by = "timepoint") + ggtitle("timepoint")) /
      (DimPlot(comb, group.by = "genotype")  + ggtitle("genotype")))
dev.off()

cat(sprintf("clusters: %d\n", length(unique(comb$seurat_clusters))))
cat("=== DONE combined ===\n")
