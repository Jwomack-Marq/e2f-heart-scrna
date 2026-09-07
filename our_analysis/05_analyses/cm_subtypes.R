#!/usr/bin/env Rscript
# Cardiomyocyte SUBTYPING (the prior pipeline's CM subtypes; our annotation only had
# one "Cardiomyocyte" class). Marker-module argmax per CM cluster -> ventricular /
# atrial / trabecular / compact / cycling. Then DESCRIPTIVE subtype composition and
# cycling fraction KO vs WT per timepoint. Panels from their _env.R CARDIAC_MARKERS.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) }))

CM_SUBTYPES <- list(
  Ventricular = c("Myl2","Myh7"),
  Atrial      = c("Myl7","Sln","Nppa"),
  Trabecular  = c("Bmp10","Nppa","Hey2"),
  Compact     = c("Hey2","Irx3","Tbx20"),
  Cycling     = c("Mki67","Top2a","Ccnb1","Aurkb","Cdca8"))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(comb) <- "RNA"; comb <- NormalizeData(comb, verbose = FALSE)
comb$genotype <- genotype_of(comb$orig.ident)
cm <- comb[, comb$celltype == "Cardiomyocyte"]
cm$cl <- droplevels(factor(cm$seurat_clusters))

ms.cols <- c()
for (s in names(CM_SUBTYPES)) {
  mk <- intersect(CM_SUBTYPES[[s]], rownames(cm)); if (length(mk) < 1) next
  cm <- AddModuleScore(cm, features = list(mk), name = paste0("cms_", s))
  ms.cols[s] <- paste0("cms_", s, "1")
}
cl.mean <- sapply(ms.cols, function(col) tapply(cm[[col]][,1], cm$cl, mean))
sub.lab <- colnames(cl.mean)[max.col(cl.mean, ties.method = "first")]; names(sub.lab) <- rownames(cl.mean)
cm$cm_subtype <- unname(sub.lab[as.character(cm$cl)])
cat("CM subtype per cluster:\n"); print(data.frame(cluster = rownames(cl.mean), subtype = sub.lab), row.names = FALSE)

## composition (descriptive) -------------------------------------------------
tab <- as.data.frame.matrix(table(cm$cm_subtype, cm$orig.ident))
prop <- sweep(tab, 2, colSums(tab), "/")
comp <- data.frame(cm_subtype = rownames(prop), round(prop, 4), check.names = FALSE)
for (tp in TIMEPOINTS) { wt <- paste0(tp,"WT"); ko <- paste0(tp,"KO")
  if (all(c(wt,ko) %in% colnames(prop))) comp[[paste0(tp,"_log2FC_KOvsWT")]] <- round(log2((prop[[ko]]+1e-4)/(prop[[wt]]+1e-4)),3) }
write.csv(comp, file.path(OUTTAB, "cm_subtype_composition.csv"), row.names = FALSE)

## cycling fraction by subtype (descriptive) ---------------------------------
cm$cycling <- cm$Phase %in% c("S","G2M")
cyc <- cm@meta.data |> group_by(timepoint, cm_subtype, genotype) |>
  summarise(n = dplyr::n(), pct_cycling = round(100*mean(cycling),1), .groups = "drop")
write.csv(cyc, file.path(OUTTAB, "cm_subtype_cycling.csv"), row.names = FALSE)

## figures -------------------------------------------------------------------
if ("umap" %in% names(cm@reductions)) {
  png(file.path(OUTFIG, "cm_subtypes_umap.png"), 900, 700)
  print(DimPlot(cm, group.by = "cm_subtype", reduction = "umap", label = TRUE) + ggtitle("Cardiomyocyte subtypes"))
  dev.off()
}
flat <- unique(unlist(lapply(CM_SUBTYPES, function(g) intersect(g, rownames(cm)))))
png(file.path(OUTFIG, "cm_subtypes_dotplot.png"), 950, 500)
print(DotPlot(cm, features = flat, group.by = "cm_subtype") + RotatedAxis() + ggtitle("CM subtype markers"))
dev.off()
longp <- pivot_longer(data.frame(cm_subtype = rownames(prop), prop, check.names = FALSE), -cm_subtype,
                      names_to = "group", values_to = "prop")
png(file.path(OUTFIG, "cm_subtype_composition.png"), 850, 600)
print(ggplot(longp, aes(group, prop, fill = cm_subtype)) + geom_col() + theme_bw() +
        ylab("proportion of cardiomyocytes") + ggtitle("CM subtype composition (descriptive, n=1)"))
dev.off()

cat("\n--- CM subtype composition (KO/WT log2FC) ---\n"); print(comp, row.names = FALSE)
cat("\n--- cycling fraction by CM subtype ---\n"); print(as.data.frame(cyc), row.names = FALSE)
cat("\n=== DONE cm_subtypes ===\n")
