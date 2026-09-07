#!/usr/bin/env Rscript
# Trajectory / maturation pseudotime on cardiomyocytes (DESCRIPTIVE, n=1, sex-confounded).
# Slingshot on the Harmony embedding -> orders CMs along an immature->mature axis; we
# validate the direction with maturation markers and compare KO vs WT pseudotime.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
if (!requireNamespace("slingshot", quietly = TRUE)) {
  writeLines("slingshot not installed; trajectory skipped.", file.path(OUTTAB, "trajectory_blocked.txt"))
  cat("slingshot missing.\n=== DONE trajectory_slingshot ===\n"); quit(save = "no")
}
suppressWarnings(suppressMessages({ library(slingshot); library(ggplot2) }))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(comb) <- "RNA"; comb <- NormalizeData(comb, verbose = FALSE)
cm <- comb[, comb$celltype == "Cardiomyocyte"]
cm$genotype <- genotype_of(cm$orig.ident)
cm <- AddModuleScore(cm, list(intersect(CM_MATURE, rownames(cm))),   name = "Mature")
cm <- AddModuleScore(cm, list(intersect(CM_IMMATURE, rownames(cm))), name = "Immature")
cm$cycling <- cm$Phase %in% c("S","G2M")

emb <- Embeddings(cm, "harmony")[, 1:30]
cl  <- droplevels(factor(cm$seurat_clusters))
# root = CM cluster with highest immaturity (immature score + cycling fraction, z-scored)
imm <- tapply(cm$Immature1, cl, mean); cyc <- tapply(cm$cycling, cl, mean)
score <- scale(imm) + scale(cyc)
root <- names(score[,1])[which.max(score[,1])]
cat("root (most immature/cycling) cluster:", root, "\n")

sds <- slingshot::slingshot(emb, clusterLabels = as.character(cl), start.clus = root)
pt  <- slingshot::slingPseudotime(sds)
cm$pseudotime <- rowMeans(pt, na.rm = TRUE)              # average across lineages
# orient so higher pseudotime = more mature (flip if anti-correlated with Mature score)
if (cor(cm$pseudotime, cm$Mature1, use = "complete.obs") < 0) cm$pseudotime <- max(cm$pseudotime, na.rm=TRUE) - cm$pseudotime

## validation: pseudotime vs maturation/immaturity/cycling --------------------
val <- data.frame(
  metric = c("Mature score","Immature score","cycling fraction"),
  cor_with_pseudotime = c(round(cor(cm$pseudotime, cm$Mature1, use="complete.obs"),3),
                          round(cor(cm$pseudotime, cm$Immature1, use="complete.obs"),3),
                          round(cor(cm$pseudotime, as.numeric(cm$cycling), use="complete.obs"),3)))
write.csv(val, file.path(OUTTAB, "pseudotime_validation.csv"), row.names = FALSE)
cat("validation (expect Mature +, Immature -, cycling -):\n"); print(val, row.names = FALSE)

## KO vs WT pseudotime per timepoint (descriptive) ---------------------------
md <- cm@meta.data
summ <- md |> dplyr::group_by(timepoint, genotype) |>
  dplyr::summarise(n = dplyr::n(), median_pt = round(median(pseudotime, na.rm=TRUE),3),
                   mean_pt = round(mean(pseudotime, na.rm=TRUE),3),
                   pct_late = round(100*mean(pseudotime > median(md$pseudotime, na.rm=TRUE), na.rm=TRUE),1),
                   .groups="drop")
summ$NOTE <- "descriptive (n=1, sex-confounded)"
write.csv(summ, file.path(OUTTAB, "pseudotime_KOvsWT.csv"), row.names = FALSE)
write.csv(data.frame(cell=colnames(cm), pseudotime=cm$pseudotime, timepoint=cm$timepoint,
                     genotype=cm$genotype, Phase=cm$Phase),
          file.path(OUTTAB, "pseudotime_per_cell.csv"), row.names = FALSE)
cat("\nKO vs WT pseudotime:\n"); print(as.data.frame(summ), row.names = FALSE)

## figures -------------------------------------------------------------------
um <- as.data.frame(Embeddings(cm, "umap")); names(um) <- c("UMAP1","UMAP2")
um$pseudotime <- cm$pseudotime; um$genotype <- cm$genotype; um$timepoint <- cm$timepoint
png(file.path(OUTFIG, "traj_umap_pseudotime.png"), 900, 700)
print(ggplot(um, aes(UMAP1, UMAP2, color = pseudotime)) + geom_point(size=.3) +
        scale_color_viridis_c() + theme_bw() + ggtitle("Cardiomyocyte maturation pseudotime (Slingshot)"))
dev.off()
png(file.path(OUTFIG, "traj_pseudotime_KOvsWT.png"), 950, 500)
print(ggplot(um, aes(pseudotime, fill = genotype)) + geom_density(alpha=.5) +
        facet_wrap(~ timepoint) + scale_fill_manual(values=c(WT="steelblue",KO="firebrick")) +
        theme_bw() + ggtitle("CM pseudotime KO vs WT (descriptive, n=1)"))
dev.off()
png(file.path(OUTFIG, "traj_marker_validation.png"), 950, 450)
vdf <- data.frame(pseudotime=cm$pseudotime, Mature=cm$Mature1, Immature=cm$Immature1)
vdf <- tidyr::pivot_longer(vdf, -pseudotime, names_to="score", values_to="value")
print(ggplot(vdf, aes(pseudotime, value)) + geom_smooth(method="loess", se=FALSE) +
        facet_wrap(~score, scales="free_y") + theme_bw() +
        ggtitle("Pseudotime vs maturation markers (direction check)"))
dev.off()
cat("\n=== DONE trajectory_slingshot ===\n")
