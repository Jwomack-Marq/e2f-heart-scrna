#!/usr/bin/env Rscript
# Fuller cell-cycle analysis (DESCRIPTIVE, n=1, sex-confounded).
#   1) CellCycleScoring on the combined annotated object (one run across all 4 groups)
#   2) cycling fraction (S/G2M) KO vs WT per cell type x timepoint  (not just CMs)
#   3) Phase (G1/S/G2M) composition per group
#   4) split cardiomyocytes into cycling vs non-cycling and test that compartment
#      KO vs WT with propeller -- the cleanest "failed cell-cycle exit" readout
# Reads seurat.combined.annotated.rds; writes results/{tables,figures}/.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(ggplot2); library(tidyr); library(dplyr) }))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
DefaultAssay(comb) <- "RNA"; comb <- NormalizeData(comb, verbose = FALSE)
comb$genotype <- genotype_of(comb$orig.ident)

## 1) cell-cycle scoring (mouse-mapped cc.genes) ------------------------------
cc <- cc_lists()
comb <- CellCycleScoring(comb, s.features = intersect(cc$S, rownames(comb)),
                         g2m.features = intersect(cc$G2M, rownames(comb)), set.ident = FALSE)
comb$cycling <- comb$Phase %in% c("S", "G2M")
saveRDS(comb, file.path(PROC, "seurat.combined.annotated.rds"))   # persist Phase/scores

md <- comb@meta.data

## 2) cycling fraction by cell type x timepoint x genotype --------------------
cyc <- md |>
  group_by(timepoint, celltype, genotype) |>
  summarise(n = dplyr::n(),
            pct_S    = round(100*mean(Phase == "S"), 1),
            pct_G2M  = round(100*mean(Phase == "G2M"), 1),
            pct_cycling = round(100*mean(cycling), 1), .groups = "drop")
write.csv(cyc, file.path(OUTTAB, "cellcycle_fraction_by_celltype.csv"), row.names = FALSE)

p1 <- ggplot(cyc, aes(reorder(celltype, pct_cycling), pct_cycling, fill = genotype)) +
  geom_col(position = "dodge") + coord_flip() + facet_wrap(~ timepoint) +
  scale_fill_manual(values = c(WT = "steelblue", KO = "firebrick")) +
  labs(title = "Cycling fraction (S/G2M) KO vs WT by cell type (descriptive, n=1)",
       x = NULL, y = "% cycling") + theme_bw()
ggsave(file.path(OUTFIG, "cellcycle_fraction_by_celltype.png"), p1, width = 10, height = 6, dpi = 120)

## 3) Phase composition per group --------------------------------------------
ph <- md |> count(orig.ident, timepoint, genotype, Phase) |>
  group_by(orig.ident) |> mutate(frac = n/sum(n)) |> ungroup()
write.csv(ph, file.path(OUTTAB, "cellcycle_phase_composition.csv"), row.names = FALSE)
p2 <- ggplot(ph, aes(orig.ident, frac, fill = Phase)) + geom_col() +
  labs(title = "Cell-cycle phase composition per group (descriptive)", x = NULL, y = "fraction") + theme_bw()
ggsave(file.path(OUTFIG, "cellcycle_phase_composition.png"), p2, width = 7, height = 5, dpi = 120)

## 4) cycling-CM compartment KO vs WT (the cell-cycle-exit test) --------------
md$celltype_cc <- ifelse(md$celltype == "Cardiomyocyte",
                         ifelse(md$cycling, "Cardiomyocyte_cycling", "Cardiomyocyte_G1"),
                         md$celltype)
comb$celltype_cc <- md$celltype_cc
# cycling-CM as % of ALL cells, and as % of CMs, KO vs WT per timepoint
cm_summary <- md |> filter(celltype == "Cardiomyocyte") |>
  group_by(timepoint, genotype) |>
  summarise(n_CM = dplyr::n(), pct_CM_cycling = round(100*mean(cycling),1),
            cyclingCM_pct_of_all = round(100*sum(cycling)/sum(md$timepoint==timepoint[1]),2), .groups="drop")
write.csv(cm_summary, file.path(OUTTAB, "cellcycle_cyclingCM_summary.csv"), row.names = FALSE)

prop <- tryCatch({
  md$sample_lane <- paste(md$orig.ident, md$lane, sep = "_")
  res <- do.call(rbind, lapply(TIMEPOINTS, function(tp) {
    k <- md$timepoint == tp
    r <- speckle::propeller(clusters = md$celltype_cc[k], sample = md$sample_lane[k], group = md$genotype[k])
    r$timepoint <- tp; r$celltype <- rownames(r); rownames(r) <- NULL; r
  }))
  res$NOTE <- "descriptive: propeller 'samples' are technical lanes (n=1 biological)"
  write.csv(res, file.path(OUTTAB, "cellcycle_propeller_CMsplit.csv"), row.names = FALSE)
  res[grepl("Cardiomyocyte", res$celltype), ]
}, error = function(e) { message("propeller on CM-split failed: ", conditionMessage(e)); NULL })

cat("\n--- cycling fraction by cell type ---\n"); print(as.data.frame(cyc), row.names = FALSE)
cat("\n--- cycling cardiomyocytes KO vs WT ---\n"); print(cm_summary, row.names = FALSE)
if (!is.null(prop)) { cat("\n--- propeller: cycling vs non-cycling CM compartment ---\n")
  print(prop[, intersect(c("timepoint","celltype","PropMean.WT","PropMean.KO","PropRatio","P.Value","FDR"), names(prop))], row.names = FALSE) }
cat("\n=== DONE cell_cycle ===\n")
