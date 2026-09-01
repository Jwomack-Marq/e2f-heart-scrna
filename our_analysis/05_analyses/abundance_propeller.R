#!/usr/bin/env Rscript
# Parity with their 07_abundance.qmd: speckle::propeller test of cell-type
# proportions KO vs WT, per timepoint. NOTE (n=1 caveat): the only within-group
# replicate unit available is the technical sequencing lane (lane1/lane6 = same
# library), so propeller's "samples" here are technical, not biological. propeller
# tolerates small n but remains underpowered -- treat as DESCRIPTIVE.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(speckle); library(ggplot2) }))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
md <- comb@meta.data
md$genotype  <- genotype_of(md$orig.ident)
md$sample_lane <- paste(md$orig.ident, md$lane, sep = "_")   # technical-lane "replicates"

run_prop <- function(tp) {
  k <- md$timepoint == tp
  pr <- speckle::propeller(clusters = md$celltype[k], sample = md$sample_lane[k], group = md$genotype[k])
  pr$timepoint <- tp; pr$celltype <- rownames(pr); rownames(pr) <- NULL
  pr
}
res <- do.call(rbind, lapply(TIMEPOINTS, run_prop))
res$NOTE <- "descriptive: propeller 'samples' are technical lanes (n=1 biological)"
write.csv(res, file.path(OUTTAB, "propeller_results.csv"), row.names = FALSE)

cat("\n--- propeller (cell-type proportion KO vs WT, per timepoint) ---\n")
keep <- intersect(c("timepoint","celltype","Tstatistic","PropMean.WT","PropMean.KO","P.Value","FDR",
                    "PropRatio","Tstatistic"), names(res))
print(res[order(res$timepoint, res$FDR), unique(c("timepoint","celltype", grep("Prop|P.Value|FDR|stat", names(res), value=TRUE, ignore.case=TRUE)))], row.names = FALSE)

# barplot of -log10 FDR per cell type (descriptive)
res$neglogFDR <- -log10(pmax(res$FDR, 1e-300))
p <- ggplot(res, aes(reorder(celltype, neglogFDR), neglogFDR, fill = timepoint)) +
  geom_col(position = "dodge") + coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "propeller KO-vs-WT proportion shift (DESCRIPTIVE; lane=technical rep)",
       x = NULL, y = "-log10(FDR) [underpowered, n=1 biological]") + theme_bw()
ggsave(file.path(OUTFIG, "abundance_propeller.png"), p, width = 900/96, height = 600/96, dpi = 96)
cat("=== DONE abundance_propeller ===\n")
