#!/usr/bin/env Rscript
# Marker-QC plots for parity with their 00b/tierA_plots, from our existing
# results/markers/*.allmarkers.csv: pct.1-vs-pct.2 specificity scatter,
# n-significant-markers per cluster, marker-overlap Jaccard, top-5 marker heatmaps.
# Cosmetic completeness (the underlying data already exists).
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(ggplot2); library(dplyr) }))
MK <- file.path(RESULTS, "markers")
conds <- c("P0WT","P0KO","P7WT","P7KO")

load1 <- function(cc) {
  f <- file.path(MK, paste0(cc, ".allmarkers.csv")); if (!file.exists(f)) return(NULL)
  d <- read.csv(f); names(d)[names(d) == "" | is.na(names(d))] <- "geneid"; d$condition <- cc; d
}
all <- do.call(rbind, lapply(conds, load1))
all <- all[is.finite(all$avg_log2FC), ]
all$sig <- all$p_val_adj < 0.05 & all$avg_log2FC > 0.25

# 1) specificity scatter (pct.1 vs pct.2)
png(file.path(OUTFIG, "markerqc_specificity_scatter.png"), 1000, 800)
print(ggplot(subset(all, sig), aes(pct.2, pct.1)) + geom_point(size = .3, alpha = .3) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        facet_wrap(~ condition) + theme_bw() +
        labs(title = "Marker specificity (pct.1 vs pct.2; above diagonal = specific)", x = "pct.2 (other clusters)", y = "pct.1 (in cluster)"))
dev.off()

# 2) n significant markers per cluster
nsig <- all |> filter(sig) |> count(condition, cluster)
png(file.path(OUTFIG, "markerqc_n_sig_markers.png"), 1100, 600)
print(ggplot(nsig, aes(factor(cluster), n)) + geom_col(fill = "steelblue") +
        facet_wrap(~ condition, scales = "free_x") + theme_bw() +
        labs(title = "Significant markers per cluster (p_adj<0.05, log2FC>0.25)", x = "cluster", y = "# markers"))
dev.off()

# 3) marker-overlap Jaccard between conditions (significant gene sets)
sets <- lapply(conds, function(cc) unique(all$gene[all$condition == cc & all$sig]))
names(sets) <- conds
J <- outer(seq_along(sets), seq_along(sets), Vectorize(function(i,j)
  round(length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]])), 2)))
dimnames(J) <- list(conds, conds)
Jl <- as.data.frame(as.table(J)); names(Jl) <- c("A","B","jaccard")
png(file.path(OUTFIG, "markerqc_overlap_jaccard.png"), 650, 550)
print(ggplot(Jl, aes(A, B, fill = jaccard, label = jaccard)) + geom_tile() + geom_text() +
        scale_fill_gradient(low = "white", high = "firebrick") + theme_bw() +
        labs(title = "Marker-set overlap (Jaccard) between conditions", x = NULL, y = NULL))
dev.off()

# 4) top-5 markers/cluster heatmap (avg_log2FC), one panel per condition
for (cc in conds) {
  d <- subset(all, condition == cc & sig)
  if (!nrow(d)) next
  top <- d |> group_by(cluster) |> slice_max(avg_log2FC, n = 5, with_ties = FALSE) |> ungroup()
  png(file.path(OUTFIG, paste0("markerqc_top5_heatmap_", cc, ".png")), 800, 900)
  print(ggplot(top, aes(factor(cluster), gene, fill = avg_log2FC)) + geom_tile() +
          scale_fill_gradient(low = "white", high = "darkred") + theme_bw(base_size = 7) +
          labs(title = paste0(cc, " top-5 markers/cluster"), x = "cluster", y = NULL))
  dev.off()
}
cat("marker-QC plots written (specificity, n-sig, jaccard, top5 x", length(conds), ")\n")
cat("=== DONE marker_qc_plots ===\n")
