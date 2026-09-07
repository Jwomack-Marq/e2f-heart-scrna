#!/usr/bin/env Rscript
# Cardiomyocyte CYCLING-CLUSTER investigation (advisor follow-up, 2026-06).
# Question raised: at res 0.2 / 0.3 the two "cycling" CM subclusters (e.g. CM4 & CM5,
# or CM4 & CM6) look strange -- the identity-marker heatmap shows their markers
# OVERLAP, and one of them appears to MIX S- and G2/M-phase CMs. Are these two real,
# distinct subclusters, or one proliferating population that low-resolution clustering
# split (and that is better explained by phase and/or timepoint than by identity)?
#
# For EACH requested resolution this writes, namespaced by resX.Y:
#   1. cm_cycling_resX.Y_phase_by_split.csv -- subcluster x Phase x timepoint x genotype
#      counts + within-(subcluster,timepoint,genotype) Phase fractions. Tests "is the
#      first cluster mixing S and G2/M" and "does the split track P0/P7 or WT/KO".
#   2. cm_cycling_resX.Y_subcluster_scores.csv -- per-subcluster mean S.Score / G2M.Score
#      / cycling-marker score / size + a cycling flag, to pick the cycling subclusters.
#   3. cm_cycling_resX.Y_marker_jaccard.csv -- Jaccard overlap of top one-vs-rest markers
#      between every subcluster pair (quantifies the heatmap "smearing").
#   4. figures: phase stacked bars faceted timepoint x genotype, and an S.Score vs
#      G2M.Score scatter coloured by subcluster (a single cluster spanning both arms =
#      one cycling population caught mid S->G2M, not two identities).
#
#   Rscript cm_cycling_investigate.R            # default resolutions 0.2 0.3
#   Rscript cm_cycling_investigate.R 0.2 0.3 0.4
# Reads processing/seurat.cm.subclustered.rds; writes results/{tables,figures}/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(ggplot2); library(Matrix) }))

.cli <- commandArgs(trailingOnly = TRUE)
RESOLUTIONS <- if (length(.cli) >= 1 && all(nzchar(.cli))) as.numeric(.cli) else c(0.2, 0.3)
TOPN_MARK   <- 25     # top one-vs-rest markers per subcluster for the Jaccard overlap
PADJ_MARK   <- 0.05   # marker significance cutoff (descriptive; ranking aid only)
CYCLE_MARKERS <- c("Mki67","Top2a","Ccnb1","Ccnb2","Aurkb","Cdca8","Cdk1","Birc5","Cenpa","Bub1")

getdata <- function(o, assay, what) tryCatch(
  SeuratObject::GetAssayData(o, assay = assay, slot = what),
  error = function(e) SeuratObject::GetAssayData(o, assay = assay, layer = what))

say <- function(...) cat(sprintf("[cm_cycling] %s\n", paste0(...)))

CM <- file.path(PROC, "seurat.cm.subclustered.rds")
stopifnot(file.exists(CM))
say("loading ", CM, " ...")
obj <- readRDS(CM)
md  <- obj@meta.data
say(sprintf("CM object: %d cells", ncol(obj)))

# resolve genotype / timepoint / phase columns defensively
gcol <- if ("genotype"  %in% names(md)) "genotype"  else { md$genotype  <- genotype_of(md$orig.ident); "genotype" }
tcol <- if ("timepoint" %in% names(md)) "timepoint" else { md$timepoint <- ifelse(grepl("P0", md$orig.ident), "P0", "P7"); "timepoint" }
if (!"Phase" %in% names(md)) stop("No 'Phase' column on the CM object -- run CellCycleScoring first.")
has_scores <- all(c("S.Score","G2M.Score") %in% names(md))
cyc_genes  <- intersect(CYCLE_MARKERS, rownames(obj))
say(sprintf("cycling markers present: %s", paste(cyc_genes, collapse = ", ")))
# per-cell cycling-marker score = mean log-norm SCT expression of the cycling set
sct <- getdata(obj, "SCT", "data")
md$cyc_score <- if (length(cyc_genes)) Matrix::colMeans(sct[cyc_genes, , drop = FALSE]) else NA_real_

# top one-vs-rest markers per subcluster (presto if available, else Seurat)
top_markers_by_cluster <- function(ident) {
  Idents(obj) <- ident
  if (requireNamespace("presto", quietly = TRUE)) {
    w <- presto::wilcoxauc(sct, as.character(ident))
    w <- w[w$padj < PADJ_MARK & w$logFC > 0, ]
    split(w$feature[order(-w$logFC)], w$group[order(-w$logFC)]) |>
      lapply(function(g) head(g, TOPN_MARK))
  } else {
    m <- Seurat::FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25, verbose = FALSE)
    m <- m[m$p_val_adj < PADJ_MARK, ]
    split(m$gene[order(-m$avg_log2FC)], m$cluster[order(-m$avg_log2FC)]) |>
      lapply(function(g) head(g, TOPN_MARK))
  }
}
jaccard <- function(a, b) if (!length(a) && !length(b)) NA_real_ else
  length(intersect(a, b)) / length(union(a, b))

for (r in RESOLUTIONS) {
  rc <- sprintf("SCT_snn_res.%s", r)
  if (!rc %in% names(md)) { say(sprintf("res %s (%s) not on object -- skipping", r, rc)); next }
  tag <- sub("\\.", ".", sprintf("res%s", r))
  say(sprintf("=== resolution %s (%s) ===", r, rc))
  sub <- factor(paste0("CM", md[[rc]]))
  sub <- factor(sub, levels = paste0("CM", sort(as.integer(sub("CM", "", levels(sub))))))
  d <- data.frame(sub = sub, Phase = factor(md$Phase, levels = c("G1","S","G2M")),
                  genotype = md[[gcol]], timepoint = md[[tcol]], cyc_score = md$cyc_score,
                  S.Score = if (has_scores) md$S.Score else NA_real_,
                  G2M.Score = if (has_scores) md$G2M.Score else NA_real_)

  # (2) per-subcluster summary -> identify the cycling subclusters --------------
  agg <- do.call(rbind, lapply(split(d, d$sub), function(s) data.frame(
    sub = s$sub[1], n = nrow(s),
    pct_S   = 100 * mean(s$Phase == "S",   na.rm = TRUE),
    pct_G2M = 100 * mean(s$Phase == "G2M", na.rm = TRUE),
    pct_G1  = 100 * mean(s$Phase == "G1",  na.rm = TRUE),
    mean_S.Score   = mean(s$S.Score,   na.rm = TRUE),
    mean_G2M.Score = mean(s$G2M.Score, na.rm = TRUE),
    mean_cyc_score = mean(s$cyc_score, na.rm = TRUE))))
  agg$cycling_fraction <- (agg$pct_S + agg$pct_G2M) / 100
  agg$is_cycling <- agg$cycling_fraction >= 0.5 | agg$mean_cyc_score >= median(agg$mean_cyc_score, na.rm = TRUE) + 0.5
  write.csv(agg, file.path(OUTTAB, sprintf("cm_cycling_%s_subcluster_scores.csv", tag)), row.names = FALSE)
  cyc_subs <- as.character(agg$sub[agg$is_cycling])
  say(sprintf("cycling subclusters: %s", paste(cyc_subs, collapse = ", ")))

  # (1) phase composition split by timepoint x genotype ------------------------
  tab <- as.data.frame(table(d$sub, d$Phase, d$timepoint, d$genotype))
  names(tab) <- c("sub","Phase","timepoint","genotype","n")
  tab <- do.call(rbind, lapply(split(tab, list(tab$sub, tab$timepoint, tab$genotype), drop = TRUE),
    function(s) { s$frac <- if (sum(s$n)) s$n / sum(s$n) else 0; s }))
  write.csv(tab, file.path(OUTTAB, sprintf("cm_cycling_%s_phase_by_split.csv", tag)), row.names = FALSE)

  # (3) marker overlap (Jaccard) between subclusters ---------------------------
  mk <- top_markers_by_cluster(d$sub)
  lev <- levels(d$sub); lev <- lev[lev %in% names(mk)]
  jm <- outer(lev, lev, Vectorize(function(a, b) jaccard(mk[[a]], mk[[b]])))
  dimnames(jm) <- list(lev, lev)
  write.csv(data.frame(subcluster = lev, round(jm, 3), check.names = FALSE),
            file.path(OUTTAB, sprintf("cm_cycling_%s_marker_jaccard.csv", tag)), row.names = FALSE)
  if (length(cyc_subs) >= 2) {
    prs <- combn(intersect(cyc_subs, lev), 2)
    for (i in seq_len(ncol(prs)))
      say(sprintf("  marker Jaccard %s vs %s = %.2f (shared: %s)",
        prs[1,i], prs[2,i], jm[prs[1,i], prs[2,i]],
        paste(head(intersect(mk[[prs[1,i]]], mk[[prs[2,i]]]), 12), collapse = ", ")))
  }

  # (4a) phase stacked bars, faceted timepoint x genotype ----------------------
  ggplot(tab, aes(sub, frac, fill = Phase)) + geom_col() +
    facet_grid(genotype ~ timepoint) +
    scale_fill_manual(values = c(G1 = "#bdbdbd", S = "#1565c0", G2M = "#c62828"), na.value = "grey90") +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "subcluster", y = "fraction of cells", title = sprintf("CM phase by subcluster -- res %s", r))
  ggsave(file.path(OUTFIG, sprintf("cm_cycling_%s_phase_by_split.png", tag)), width = 9, height = 7, dpi = 150)

  # (4b) S vs G2M score scatter coloured by subcluster (cycling subs emphasised)
  if (has_scores) {
    d2 <- d; d2$hl <- ifelse(as.character(d2$sub) %in% cyc_subs, as.character(d2$sub), "other")
    ggplot(d2, aes(S.Score, G2M.Score, color = hl)) +
      geom_point(size = .5, alpha = .5) +
      geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
      theme_minimal(base_size = 12) +
      labs(color = "subcluster", title = sprintf("Cell-cycle scores by CM subcluster -- res %s", r),
           subtitle = "a single cluster spanning both arms = one cycling pool caught mid S->G2M")
    ggsave(file.path(OUTFIG, sprintf("cm_cycling_%s_score_scatter.png", tag)), width = 8, height = 6.5, dpi = 150)
  }
}

say("=== DONE cm_cycling_investigate ===")
say("Read scores + jaccard CSVs together: if the two cycling subclusters have high")
say("mutual marker Jaccard AND one carries both S and G2/M cells, they are one")
say("proliferating population split by phase/timepoint, not two distinct identities.")
