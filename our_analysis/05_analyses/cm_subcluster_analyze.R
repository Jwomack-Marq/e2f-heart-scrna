#!/usr/bin/env Rscript
# Cardiomyocyte SUBCLUSTERING -- analyze step (cheap; re-run per resolution).
# Loads seurat.cm.subclustered.rds and, for EACH resolution, writes:
#   (3) DISTINGUISHED marker heatmap -- top markers/subcluster, genes on y, subcluster
#       blocks labelled on the RIGHT. Clean diagonal blocks => distinct subclusters;
#       smeared/overlapping blocks => lower resolution (fewer subclusters needed).
#   (2) KO-vs-WT DESCRIPTIVE DE within EACH subcluster (pseudobulk, ~timepoint+condition,
#       apeglm-shrunken LFC). n=1, sex-confounded => NO valid p-values; thin subclusters
#       skipped; sex/ROSA26 confounder genes flagged.
#   (4) G2/M/S phase across subclusters: phase composition + a cell-cycle marker heatmap.
# All outputs are namespaced by resolution (resX.Y), so resolutions coexist on disk.
#   Default run does BOTH kept resolutions 0.1 and 0.2; override with CLI args, e.g.
#   `Rscript cm_subcluster_analyze.R 0.3`  or  `Rscript cm_subcluster_analyze.R 0.1 0.2 0.3`.
# Reads processing/seurat.cm.subclustered.rds; writes results/{figures,tables}/.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(DESeq2); library(ggplot2); library(dplyr); library(tidyr) }))

.cli <- commandArgs(trailingOnly = TRUE)
RESOLUTIONS <- if (length(.cli) >= 1 && all(nzchar(.cli))) as.numeric(.cli) else c(0.1, 0.2)
TOPN             <- 10           # top markers per subcluster for the heatmap
MIN_CELLS_PER_PB <- 20           # min cells per (group,lane) pseudobulk sample (mirrors per_celltype_de.R)
CONFOUND <- c("Eif2s3y","Kdm5d","Uty","Ddx3y","Xist","Tsix","Gt(ROSA)26Sor")  # sex + ROSA26 construct
# Known CM-subtype vocabulary (from cm_subtypes.R) to annotate data-driven subclusters
CM_SUBTYPES <- list(Ventricular=c("Myl2","Myh7"), Atrial=c("Myl7","Sln","Nppa"),
                    Trabecular=c("Bmp10","Nppa","Hey2"), Compact=c("Hey2","Irx3","Tbx20"),
                    Cycling=c("Mki67","Top2a","Ccnb1","Aurkb","Cdca8"))

# cross-version GetAssayData (slot in v4, layer in v5; slot still accepted in v5)
getdata <- function(o, assay, what) tryCatch(
  SeuratObject::GetAssayData(o, assay = assay, slot = what),
  error = function(e) SeuratObject::GetAssayData(o, assay = assay, layer = what))

# average log-norm SCT expression per subcluster, z-scored across subclusters, capped
avg_z <- function(obj, genes, ident, cap = 2.5) {
  genes <- intersect(genes, rownames(obj))
  dat <- getdata(obj, "SCT", "data")[genes, , drop = FALSE]
  avg <- sapply(levels(ident), function(cl) Matrix::rowMeans(dat[, ident == cl, drop = FALSE]))
  m <- t(scale(t(as.matrix(avg)))); m[is.na(m)] <- 0
  m[m > cap] <- cap; m[m < -cap] <- -cap
  list(mat = m, genes = genes)
}

# marker-block heatmap: rows=genes (split by the block they mark), block labels on RIGHT.
# ComplexHeatmap preferred -> pheatmap -> ggplot(facet_grid: right-side strips) fallback.
block_heatmap <- function(mat, block, file, title, w = 850, h = 1000, legend = "z-score") {
  block <- droplevels(factor(block, levels = levels(block)))
  pal <- setNames(grDevices::hcl.colors(nlevels(block), "Dark 3"), levels(block))
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)) {
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("steelblue", "white", "firebrick"))
    ha <- ComplexHeatmap::rowAnnotation(block = block, col = list(block = pal),
            show_legend = FALSE, show_annotation_name = FALSE)
    ht <- ComplexHeatmap::Heatmap(mat, name = legend, col = col_fun,
            cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
            row_split = block, row_title_side = "right", row_title_rot = 0,
            right_annotation = ha, row_names_gp = grid::gpar(fontsize = 7),
            column_names_gp = grid::gpar(fontsize = 10), column_names_rot = 45,
            column_title = title)
    png(file, w, h); ComplexHeatmap::draw(ht, heatmap_legend_side = "left"); dev.off()
    return("ComplexHeatmap")
  }
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann <- data.frame(block = block); rownames(ann) <- rownames(mat)
    gaps <- head(cumsum(as.integer(table(block))), -1)
    pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
      gaps_row = gaps, fontsize_row = 7, main = title, annotation_colors = list(block = pal),
      color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
      filename = file, width = w/100, height = h/100)
    return("pheatmap")
  }
  blk <- setNames(as.character(block), rownames(mat))
  df <- as.data.frame(as.table(mat)); names(df) <- c("gene", "subcluster", "z")
  df$gene <- factor(df$gene, levels = rev(rownames(mat)))
  df$block <- factor(blk[as.character(df$gene)], levels = levels(block))
  p <- ggplot(df, aes(subcluster, gene, fill = z)) + geom_tile() +
    facet_grid(block ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y.right = element_text(angle = 0)) + labs(title = title, x = NULL, y = NULL)
  ggsave(file, p, width = w/100, height = h/100, dpi = 120, limitsize = FALSE)
  "ggplot"
}

## --- load object ONCE (run_res works on a local copy per resolution) --------
cm0 <- readRDS(file.path(PROC, "seurat.cm.subclustered.rds"))
cm0$genotype <- genotype_of(cm0$orig.ident)

run_res <- function(RES, cm) {
  res_col <- paste0("SCT_snn_res.", RES)
  if (!res_col %in% colnames(cm@meta.data)) { message("resolution ", RES, " (", res_col, ") not in object -- skipped."); return(invisible()) }
  ids <- as.character(cm[[res_col]][,1])
  lev <- paste0("CM", sort(as.integer(unique(ids))))
  cm$cm_subcluster <- factor(paste0("CM", ids), levels = lev)
  Idents(cm) <- "cm_subcluster"
  cat(sprintf("\n########## RESOLUTION %.1f -> %d subclusters: %s ##########\n",
              RES, length(lev), paste(lev, collapse = ", ")))
  print(table(cm$cm_subcluster, cm$orig.ident))

  ## === (3) markers + distinguished heatmap =================================
  DefaultAssay(cm) <- "SCT"
  cm <- PrepSCTFindMarkers(cm, verbose = FALSE)
  markers <- FindAllMarkers(cm, assay = "SCT", only.pos = TRUE,
                            min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  markers$cluster <- factor(markers$cluster, levels = lev)
  write.csv(markers, file.path(OUTTAB, sprintf("cm_subcluster_markers_res%s.csv", RES)), row.names = FALSE)

  top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = TOPN, with_ties = FALSE) %>%
    ungroup() %>% group_by(gene) %>% slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>% ungroup() %>%
    arrange(cluster, desc(avg_log2FC))
  write.csv(top, file.path(OUTTAB, sprintf("cm_subcluster_top_markers_res%s.csv", RES)), row.names = FALSE)

  az <- avg_z(cm, top$gene, cm$cm_subcluster)
  blk <- droplevels(factor(top$cluster[match(az$genes, top$gene)], levels = lev))
  hbk <- block_heatmap(az$mat[az$genes, , drop = FALSE], blk,
    file.path(OUTFIG, sprintf("cm_subcluster_marker_heatmap_res%s.png", RES)),
    sprintf("CM subcluster markers (res %.1f) -- top %d/subcluster", RES, TOPN),
    w = 850, h = max(700, 22 * length(az$genes)))
  cat(sprintf("marker heatmap drawn via %s\n", hbk))

  # annotate data-driven subclusters with the known CM-subtype vocabulary (descriptive)
  DefaultAssay(cm) <- "RNA"; cm <- NormalizeData(cm, verbose = FALSE)
  ms.cols <- c()
  for (s in names(CM_SUBTYPES)) {
    mk <- intersect(CM_SUBTYPES[[s]], rownames(cm)); if (length(mk) < 1) next
    cm <- AddModuleScore(cm, features = list(mk), name = paste0("cms_", s)); ms.cols[s] <- paste0("cms_", s, "1")
  }
  cl.mean <- sapply(ms.cols, function(col) tapply(cm[[col]][,1], cm$cm_subcluster, mean))
  sub.lab <- colnames(cl.mean)[max.col(cl.mean, ties.method = "first")]; names(sub.lab) <- rownames(cl.mean)
  write.csv(data.frame(subcluster = names(sub.lab), nearest_CM_subtype = unname(sub.lab)),
            file.path(OUTTAB, sprintf("cm_subcluster_subtype_map_res%s.csv", RES)), row.names = FALSE)

  ## === (2) KO-vs-WT descriptive DE within each subcluster ==================
  # NB: AggregateExpression replaces '_' in group names with '-', so we derive the
  # per-sample genotype/timepoint from colnames(pb) (orig.ident is '-'/'_'-free)
  # instead of indexing by the raw pb.sample names (which would mismatch).
  cm$pb.sample <- paste(cm$orig.ident, cm$lane, sep = "_")
  summ <- list()
  for (scl in lev) {
    sub <- cm[, cm$cm_subcluster == scl]
    pb <- AggregateExpression(sub, assays = "RNA", group.by = "pb.sample", slot = "counts")$RNA
    samp <- colnames(pb)                                          # e.g. "P0KO-lane1"
    pbc <- table(sub$pb.sample); names(pbc) <- gsub("_", "-", names(pbc))
    ncell <- as.integer(pbc[samp])                                # cells per pseudobulk sample
    oid   <- sub("-lane[0-9]+$", "", samp)                        # orig.ident, e.g. "P0KO"
    gtype <- setNames(ifelse(grepl("KO$", oid), "KO", "WT"), samp)
    tpv   <- setNames(substr(oid, 1, 2), samp)
    keep  <- samp[ncell >= MIN_CELLS_PER_PB]
    nKO <- sum(gtype[keep] == "KO"); nWT <- sum(gtype[keep] == "WT")
    rowbase <- data.frame(subcluster = scl, n_cells = ncol(sub), n_KO_samp = nKO, n_WT_samp = nWT)
    if (nKO < 2 || nWT < 2) {
      summ[[scl]] <- cbind(rowbase, design = NA, status = "skipped_too_few_or_unbalanced",
                           n_DE_absLFC_gt1 = NA, top_KO_up = NA, top_KO_down = NA)
      cat(sprintf("  %-5s n=%5d  KO_samp=%d WT_samp=%d  -> SKIP\n", scl, ncol(sub), nKO, nWT)); next
    }
    pb <- pb[, keep, drop = FALSE]
    cd <- data.frame(row.names = keep,
                     condition = relevel(factor(gtype[keep]), ref = "WT"),
                     timepoint = factor(tpv[keep]))
    design <- if (nlevels(cd$timepoint) > 1) ~ timepoint + condition else ~ condition
    dds <- DESeqDataSetFromMatrix(round(as.matrix(pb)), colData = cd, design = design)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) { message("    DESeq failed ", scl, ": ", conditionMessage(e)); NULL })
    if (is.null(dds)) {
      summ[[scl]] <- cbind(rowbase, design = paste(deparse(design), collapse=""), status = "DESeq_failed",
                           n_DE_absLFC_gt1 = NA, top_KO_up = NA, top_KO_down = NA); next
    }
    coef <- grep("condition_KO_vs_WT", resultsNames(dds), value = TRUE)
    res <- as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm"))
    res$gene <- rownames(res); res$confounder <- res$gene %in% CONFOUND
    res$NOTE <- "descriptive_n1_sexconfound_no_valid_pvalues"
    res <- res[order(-abs(res$log2FoldChange)), ]
    write.csv(res, file.path(OUTTAB, sprintf("cm_subcluster_res%s_KOvsWT_%s.descriptive.DE.csv", RES, scl)), row.names = FALSE)
    bio <- res[!res$confounder, ]
    up <- head(bio$gene[bio$log2FoldChange > 1], 6); dn <- head(bio$gene[bio$log2FoldChange < -1], 6)
    summ[[scl]] <- cbind(rowbase, design = paste(deparse(design), collapse=""), status = "ok",
      n_DE_absLFC_gt1 = sum(abs(bio$log2FoldChange) > 1, na.rm = TRUE),
      top_KO_up = paste(up, collapse = ","), top_KO_down = paste(dn, collapse = ","))
    cat(sprintf("  %-5s n=%5d  KO-up: %s\n", scl, ncol(sub), paste(up, collapse = ", ")))
  }
  summ <- do.call(rbind, summ)
  write.csv(summ, file.path(OUTTAB, sprintf("cm_subcluster_res%s_KOvsWT_summary.csv", RES)), row.names = FALSE)

  ## === (4) G2/M/S phase across subclusters =================================
  if (!"Phase" %in% colnames(cm@meta.data)) {
    cc <- cc_lists()
    cm <- CellCycleScoring(cm, s.features = intersect(cc$S, rownames(cm)),
                           g2m.features = intersect(cc$G2M, rownames(cm)), set.ident = FALSE)
  }
  cm$cycling <- cm$Phase %in% c("S", "G2M")
  lv <- c("G1", "S", "G2M"); lv <- c(lv, setdiff(unique(as.character(cm$Phase)), lv))
  cm$Phase <- factor(cm$Phase, levels = lv)

  cyc <- cm@meta.data |> dplyr::group_by(cm_subcluster, genotype) |>
    dplyr::summarise(n = dplyr::n(), pct_S = round(100*mean(Phase == "S"), 1),
              pct_G2M = round(100*mean(Phase == "G2M"), 1),
              pct_cycling = round(100*mean(cycling), 1), .groups = "drop")
  write.csv(cyc, file.path(OUTTAB, sprintf("cm_subcluster_res%s_cellcycle.csv", RES)), row.names = FALSE)

  # NB: dplyr::count is masked by plyr (pulled in via DESeq2's deps) -> use group_by/summarise.
  ph <- cm@meta.data |> dplyr::group_by(cm_subcluster, genotype, Phase) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") |>
    dplyr::mutate(frac = n / sum(n)) |> dplyr::ungroup()
  p_ph <- ggplot(ph, aes(cm_subcluster, frac, fill = Phase)) + geom_col() + facet_wrap(~ genotype) +
    theme_bw() + labs(title = sprintf("Cell-cycle phase composition by CM subcluster (res %.1f, descriptive n=1)", RES),
                      x = NULL, y = "fraction") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUTFIG, sprintf("cm_subcluster_res%s_cellcycle_phase_composition.png", RES)), p_ph, width = 10, height = 5, dpi = 120)

  cc <- cc_lists(); sgen <- intersect(cc$S, rownames(cm)); ggen <- setdiff(intersect(cc$G2M, rownames(cm)), sgen)
  DefaultAssay(cm) <- "SCT"
  azcc <- avg_z(cm, c(sgen, ggen), cm$cm_subcluster)
  blkcc <- factor(ifelse(azcc$genes %in% sgen, "S", "G2M"), levels = c("S", "G2M"))
  hcc <- block_heatmap(azcc$mat[azcc$genes, , drop = FALSE], blkcc,
    file.path(OUTFIG, sprintf("cm_subcluster_res%s_cellcycle_marker_heatmap.png", RES)),
    sprintf("S / G2M phase markers by CM subcluster (res %.1f)", RES), w = 800, h = max(700, 16 * length(azcc$genes)))
  cat(sprintf("cell-cycle marker heatmap drawn via %s\n", hcc))

  ## --- subcluster + phase UMAP ----------------------------------------------
  png(file.path(OUTFIG, sprintf("cm_subcluster_res%s_umap.png", RES)), 1400, 600)
  print((DimPlot(cm, group.by = "cm_subcluster", reduction = "umap", label = TRUE) + ggtitle(sprintf("CM subclusters (res %.1f)", RES))) |
        (DimPlot(cm, group.by = "Phase", reduction = "umap") + ggtitle("Cell-cycle phase")))
  dev.off()

  cat("\n--- KO-vs-WT per subcluster (biological genes; sex/construct flagged) ---\n")
  print(summ[, c("subcluster","n_cells","status","n_DE_absLFC_gt1","top_KO_up")], row.names = FALSE)
  cat("\n--- cycling fraction by subcluster ---\n"); print(as.data.frame(cyc), row.names = FALSE)
  cat(sprintf("=== DONE resolution %.1f ===\n", RES))
  invisible()
}

for (RES in RESOLUTIONS) { run_res(RES, cm0); gc(verbose = FALSE) }
cat(sprintf("\n=== DONE cm_subcluster_analyze (resolutions: %s) ===\n", paste(RESOLUTIONS, collapse = ", ")))
