# 03_excel.R -- the two workbooks the collaborator asked for.
# ---------------------------------------------------------------------------
# Reads de_tables.rds (01) and enrich.rds (02) and writes:
#   P7_KO_vs_WT_by_CM_subcluster.xlsx   part 1
#   P7WT_vs_P0WT.xlsx                   part 2
#
# There is no Excel code anywhere else in this repo, so nothing to reuse here.
# Every sheet is frozen at row 1 with an autofilter; the README sheet carries the
# caveats, because a table handed over as a file loses whatever context the
# website put around it.
# ---------------------------------------------------------------------------

suppressMessages(library(openxlsx))

# Conventions shared with the app (shiny_app/download_helpers.R): the caveat text,
# the 3-significant-figure rule, the column widths and the 31-char sheet-name
# guard. Sourced rather than copied so an emailed workbook and one downloaded
# from the app cannot end up saying different things. The sheet-building below
# stays bespoke -- it does conditional formatting and per-sheet freeze columns
# that the app's generic writer does not need.
HELPERS <- Sys.getenv("DL_HELPERS", "/repo/shiny_app/download_helpers.R")
if (!file.exists(HELPERS)) HELPERS <- "shiny_app/download_helpers.R"
stopifnot(file.exists(HELPERS))
source(HELPERS)

OUT <- "/out"
de  <- readRDS(file.path(OUT, "de_tables.rds"))
en  <- readRDS(file.path(OUT, "enrich.rds"))
ms  <- if (file.exists(file.path(OUT, "mito_sensitivity.rds")))
         readRDS(file.path(OUT, "mito_sensitivity.rds")) else NULL
P   <- de$params
man <- de$manifest

HDR <- createStyle(textDecoration = "bold", fgFill = "#EEEEEE", halign = "left",
                   border = "bottom", borderColour = "#999999")
WRAP <- createStyle(wrapText = TRUE, valign = "top")

add_sheet <- function(wb, name, df, widths = "auto", freeze_col = 1) {
  stopifnot(nchar(name) <= 31)
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = HDR, withFilter = TRUE)
  freezePane(wb, name, firstActiveRow = 2, firstActiveCol = freeze_col + 1)
  # The per-sheet width vectors are written against the expected column set. A
  # clusterProfiler or fgsea version that drops or adds a column would otherwise
  # abort the whole workbook over cosmetics, so pad/trim instead of failing.
  n <- max(1, ncol(df))
  if (length(widths) != 1 && length(widths) != n)
    widths <- if (length(widths) > n) widths[seq_len(n)] else c(widths, rep("auto", n - length(widths)))
  setColWidths(wb, name, cols = seq_len(n), widths = widths)
  invisible(NULL)
}
add_notes <- function(wb, name, lines) {
  addWorksheet(wb, name)
  writeData(wb, name, data.frame(x = lines), colNames = FALSE)
  setColWidths(wb, name, cols = 1, widths = 120)
  addStyle(wb, name, WRAP, rows = seq_along(lines), cols = 1, gridExpand = TRUE)
  invisible(NULL)
}
# colour the direction column so up/down is readable without filtering
colour_direction <- function(wb, name, df) {
  j <- match("direction", names(df)); if (is.na(j)) return(invisible(NULL))
  up <- unique(df$direction[grepl("_up$", df$direction)])
  dn <- setdiff(unique(df$direction), c(up, "ns"))
  for (v in up) conditionalFormatting(wb, name, cols = j, rows = 2:(nrow(df) + 1),
      rule = sprintf('=="%s"', v), type = "expression",
      style = createStyle(bgFill = "#FADBD8", fontColour = "#922B21"))
  for (v in dn) conditionalFormatting(wb, name, cols = j, rows = 2:(nrow(df) + 1),
      rule = sprintf('=="%s"', v), type = "expression",
      style = createStyle(bgFill = "#D6EAF8", fontColour = "#1B4F72"))
}
# from download_helpers.R: explicit widths (openxlsx's "auto" scans every cell,
# which across sheets of ~20,000 rows is minutes of work for a cosmetic result)
# and the 3-significant-figure rule.
de_widths <- dl_widths
sig3      <- dl_signif
go_cols <- c("cluster","stratum","direction","ontology","ID","Description","GeneRatio",
             "BgRatio","FoldEnrichment","pvalue","p.adjust","qvalue","Count","geneID",
             "n_input","n_universe","input_rule")
gsea_cols <- c("cluster","stratum","pathway","NES","pval","padj","size","leadingEdge")

pick <- function(df, cols) df[, intersect(cols, names(df)), drop = FALSE]

# The five shared caveats, plus the two specific to this offline pipeline's
# matrix choices. Kept in this order so the numbering in the README sheet is stable.
CAVEATS <- c(DL_CAVEATS,
  "Differential expression runs on the bundle's genome-wide matrix (24,221 genes), which is downsampled to 8,026 cells. That is the only genome-wide matrix available without the upstream Seurat object, and GO needs a genome-wide gene space. The per-arm cell counts are in the Summary sheet; CM5, CM7 and CM8 are thin (39-96 cells per arm).",
  "The 'lfc_fullcells' / 'padj_fullcells' columns re-run the SAME contrast on the curated 2,181-gene panel, which retains all cells at full depth. Where a gene is on that panel, those columns are the better-powered estimate; they are blank for the other ~22,000 genes.")

# ---------------------------------------------------------------------------
build_workbook <- function(contrast_key, title, extra_notes) {
  m <- man[man$contrast == contrast_key, , drop = FALSE]
  m <- m[order(m$cluster != "AllCM", m$cluster, !m$is_primary), , drop = FALSE]
  go   <- if (is.null(en$go))   NULL else en$go[en$go$contrast == contrast_key, , drop = FALSE]
  gsea <- if (is.null(en$gsea)) NULL else en$gsea[en$gsea$contrast == contrast_key, , drop = FALSE]
  aud  <- en$audit[en$audit$contrast == contrast_key, , drop = FALSE]

  wb <- createWorkbook()
  de_index <- list()

  # ---- README
  up_lab <- m$up_label[1]; down_lab <- m$down_label[1]
  prim <- unique(m$stratum[m$is_primary])
  add_notes(wb, "README", c(
    title, "",
    sprintf("Generated %s from shiny_app/app_data.rds by analysis/2026-08-21_email/ (01_de.R -> 02_enrich.R -> 03_excel.R).", P$built),
    "",
    "WHAT THE CONTRAST IS",
    sprintf("  %s: log2FoldChange > 0 means UP in %s; < 0 means UP in %s.", contrast_key, m$arm_A[1], m$arm_B[1]),
    sprintf("  Significance: padj < %.2g AND |log2FoldChange| >= %.2g  ->  the 'direction' column (%s / %s / ns).",
            P$sig_padj, P$sig_lfc, up_lab, down_lab),
    sprintf("  Primary stratum: %s cells.", prim),
    extra_notes, "",
    "METHOD",
    "  Cell-level Wilcoxon rank-sum via presto::wilcoxauc on the log-normalised matrix -- the same test",
    "  the website's 'Subset & DEGs' tab and its four-group DE grid use, so these numbers are directly",
    "  comparable to what you have already seen there. The difference is that NO row filter is applied:",
    "  every gene in the matrix is reported, which is what makes these the complete lists.",
    sprintf("  Matrix: %s", P$matrix),
    sprintf("  Cross-check matrix: %s", P$crosscheck_matrix),
    sprintf("  GO: clusterProfiler::enrichGO, org.Mm.eg.db, keyType=SYMBOL, BH-adjusted, p<%.2g, q<%.2g, gene-set size 10-500.",
            en$params$go_pcut, en$params$go_qcut),
    sprintf("  GO universe: genes expressed in >= %d%% of at least one arm IN THAT CLUSTER (not all 24,221 genes,", en$params$min_pct_universe),
    "  which would inflate every fold enrichment, and not the significant list).",
    sprintf("  GO is run SEPARATELY on the %s and the %s list, as requested.", up_lab, down_lab),
    sprintf("  GSEA: fgsea over %d MSigDB mouse gene sets (Hallmark + C2 CP:KEGG_LEGACY), ranked by", en$params$n_pathways),
    "  sign(log2FC) * -log10(p) over all non-confounder genes.",
    "",
    "CAVEATS - PLEASE READ BEFORE INTERPRETING", paste0("  ", seq_along(CAVEATS), ". ", CAVEATS),
    "",
    "COLUMN DICTIONARY (DE sheets)",
    "  gene                 mouse gene symbol",
    sprintf("  log2FoldChange       log2 fold change, %s vs %s (>0 = up in %s)", m$arm_A[1], m$arm_B[1], m$arm_A[1]),
    "  auc                  Wilcoxon AUC; 0.5 = no separation, >0.5 = higher in the first arm. Scale-free,",
    "                       so it is the fairer effect size when comparing a high- and a low-dynamic-range gene.",
    "  pvalue / padj        Wilcoxon p and its BH adjustment across all tested genes. Pseudoreplicated - see caveat 1.",
    "  pct_* / mean_*       percent of cells expressing, and mean log-normalised expression, in each arm",
    "  n_*                  cells per arm in the genome-wide (downsampled) matrix",
    "  n_*_fullcells        cells per arm in the full-depth curated-panel cross-check",
    "  confounder           TRUE for the sex/construct genes listed in caveat 2",
    "  sig_sets             membership in the app's curated programs (maturation, mat_mature, mat_immature,",
    "                       glycolysis, faox, prolif, cytokinesis, ccexit, and the Seurat S / G2M phase lists)",
    "  lfc_fullcells        same contrast on the curated panel at full cell depth; blank if off-panel",
    "  lfc_pooled_website   for reference: the value the website's 'KO-vs-WT DE (per subgroup)' tab shows,",
    "                       which is KO vs WT POOLED ACROSS P0 AND P7 - not this P7-specific contrast",
    "  direction            up / down / ns, per the significance rule above",
    "  mito_encoded         TRUE for mitochondrially-encoded genes - see the Mito_check_READ_ME sheet",
    "",
    "SHEETS",
    "  Summary              one row per cluster: arm sizes, gene counts, top genes, top GO term each direction",
    "  DE_<cluster>...      complete gene lists (every tested gene), sorted by padj then |log2FC|",
    "  GO_BP_<cluster>      GO Biological Process, with a 'direction' column separating the two analyses",
    "  GO_MF_CC             GO Molecular Function and Cellular Component, all clusters",
    "  GSEA                 Hallmark + KEGG, all clusters",
    "  GO_audit             per cluster x direction x ontology: input size, universe size, terms found, and the",
    "                       selection rule used. Check here before reading an empty GO sheet as 'nothing enriched'.",
    "  Mito_check_READ_ME   why the mitochondrially-encoded genes need care before the KO-up result is read",
    "  Mito_GO_comparison   GO terms found with vs without the mt- genes, per cluster and direction",
    "  Mito_fraction        measured mitochondrial share of signal per cluster x group",
    "  GO_BP_no_mt          GO BP re-run with mt- genes removed from both the input list and the universe",
    "  Secondary_tables     which tables are in this workbook and which ship as CSV only, with the path",
    "  Skipped              cluster/arm combinations that fell below the 10-cell floor",
    "",
    "WHAT IS IN THE SHEETS vs WHAT IS IN THE CSVs",
    "  The DE sheets hold the primary stratum (plus both strata for AllCM) and drop genes detected in zero",
    "  cells on BOTH sides, purely to keep the file emailable. The csv/ folder next to this workbook holds",
    "  every contrast, every stratum, every gene, unfiltered. Secondary_tables maps one to the other."))

  # ---- Summary
  summ <- do.call(rbind, lapply(seq_len(nrow(m)), function(i) {
    k <- m$key[i]; d <- de$tables[[k]]
    nc <- !d$confounder
    tp <- function(sel, n = 10) paste(head(d$gene[sel & nc], n), collapse = ", ")
    gg <- if (is.null(go)) NULL else go[go$cluster == m$cluster[i] & go$stratum == m$stratum[i] & go$ontology == "BP", ]
    topterm <- function(dirn) {
      if (is.null(gg)) return("")
      x <- gg[gg$direction == dirn, ]
      if (!nrow(x)) "" else x$Description[which.min(x$p.adjust)]
    }
    dd <- d[order(-abs(d$log2FoldChange)), ]
    ddn <- !dd$confounder
    data.frame(cluster = m$cluster[i], stratum = m$stratum[i], primary = m$is_primary[i],
      n_A = m$n_A[i], n_B = m$n_B[i], arm_A = m$arm_A[i], arm_B = m$arm_B[i],
      n_A_fullcells = m$n_A_fullcells[i], n_B_fullcells = m$n_B_fullcells[i],
      n_genes_tested = m$n_genes_tested[i], n_up = m$n_up[i], n_down = m$n_down[i],
      top_up_by_padj = tp(d$direction == m$up_label[i]),
      top_down_by_padj = tp(d$direction == m$down_label[i]),
      top_up_by_lfc = paste(head(dd$gene[dd$direction == m$up_label[i] & ddn], 10), collapse = ", "),
      top_down_by_lfc = paste(head(dd$gene[dd$direction == m$down_label[i] & ddn], 10), collapse = ", "),
      top_GO_BP_up = topterm(m$up_label[i]), top_GO_BP_down = topterm(m$down_label[i]),
      stringsAsFactors = FALSE)
  }))
  add_sheet(wb, "Summary", summ, widths = c(rep("auto", 12), rep(40, 6)), freeze_col = 2)

  # ---- DE sheets
  # A workbook with every stratum at 24,221 rows a sheet runs to ~7M cells, which is
  # slow to write and too big to email. The primary stratum answers the question and
  # goes in the book; AllCM keeps both strata because in part 2 the gap between them
  # IS the sort artefact. Undetected genes (zero cells expressing on both sides) are
  # dropped from the sheets only -- the CSVs under csv/ stay complete, every gene,
  # every stratum, and the Secondary_tables sheet names them.
  in_book <- m$is_primary | m$cluster == "AllCM"
  for (i in which(in_book)) {
    k <- m$key[i]
    nm <- sprintf("DE_%s_%s", m$cluster[i], m$stratum[i])
    d <- de$tables[[k]]
    ndrop <- sum(!d$detected)
    d <- d[d$detected, setdiff(names(d), "detected"), drop = FALSE]
    # mt- genes are up in KO in every cluster and down in none; flag them so a
    # reader can see at a glance which rows carry that block (see the Mito_check sheet).
    d$mito_encoded <- startsWith(d$gene, "mt-")
    d <- sig3(d)
    add_sheet(wb, nm, d, widths = de_widths(d))
    colour_direction(wb, nm, d)
    de_index[[length(de_index) + 1]] <- data.frame(
      sheet = nm, cluster = m$cluster[i], stratum = m$stratum[i], primary = m$is_primary[i],
      rows_in_sheet = nrow(d), genes_undetected_dropped = ndrop,
      csv = paste0("csv/", k, ".csv"), stringsAsFactors = FALSE)
  }
  for (i in which(!in_book))
    de_index[[length(de_index) + 1]] <- data.frame(
      sheet = "(CSV only)", cluster = m$cluster[i], stratum = m$stratum[i], primary = FALSE,
      rows_in_sheet = NA_integer_, genes_undetected_dropped = NA_integer_,
      csv = paste0("csv/", m$key[i], ".csv"), stringsAsFactors = FALSE)

  # ---- GO sheets
  if (!is.null(go) && nrow(go)) {
    bp <- go[go$ontology == "BP", ]
    for (cl in unique(bp$cluster)) {
      x <- pick(bp[bp$cluster == cl, ], go_cols)
      x <- x[order(x$direction, x$p.adjust), ]
      add_sheet(wb, sprintf("GO_BP_%s", cl), sig3(x),
                widths = c(rep("auto", 5), 55, rep("auto", 7), 60, rep("auto", 2), 50))
    }
    mfcc <- pick(go[go$ontology != "BP", ], go_cols)
    if (nrow(mfcc)) add_sheet(wb, "GO_MF_CC", sig3(mfcc[order(mfcc$cluster, mfcc$ontology, mfcc$p.adjust), ]),
                              widths = c(rep("auto", 5), 55, rep("auto", 7), 60, rep("auto", 2), 50))
  } else add_notes(wb, "GO_BP_none", "No GO term reached significance in any cluster or direction. See the GO_audit sheet for the input and universe sizes behind each test.")

  if (!is.null(gsea) && nrow(gsea)) {
    g <- pick(gsea, gsea_cols)
    add_sheet(wb, "GSEA", sig3(g[order(g$cluster, g$padj), ]),
              widths = c(rep("auto", 2), 45, rep("auto", 4), 80))
  }
  if (!is.null(ms)) {
    cm <- ms$comparison[ms$comparison$contrast == contrast_key, , drop = FALSE]
    add_notes(wb, "Mito_check_READ_ME", c(
      "MITOCHONDRIALLY-ENCODED GENES: A CAVEAT THAT CHANGES HOW THE KO-UP RESULT READS",
      "",
      "Every one of the seven subclusters has 3-7 mitochondrially-encoded genes (mt-Nd1, mt-Nd2, mt-Nd4,",
      "mt-Cytb, mt-Co1, mt-Rnr1/2, mt-Atp8) in its KO-UP list, and NOT ONE has an mt- gene on its KO-down",
      "side. Seven independently clustered cell populations do not agree that perfectly by biology. A",
      "one-directional shift of the whole mt- block, present in every cluster, is the classic signature of",
      "a difference in mitochondrial read fraction between the two libraries - a QC covariate, not a pathway.",
      "",
      "This matters because those genes ARE the 'oxidative phosphorylation', 'electron transport chain' and",
      "'ATP synthesis' GO terms. With them in, the headline KO-up story in several clusters is largely them.",
      "",
      "Measured mitochondrial share of log-normalised signal per arm is on the Mito_fraction sheet.",
      "GO BP re-run with the mt- genes removed is on the GO_BP_no_mt sheet. The table below counts what",
      "each direction loses. A term that survives the removal is about nuclear genes and can be read",
      "normally; a term that disappears was being carried by the mt- block.",
      "",
      "DOES THIS ALSO CONTAMINATE THE KO-DOWN SIDE? NO - and here is the arithmetic.",
      "If the mitochondrial share rises, every other gene's share must fall by the same total amount, so it",
      "is fair to ask whether the KO-down lists are just that squeeze. They are not, by two orders of",
      "magnitude. Taking CM1: the mt- share goes 1.35% -> 1.71% of signal, a ratio of 1.27, i.e. +0.34 in",
      "log2 - which is the size of the mt- fold changes actually observed (mt-Nd1 +0.54, mt-Nd2 +0.45,",
      "mt-Cytb +0.44, and similar in every cluster). The reciprocal effect on everything else is",
      "98.65% -> 98.29%, a ratio of 0.9964, or -0.005 in log2. That is fifty times smaller than the 0.25",
      "significance threshold. The KO-down results stand on their own, and the comparison table bears that",
      "out: removing the mt- genes costs the KO-down lists 0-2 GO terms out of 37-425, while it costs the",
      "KO-up lists nearly everything.",
      "",
      "This is a caveat about interpretation, not an error in the numbers: the mt- rows in the DE sheets",
      "are correct, they are flagged in the 'mito_encoded' column, and nothing has been deleted."))
    if (nrow(cm)) add_sheet(wb, "Mito_GO_comparison", cm, widths = c(rep("auto", 8), 60))
    add_sheet(wb, "Mito_fraction", ms$mito_wide)
    if (!is.null(ms$go_no_mt)) {
      g <- ms$go_no_mt[ms$go_no_mt$contrast == contrast_key, , drop = FALSE]
      if (nrow(g)) add_sheet(wb, "GO_BP_no_mt", sig3(pick(g, c(go_cols, "n_input_no_mt", "n_mt_removed"))),
                             widths = c(rep("auto", 5), 55, rep("auto", 7), 60, rep("auto", 4), 50))
    }
  }
  add_sheet(wb, "Secondary_tables", do.call(rbind, de_index))
  add_sheet(wb, "GO_audit", aud, widths = c(rep("auto", 9), 60))
  add_sheet(wb, "Skipped", if (nrow(de$skipped)) de$skipped else
            data.frame(note = "No cluster or arm fell below the 10-cell floor."))
  wb
}

# ---------------------------------------------------------------------------
cat("Building part 1 workbook ...\n")
wb1 <- build_workbook("P7_KO_vs_WT",
  title = "P7 KO vs P7 WT, within cardiomyocyte subclusters CM1, CM2, CM3, CM4, CM5, CM7, CM8",
  extra_notes = c(
    "  Both arms are P7, so the FACS cycling-enrichment applied at P7 is matched between them and the",
    "  all-cells stratum is the primary read. The G1-only sheets are a phase-matched sanity check;",
    "  CM4 has no G1 sheet because that subcluster contains no G1 cells at all (it is entirely S/G2M).",
    "  CM0 and CM6 are absent by design: CM6 contains ZERO P7 cells and CM0 only ~82 WT-P7 / 167 KO-P7,",
    "  so neither supports a P7 KO-vs-WT contrast."))

cat("Building part 2 workbook ...\n")
wb2 <- build_workbook("WT_P7_vs_P0",
  title = "P7 WT vs P0 WT, over all cardiomyocytes and within CM1, CM2, CM3, CM4, CM5, CM7, CM8",
  extra_notes = c(
    "  ON THE CELL-SORT CONFOUND, AND HOW BIG IT ACTUALLY IS.",
    "  P0 and P7 were not sorted identically, so in principle a P0-vs-P7 comparison can read out the sort",
    "  rather than development. Both a G1-matched and a raw version of every table therefore ship here.",
    "  We measured the size of the problem instead of assuming it, and it is small IN THIS COMPARISON:",
    "    - within cardiomyocytes the cycling fraction is 16.3% at P0-WT and 25.0% at P7-WT, a ratio of 1.53x",
    "      (the 4.5-5.2x enrichment figure quoted in the project notes is not the P0-vs-P7 ratio inside the",
    "      cardiomyocyte compartment of this dataset);",
    "    - G1-matched and raw log2 fold changes correlate at r = 0.99 with a median absolute difference of",
    "      0.001, and the two agree on ~450 of ~500 significant genes in all cardiomyocytes;",
    "    - cell-cycle genes (Mki67, Top2a, Cdk1, Ccnb1, Birc5, Aurkb, Ect2) move by less than 0.07 log2 units",
    "      between the two, and none of them is significant in either.",
    "  The G1-matched sheets are still the primary read - they are free insurance - but you are not choosing",
    "  between two different biological answers, and a gene that looks interesting in one will look",
    "  interesting in the other. CM4 has no G1 sheet: that subcluster contains no G1 cells at all.",
    "",
    "  For the maturation / metabolism / cell-cycle question, the 'sig_sets' column tags each gene with the",
    "  curated programs it belongs to; filter on it to pull those genes out directly. Note that Cpt1b is not",
    "  in the matrix at all, so it is absent rather than unchanged."))

saveWorkbook(wb1, file.path(OUT, "P7_KO_vs_WT_by_CM_subcluster.xlsx"), overwrite = TRUE)
saveWorkbook(wb2, file.path(OUT, "P7WT_vs_P0WT.xlsx"), overwrite = TRUE)

for (f in c("P7_KO_vs_WT_by_CM_subcluster.xlsx", "P7WT_vs_P0WT.xlsx")) {
  p <- file.path(OUT, f)
  cat(sprintf("%-42s %6.1f MB  sheets: %s\n", f, file.size(p) / 1e6,
              paste(getSheetNames(p), collapse = ", ")))
}
cat("DONE 03_excel.R\n")
