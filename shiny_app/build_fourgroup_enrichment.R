# build_fourgroup_enrichment.R
# ---------------------------------------------------------------------------
# GO + GSEA for the four-group CM contrasts, written back into app_data.rds as
# app$enrich$fourgroup. Drives the "Enrichment" panel of the Four-group tab.
#
# WHY THIS EXISTS. The app already had enrichment, but only on the KO-vs-WT
# tables POOLED across P0 and P7 (app$enrich$sub, from build_subcluster_enrichment.R).
# A collaborator asked for GO on the P7-specific KO-vs-WT lists per subcluster,
# and on P7-WT vs P0-WT -- neither of which the pooled tables answer. The DE for
# those contrasts has existed in app$fourgroup$de since build_fourgroup.R; GO was
# simply never run on it.
#
# It covers ALL four contrasts x both strata x every subcluster (77 tables), not
# just the two that were asked about, because the marginal cost is build time and
# the alternative is doing this again for the next question.
#
#   Run LOCALLY in a real R session (NOT the bash sandbox -- it segfaults on the
#   ~100 MB gzip rds), AFTER build_fourgroup.R:
#     Rscript shiny_app/build_fourgroup_enrichment.R --probe   # sizes only, writes nothing
#     Rscript shiny_app/build_fourgroup_enrichment.R           # compute + save
#
#   Options:
#     --probe        report input/universe sizes and estimated bundle growth, exit
#     --ont=BP       restrict ontologies (default BP,MF,CC)
#     --contrast=K   restrict to one contrast key (e.g. P7_KO_vs_WT)
#     --grid=de2     use the curated-panel DE grid instead of the broad one
#
# THE UNIVERSE, AND WHY IT IS NOT THE DE TABLE.
# app$fourgroup$de is ROW-GATED by build_fourgroup.R: a gene is kept only if it
# is expressed in >= 5% of one arm AND (padj < 0.05 OR |log2FC| >= 0.5). Handing
# that table to enrichGO as the universe would be badly wrong -- the gate has
# already removed most of the expressed-but-unchanging genes, which is exactly
# the background a hypergeometric test needs. So the universe is recomputed here
# from the expression matrix directly: genes detected in >= 5% of the cells of at
# least one arm of that specific contrast, in that specific cluster and stratum.
# The significant lists are unaffected -- the gate keeps everything with
# padj < 0.05, so the input genes are all present.
#
# DESCRIPTIVE ONLY. n = 1 animal per genotype x timepoint. Everything downstream
# of a pseudoreplicated p-value is a ranking, not a test.
# ---------------------------------------------------------------------------

suppressMessages({
  library(Matrix)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(fgsea)
  library(msigdbr)
})

ARGS   <- commandArgs(trailingOnly = TRUE)
PROBE  <- "--probe" %in% ARGS
argval <- function(flag, default = "") {
  hit <- grep(paste0("^", flag, "="), ARGS, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
GO_ONTS   <- strsplit(argval("--ont", "BP,MF,CC"), ",")[[1]]
ONE_CT    <- argval("--contrast", "")
GRID      <- argval("--grid", "de")

# Same cutoffs as the offline pipeline (analysis/2026-08-21_email/02_enrich.R),
# deliberately tighter than build_subcluster_enrichment.R's permissive 0.2/0.2.
SIG_PADJ  <- 0.05
SIG_LFC   <- 0.25
GO_PCUT   <- 0.05
GO_QCUT   <- 0.2
MIN_PCT   <- 5      # universe: detected in >= this % of one arm
MIN_INPUT <- 10     # below this a list is not worth testing

if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
FG  <- app$fourgroup
stopifnot(!is.null(FG), !is.null(FG$de), !is.null(app$deg_expr), !is.null(app$deg_meta))
CONF  <- app$confound
CTAB  <- FG$built$contrasts
stopifnot(!is.null(CTAB))

DE   <- if (identical(GRID, "de2") && !is.null(FG$de2)) FG$de2 else FG$de
M    <- app$deg_expr
MD   <- app$deg_meta[match(colnames(M), app$deg_meta$cell), , drop = FALSE]
MD$fourgrp <- paste(MD$genotype, MD$timepoint, sep = "-")
if (is.null(MD$cm_subcluster)) {
  cml <- setNames(paste0("CM", app$cm$meta[[paste0("SCT_snn_res.", FG$built$res)]]), app$cm$meta$cell)
  MD$cm_subcluster <- unname(cml[MD$cell])
}
MT <- grep("^mt-", rownames(M), value = TRUE)

cat(sprintf("grid: %s   %d clusters, %d tables\n", GRID, length(DE), sum(vapply(DE, length, 0L))))
cat(sprintf("ontologies: %s\n", paste(GO_ONTS, collapse = ", ")))

# ---- gene sets for GSEA (same recipe as build_subcluster_enrichment.R) ------
PATHWAYS <- NULL
if (!PROBE) {
  cat("Fetching MSigDB gene sets (Hallmark + KEGG_LEGACY, mouse) ...\n")
  PATHWAYS <- tryCatch({
    H  <- msigdbr(species = "Mus musculus", collection = "H")
    KG <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
    lapply(split(c(H$gene_symbol, KG$gene_symbol), c(H$gs_name, KG$gs_name)), unique)
  }, error = function(e) { cat("  !! MSigDB unavailable (", conditionMessage(e), ") -- skipping GSEA.\n"); NULL })
  cat(sprintf("  %d gene sets\n", length(PATHWAYS)))
}

# ---- the universe, from the matrix rather than the gated table --------------
arm_cells <- function(cluster, group, stratum) {
  inCl <- if (cluster == "AllCM") !is.na(MD$cm_subcluster) else
          (!is.na(MD$cm_subcluster) & MD$cm_subcluster == cluster)
  ph <- if (identical(stratum, "G1")) as.character(MD$Phase) == "G1" else TRUE
  which(inCl & ph & MD$fourgrp == group)
}
detected_pct <- function(idx) {
  if (!length(idx)) return(setNames(numeric(nrow(M)), rownames(M)))
  100 * Matrix::rowSums(M[, idx, drop = FALSE] > 0) / length(idx)
}
universe_of <- function(cluster, contrast_key, stratum) {
  r <- CTAB[match(contrast_key, CTAB$key), ]
  a <- arm_cells(cluster, r$A, stratum); b <- arm_cells(cluster, r$B, stratum)
  if (!length(a) || !length(b)) return(character(0))
  keep <- pmax(detected_pct(a), detected_pct(b)) >= MIN_PCT
  rownames(M)[keep]
}

# ---- gene lists, with a fallback ladder so "small" != "nothing" -------------
pick_genes <- function(d, sign) {
  ok <- !d$confounder & is.finite(d$padj) & is.finite(d$log2FoldChange)
  cand <- list(
    list(rule = sprintf("padj<%.2g & sign*log2FC>=%.2g", SIG_PADJ, SIG_LFC),
         g = d$gene[ok & d$padj < SIG_PADJ & sign * d$log2FoldChange >= SIG_LFC]),
    list(rule = sprintf("relaxed: padj<%.2g, any log2FC this direction", SIG_PADJ),
         g = d$gene[ok & d$padj < SIG_PADJ & sign * d$log2FoldChange > 0]),
    list(rule = "relaxed: top 200 by |log2FC| this direction (NOT significance-filtered)",
         g = { sub <- d[ok & sign * d$log2FoldChange > 0, , drop = FALSE]
               head(sub$gene[order(-abs(sub$log2FoldChange))], 200) }))
  for (cc in cand) if (length(unique(cc$g)) >= MIN_INPUT) return(list(genes = unique(cc$g), rule = cc$rule))
  list(genes = unique(cand[[1]]$g),
       rule = paste0(cand[[1]]$rule, sprintf(" (below the %d-gene floor even relaxed; not tested)", MIN_INPUT)))
}

run_go <- function(genes, universe, ont) {
  if (length(genes) < MIN_INPUT || !length(universe)) return(NULL)
  eg <- tryCatch(suppressWarnings(suppressMessages(
    enrichGO(gene = unique(genes), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = ont,
             universe = unique(universe), pAdjustMethod = "BH",
             pvalueCutoff = GO_PCUT, qvalueCutoff = GO_QCUT,
             minGSSize = 10, maxGSSize = 500, readable = FALSE))),
    error = function(e) NULL)
  if (is.null(eg)) return(NULL)
  df <- as.data.frame(eg); if (!nrow(df)) NULL else df
}

run_gsea <- function(d) {
  if (is.null(PATHWAYS)) return(NULL)
  dd <- d[!d$confounder & is.finite(d$pvalue) & is.finite(d$log2FoldChange), ]
  stat <- sign(dd$log2FoldChange) * -log10(pmax(dd$pvalue, 1e-300))
  names(stat) <- dd$gene
  stat <- stat[!duplicated(names(stat))]
  stat <- sort(stat[is.finite(stat)], decreasing = TRUE)
  if (length(stat) < 15) return(NULL)
  set.seed(1)
  fg <- tryCatch(suppressWarnings(fgsea(PATHWAYS, stat, minSize = 10, maxSize = 500)),
                 error = function(e) NULL)
  if (is.null(fg) || !nrow(fg)) return(NULL)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, "", collapse = ", ")
  as.data.frame(fg)
}

# ---- walk every table -------------------------------------------------------
keys <- CTAB$key
if (nzchar(ONE_CT)) keys <- intersect(keys, ONE_CT)
clusters <- names(DE)

go_rows <- list(); gsea_rows <- list(); audit <- list(); n_done <- 0L
for (cl in clusters) for (ck in keys) for (st in c("all", "G1")) {
  d <- DE[[cl]][[paste0(ck, "__", st)]]
  if (is.null(d) || !nrow(d)) next
  n_done <- n_done + 1L
  r    <- CTAB[match(ck, CTAB$key), ]
  uni  <- universe_of(cl, ck, st)
  # Direction labels come from the contrast table, not hard-coded: "up in KO" is
  # simply wrong for WT_P0_vs_P7, which has no KO in it.
  dirs <- list(list(key = "A_up", lab = r$pos, sign =  1),
               list(key = "B_up", lab = r$neg, sign = -1))
  cat(sprintf("[%3d] %-6s %-14s %-3s  universe=%5d\n", n_done, cl, ck, st, length(uni)))

  for (dd in dirs) {
    pk <- pick_genes(d, dd$sign)
    n_mt <- sum(pk$genes %in% MT)
    if (PROBE) {
      audit[[length(audit) + 1]] <- data.frame(
        cluster = cl, contrast = ck, stratum = st, direction = dd$key, direction_label = dd$lab,
        ontology = NA_character_, n_input = length(pk$genes), n_mt_input = n_mt,
        n_universe = length(uni), n_terms = NA_integer_, input_rule = pk$rule,
        stringsAsFactors = FALSE)
      next
    }
    for (ont in GO_ONTS) {
      g  <- run_go(pk$genes, uni, ont)
      nr <- if (is.null(g)) 0L else nrow(g)
      if (!is.null(g)) {
        g$cluster <- cl; g$contrast <- ck; g$stratum <- st; g$grid <- GRID
        g$direction <- dd$key; g$direction_label <- dd$lab; g$ontology <- ont
        g$n_input <- length(pk$genes); g$n_universe <- length(uni); g$input_rule <- pk$rule
        g$n_mt_input <- n_mt
        go_rows[[paste(cl, ck, st, dd$key, ont)]] <- g
      }
      audit[[length(audit) + 1]] <- data.frame(
        cluster = cl, contrast = ck, stratum = st, direction = dd$key, direction_label = dd$lab,
        ontology = ont, n_input = length(pk$genes), n_mt_input = n_mt,
        n_universe = length(uni), n_terms = nr, input_rule = pk$rule, stringsAsFactors = FALSE)
      if (ont == "BP") cat(sprintf("        %-10s n=%-4d GO BP: %3d terms%s\n", dd$lab,
                                   length(pk$genes), nr,
                                   if (n_mt) sprintf("  (%d mt- genes in the input)", n_mt) else ""))
    }
  }
  if (!PROBE) {
    fg <- run_gsea(d)
    if (!is.null(fg)) {
      fg$cluster <- cl; fg$contrast <- ck; fg$stratum <- st; fg$grid <- GRID
      fg$pos_label <- r$pos; fg$neg_label <- r$neg
      gsea_rows[[paste(cl, ck, st)]] <- fg
    }
  }
}

audit <- do.call(rbind, audit)

if (PROBE) {
  cat("\n== probe: nothing written ==\n")
  cat(sprintf("tables: %d   direction-tests: %d\n", n_done, nrow(audit)))
  cat(sprintf("input list sizes  : median %d, range %d-%d\n",
              median(audit$n_input), min(audit$n_input), max(audit$n_input)))
  cat(sprintf("universe sizes    : median %d, range %d-%d\n",
              median(audit$n_universe), min(audit$n_universe), max(audit$n_universe)))
  cat(sprintf("below the %d-gene floor: %d of %d direction-tests\n",
              MIN_INPUT, sum(audit$n_input < MIN_INPUT), nrow(audit)))
  cat(sprintf("mt- genes in input: %d tests affected, up to %d genes\n",
              sum(audit$n_mt_input > 0), max(audit$n_mt_input)))
  cat(sprintf("\nestimated GO calls: %d  (x %d ontologies)\n",
              nrow(audit) * length(GO_ONTS), length(GO_ONTS)))
  cat("Re-run without --probe to compute and save.\n")
  quit(save = "no")
}

bind <- function(L) { L <- L[!vapply(L, is.null, logical(1))]
                      if (length(L)) do.call(rbind, L) else NULL }
GO   <- bind(go_rows)
GSEA <- bind(gsea_rows)

# Add a sub-slot; do NOT replace app$enrich, which already holds $go/$gsea/$tf
# (cell-type level) and $sub (per-subcluster, pooled over timepoints).
before <- names(app$enrich)
app$enrich$fourgroup <- list(
  go = GO, gsea = GSEA, audit = audit,
  params = list(sig_padj = SIG_PADJ, sig_lfc = SIG_LFC, go_pcut = GO_PCUT, go_qcut = GO_QCUT,
                min_pct_universe = MIN_PCT, min_input = MIN_INPUT, onts = GO_ONTS,
                n_pathways = length(PATHWAYS)),
  built = list(when = as.character(Sys.time()), grid = GRID, res = FG$built$res,
               n_tables = n_done, mt_genes = MT,
               universe = sprintf("detected in >= %d%% of one arm, from app$deg_expr (NOT the gated DE table)", MIN_PCT))
)
lost <- setdiff(before, names(app$enrich))
if (length(lost)) stop("!! would drop existing enrich slots: ", paste(lost, collapse = ", "))

cat(sprintf("\n== summary ==\n  GO rows   : %s\n  GSEA rows : %s\n  audit rows: %d\n",
            if (is.null(GO)) 0 else nrow(GO), if (is.null(GSEA)) 0 else nrow(GSEA), nrow(audit)))
cat(sprintf("  slot size : %.1f MB\n", as.numeric(object.size(app$enrich$fourgroup)) / 1e6))

cat("\nBacking up -> app_data.pre_fgenrich.bak.rds\n")
file.copy("app_data.rds", "app_data.pre_fgenrich.bak.rds", overwrite = TRUE)
cat("Saving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
