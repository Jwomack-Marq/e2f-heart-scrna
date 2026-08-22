# 02_enrich.R -- GO + GSEA + every figure, for the 2026-08-21 deliverable.
# ---------------------------------------------------------------------------
# Reads de_tables.rds from 01_de.R. For each DE table, the KO-up and KO-down
# lists are enriched SEPARATELY (that is what the collaborator asked for), over
# GO BP / MF / CC, plus MSigDB Hallmark + KEGG_LEGACY GSEA on the full ranking.
#
# The enrichGO call, the GSEA gene sets and the ranking statistic are the same
# ones build_subcluster_enrichment.R uses, so these results are methodologically
# continuous with the website's "Subcluster enrichment" tab. Two deliberate
# differences: the input lists use the agreed padj<0.05 & |log2FC|>=0.25 rather
# than that script's |log2FC|>=1, and the GO p/q cutoffs are 0.05/0.2 rather than
# its deliberately permissive 0.2/0.2.
#
# Plot styling is ported from app.R so the figures match what the collaborator
# has already been looking at: go_dotplot_gg() (app.R:279) and gsea_barplot_gg()
# (app.R:265), saved through the same ggsave settings as dl_ggplot() (app.R:541).
# ---------------------------------------------------------------------------

suppressMessages({
  library(Matrix); library(ggplot2); library(clusterProfiler)
  library(org.Mm.eg.db); library(fgsea); library(msigdbr)
})

OUT <- "/out"
PLOT_ROOT <- file.path(OUT, "plots")
dir.create(file.path(PLOT_ROOT, "part1"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PLOT_ROOT, "part2"), recursive = TRUE, showWarnings = FALSE)

GO_PCUT   <- 0.05
GO_QCUT   <- 0.2
GO_ONTS   <- c("BP", "MF", "CC")
MIN_PCT   <- 5      # universe gate: expressed in >=5% of at least one arm
MIN_INPUT <- 10     # enrichGO needs a list worth testing
TOPN      <- 20

de <- readRDS(file.path(OUT, "de_tables.rds"))
P  <- de$params
cat(sprintf("loaded %d DE tables (sig: padj<%.2g & |log2FC|>=%.2g)\n",
            length(de$tables), P$sig_padj, P$sig_lfc))

part_of <- function(contrast) if (contrast == "P7_KO_vs_WT") "part1" else "part2"

# ---- gene sets for GSEA (identical recipe to build_subcluster_enrichment.R) --
# The image pre-warms the msigdbr cache at R_USER_CACHE_DIR so this needs no network.
# If it fails anyway, GSEA is skipped with a loud note rather than taking the whole
# run down -- GO is the part that was actually asked for.
cat("Fetching MSigDB gene sets (Hallmark + KEGG_LEGACY, mouse) ...\n")
PATHWAYS <- tryCatch({
  H  <- msigdbr(species = "Mus musculus", collection = "H")
  KG <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  lapply(split(c(H$gene_symbol, KG$gene_symbol), c(H$gs_name, KG$gs_name)), unique)
}, error = function(e) { cat("  !! MSigDB unavailable (", conditionMessage(e), ") -- GSEA will be skipped.\n"); NULL })
cat(sprintf("  %d gene sets\n", length(PATHWAYS)))

# ---- helpers ---------------------------------------------------------------
pct_cols <- function(d) grep("^pct_", names(d), value = TRUE)

# The universe is the expressed-gene space of THAT table -- not all 24,221 genes
# (which would inflate every fold enrichment) and certainly not the significant
# list. A gene that was tested and did not move belongs in the universe.
universe_of <- function(d) {
  pc <- pct_cols(d)
  unique(d$gene[do.call(pmax, c(d[pc], list(na.rm = TRUE))) >= MIN_PCT])
}

# Fallback ladder, so a thin cluster yields "the list was too small" rather than
# a silent empty table that reads as "nothing is enriched here".
pick_genes <- function(d, dirlabel, sign) {
  ok <- !d$confounder & is.finite(d$padj) & is.finite(d$log2FoldChange)
  cand <- list(
    list(rule = sprintf("padj<%.2g & %slog2FC%s%.2g", P$sig_padj,
                        if (sign > 0) "" else "", if (sign > 0) ">=+" else "<=-", P$sig_lfc),
         g = d$gene[ok & d$padj < P$sig_padj & sign * d$log2FoldChange >= P$sig_lfc]),
    list(rule = sprintf("relaxed: padj<%.2g, any log2FC in this direction", P$sig_padj),
         g = d$gene[ok & d$padj < P$sig_padj & sign * d$log2FoldChange > 0]),
    list(rule = "relaxed: top 200 by |log2FC| in this direction (NOT significance-filtered)",
         g = head(d$gene[ok & sign * d$log2FoldChange > 0][
                    order(-abs(d$log2FoldChange[ok & sign * d$log2FoldChange > 0]))], 200)))
  for (cc in cand) if (length(unique(cc$g)) >= MIN_INPUT) return(list(genes = unique(cc$g), rule = cc$rule))
  list(genes = unique(cand[[1]]$g), rule = paste0(cand[[1]]$rule, " (below the ", MIN_INPUT,
       "-gene floor even after relaxing; not tested)"))
}

run_go <- function(genes, universe, ont) {
  if (length(genes) < MIN_INPUT) return(NULL)
  eg <- tryCatch(suppressWarnings(suppressMessages(
    enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = ont,
             universe = unique(universe), pAdjustMethod = "BH",
             pvalueCutoff = GO_PCUT, qvalueCutoff = GO_QCUT,
             minGSSize = 10, maxGSSize = 500, readable = FALSE))),
    error = function(e) NULL)
  if (is.null(eg)) return(NULL)
  df <- as.data.frame(eg)
  if (!nrow(df)) NULL else df
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

# ---- plotting (ported from app.R) ------------------------------------------
save_plot <- function(p, part, name, w = 9, h = 6) {
  for (fmt in c("png", "pdf"))
    ggsave(file.path(PLOT_ROOT, part, paste0(name, ".", fmt)), p,
           device = fmt, width = w, height = h, dpi = 300, units = "in")
}
placeholder <- function(ttl, msg) {
  ggplot() + annotate("text", 0, 0, label = paste(strwrap(msg, 60), collapse = "\n"),
                      size = 4.2, colour = "grey25") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) + labs(title = ttl)
}
# app.R:279 go_dotplot_gg
go_dotplot <- function(d, ttl, topn = TOPN) {
  d <- head(d[order(d$p.adjust), ], topn)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  ggplot(d, aes(FoldEnrichment, Description, size = Count, color = p.adjust)) +
    geom_point() + scale_color_viridis_c(option = "magma", direction = -1) +
    theme_minimal(base_size = 11) +
    labs(x = "fold enrichment", y = NULL, color = "padj", size = "genes", title = ttl)
}
go_dotplot_facet <- function(d, ttl, topn = TOPN) {
  d <- do.call(rbind, lapply(split(d, d$direction), function(x) head(x[order(x$p.adjust), ], topn)))
  d$Description <- factor(d$Description, levels = rev(unique(d$Description[order(d$direction, d$p.adjust)])))
  ggplot(d, aes(FoldEnrichment, Description, size = Count, color = p.adjust)) +
    geom_point() + scale_color_viridis_c(option = "magma", direction = -1) +
    facet_wrap(~ direction, scales = "free_y") +
    theme_minimal(base_size = 10) +
    labs(x = "fold enrichment", y = NULL, color = "padj", size = "genes", title = ttl)
}
# app.R:265 gsea_barplot_gg
gsea_barplot <- function(d, ttl, up_lab, down_lab, topn = TOPN) {
  d <- head(d[order(-abs(d$NES)), ], topn)
  d$dir <- ifelse(d$NES > 0, up_lab, down_lab)
  d$pathway <- factor(d$pathway, levels = d$pathway[order(d$NES)])
  ggplot(d, aes(NES, pathway, fill = dir)) +
    geom_col() + geom_vline(xintercept = 0, color = "grey60") +
    scale_fill_manual(values = setNames(c("#c62828", "#1565c0"), c(up_lab, down_lab))) +
    theme_minimal(base_size = 10) +
    labs(x = sprintf("NES (>0 enriched toward %s)", up_lab), y = NULL, fill = NULL, title = ttl)
}
volcano <- function(d, ttl, up_lab, down_lab, nlab = 20) {
  d$y <- -log10(pmax(d$pvalue, 1e-300))
  d$grp <- ifelse(d$confounder, "sex/construct confounder", d$direction)
  cols <- setNames(c("#c62828", "#1565c0", "grey75", "#f9a825"),
                   c(up_lab, down_lab, "ns", "sex/construct confounder"))
  lab <- d[d$direction != "ns" & !d$confounder, ]
  lab <- head(lab[order(-abs(lab$log2FoldChange) * -log10(pmax(lab$padj, 1e-300))), ], nlab)
  p <- ggplot(d, aes(log2FoldChange, y, colour = grp)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = cols, name = NULL) +
    geom_vline(xintercept = c(-P$sig_lfc, P$sig_lfc), linetype = "dotted", colour = "grey45") +
    theme_minimal(base_size = 11) +
    labs(x = sprintf("log2 fold change  (>0 = %s)", up_lab), y = "-log10(p)", title = ttl)
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab))
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                                      max.overlaps = 30, show.legend = FALSE)
  p
}
# app.R:362 lfc_heat
lfc_heat <- function(tabs, ttl, topn = 25) {
  allg <- unique(unlist(lapply(tabs, `[[`, "gene")))
  M <- vapply(tabs, function(d) d$log2FoldChange[match(allg, d$gene)], numeric(length(allg)))
  rownames(M) <- allg
  score <- apply(abs(M), 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
  score[allg %in% de$confound] <- NA
  top <- names(sort(score, decreasing = TRUE))[seq_len(min(topn, sum(!is.na(score))))]
  long <- expand.grid(gene = top, grp = colnames(M), stringsAsFactors = FALSE)
  long$lfc <- M[cbind(match(long$gene, rownames(M)), match(long$grp, colnames(M)))]
  long$gene <- factor(long$gene, levels = rev(top))
  ggplot(long, aes(grp, gene, fill = lfc)) + geom_tile(color = "grey92") +
    scale_fill_gradient2(low = "#1565c0", mid = "white", high = "#c62828", midpoint = 0) +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "log2FC", title = ttl)
}

# ---- run -------------------------------------------------------------------
go_rows <- list(); gsea_rows <- list(); audit <- list()
man <- de$manifest

for (i in seq_len(nrow(man))) {
  key <- man$key[i]; d <- de$tables[[key]]
  part <- part_of(man$contrast[i])
  uni  <- universe_of(d)
  up_lab <- man$up_label[i]; down_lab <- man$down_label[i]
  ttl_base <- sprintf("%s — %s (%s cells)", man$cluster[i], man$contrast[i], man$stratum[i])
  cat(sprintf("\n[%2d/%2d] %s   universe=%d\n", i, nrow(man), key, length(uni)))

  go_this <- list(); n_input <- c(); rule_of <- c()
  for (dirn in c(up_lab, down_lab)) {
    sgn <- if (dirn == up_lab) 1 else -1
    pk  <- pick_genes(d, dirn, sgn)
    n_input[[dirn]] <- length(pk$genes); rule_of[[dirn]] <- pk$rule
    cat(sprintf("   %-8s n=%-5d rule: %s\n", dirn, length(pk$genes), pk$rule))
    for (ont in GO_ONTS) {
      g <- run_go(pk$genes, uni, ont)
      nr <- if (is.null(g)) 0L else nrow(g)
      if (!is.null(g)) {
        g$contrast <- man$contrast[i]; g$cluster <- man$cluster[i]
        g$stratum <- man$stratum[i]; g$direction <- dirn; g$ontology <- ont
        g$n_input <- length(pk$genes); g$n_universe <- length(uni); g$input_rule <- pk$rule
        go_rows[[paste(key, dirn, ont)]] <- g
        if (ont == "BP") go_this[[dirn]] <- g
      }
      audit[[length(audit) + 1]] <- data.frame(
        key = key, contrast = man$contrast[i], cluster = man$cluster[i], stratum = man$stratum[i],
        direction = dirn, ontology = ont, n_input = length(pk$genes), n_universe = length(uni),
        n_terms = nr, input_rule = pk$rule, stringsAsFactors = FALSE)
      if (ont == "BP") cat(sprintf("      GO %s: %d terms\n", ont, nr))
    }
  }

  # ---- figures (primary stratum only; the secondary tables still ship as data)
  if (man$is_primary[i]) {
    tag <- sprintf("%s_%s", man$contrast[i], man$cluster[i])
    for (dirn in names(go_this))
      save_plot(go_dotplot(go_this[[dirn]], sprintf("GO BP — %s, %s", ttl_base, dirn)),
                part, sprintf("go_bp_%s_%s", tag, dirn), w = 9, h = 6)
    for (dirn in setdiff(c(up_lab, down_lab), names(go_this)))
      save_plot(placeholder(sprintf("GO BP — %s, %s", ttl_base, dirn),
                  sprintf("No GO BP term reached padj < %.2g / q < %.2g. Input list: %d genes, selected by [%s]. Universe: %d expressed genes. An empty panel means the test ran and found nothing, or the list was below the %d-gene floor — the GO_audit sheet says which.",
                          GO_PCUT, GO_QCUT, n_input[[dirn]], rule_of[[dirn]], length(uni), MIN_INPUT)),
                part, sprintf("go_bp_%s_%s", tag, dirn), w = 9, h = 6)
    if (length(go_this) == 2) {
      both <- rbind(go_this[[1]], go_this[[2]])
      save_plot(go_dotplot_facet(both, sprintf("GO BP — %s", ttl_base)),
                part, sprintf("go_bp_%s_BOTH", tag), w = 13, h = 6.5)
    }
    save_plot(volcano(d, sprintf("%s", ttl_base), up_lab, down_lab),
              part, sprintf("volcano_%s", tag), w = 8, h = 6.5)
  }

  # ---- GSEA
  fg <- run_gsea(d)
  if (!is.null(fg)) {
    fg$contrast <- man$contrast[i]; fg$cluster <- man$cluster[i]; fg$stratum <- man$stratum[i]
    gsea_rows[[key]] <- fg
    cat(sprintf("      GSEA: %d sets tested, %d at padj<0.05\n", nrow(fg), sum(fg$padj < 0.05, na.rm = TRUE)))
    if (man$is_primary[i]) {
      sig <- fg[is.finite(fg$padj) & fg$padj < 0.25, ]
      pl <- if (nrow(sig)) gsea_barplot(sig, sprintf("GSEA (Hallmark + KEGG) — %s", ttl_base), up_lab, down_lab)
            else placeholder(sprintf("GSEA — %s", ttl_base), "No gene set reached padj < 0.25.")
      save_plot(pl, part, sprintf("gsea_%s_%s", man$contrast[i], man$cluster[i]), w = 10, h = 6.5)
    }
  } else cat("      GSEA: not run (fewer than 15 rankable genes)\n")
}

# ---- cross-cluster overviews ----------------------------------------------
for (ck in unique(man$contrast)) {
  sel <- man[man$contrast == ck & man$is_primary & man$cluster != "AllCM", ]
  if (nrow(sel) < 2) next
  tabs <- setNames(lapply(sel$key, function(k) de$tables[[k]]), sel$cluster)
  save_plot(lfc_heat(tabs, sprintf("Top responsive genes x subcluster — %s (%s cells)",
                                   ck, sel$stratum[1])),
            part_of(ck), sprintf("lfc_heatmap_%s", ck), w = 8, h = 8)

  cnt <- data.frame(cluster = factor(sel$cluster, levels = sel$cluster),
                    up = sel$n_up, down = -sel$n_down)
  long <- rbind(data.frame(cluster = cnt$cluster, n = cnt$up,   dir = sel$up_label[1]),
                data.frame(cluster = cnt$cluster, n = cnt$down, dir = sel$down_label[1]))
  p <- ggplot(long, aes(cluster, n, fill = dir)) +
    geom_col() + geom_hline(yintercept = 0, colour = "grey40") +
    scale_fill_manual(values = setNames(c("#c62828", "#1565c0"), c(sel$up_label[1], sel$down_label[1]))) +
    scale_y_continuous(labels = abs) + theme_minimal(base_size = 12) +
    labs(x = NULL, y = "significant genes", fill = NULL,
         title = sprintf("Significant genes per subcluster — %s", ck),
         subtitle = sprintf("padj < %.2g and |log2FC| >= %.2g; %s cells", P$sig_padj, P$sig_lfc, sel$stratum[1]))
  save_plot(p, part_of(ck), sprintf("counts_%s", ck), w = 8, h = 5.5)
}

bind <- function(L) { L <- L[!vapply(L, is.null, logical(1))]
                      if (length(L)) do.call(rbind, L) else NULL }
res <- list(go = bind(go_rows), gsea = bind(gsea_rows),
            audit = do.call(rbind, audit),
            params = c(P, list(go_pcut = GO_PCUT, go_qcut = GO_QCUT, min_pct_universe = MIN_PCT,
                               min_input = MIN_INPUT, n_pathways = length(PATHWAYS))))
saveRDS(res, file.path(OUT, "enrich.rds"))
write.csv(res$audit, file.path(OUT, "csv", "_go_audit.csv"), row.names = FALSE)
if (!is.null(res$go))   write.csv(res$go,   file.path(OUT, "csv", "_go_all.csv"),   row.names = FALSE)
if (!is.null(res$gsea)) write.csv(res$gsea, file.path(OUT, "csv", "_gsea_all.csv"), row.names = FALSE)

cat(sprintf("\nGO rows: %d | GSEA rows: %d | figures: %d\n",
            if (is.null(res$go)) 0 else nrow(res$go),
            if (is.null(res$gsea)) 0 else nrow(res$gsea),
            length(list.files(PLOT_ROOT, recursive = TRUE))))
cat("DONE 02_enrich.R\n")
