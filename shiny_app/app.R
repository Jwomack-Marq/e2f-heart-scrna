# E2F7/8 mouse-heart scRNA-seq — interactive cell browser (server-side Shiny app).
# Self-contained: ggplot2 + Matrix + plotly + DT on the slim app_data.rds built by
# build_app_data.R (no Seurat needed at runtime). Previously also shipped as a
# shinylive/WebAssembly static export — that build is archived at the git tag
# `shinylive-static-archive`; this app now runs only against an R/Shiny server.
#
#   Local preview : shiny::runApp("shiny_app")
#   Deploy        : rsconnect::deployApp("shiny_app")   # shinyapps.io (account jwomackmu)
#
# DESCRIPTIVE pilot (n = 1, sex-confounded, KO not transcript-confirmed) — see About tab.

library(shiny)
library(bslib)
library(ggplot2)
library(Matrix)
library(plotly)
library(DT)
library(svglite)         # vector SVG device for ggsave() figure export
library(shinycssloaders) # loading spinners on plot outputs
library(openxlsx)        # multi-sheet .xlsx export (shared with the offline analysis scripts)

# Workbook/CSV conventions shared with analysis/<date>/03_excel.R, so the file a
# collaborator gets by email and the one they get by clicking Download match.
source("download_helpers.R", local = TRUE)

# "Figure Studio" handoff: sends a built ggplot to the companion editor app for
# publication styling. Inert (no UI, no observers) unless FIGURE_STUDIO_BASE and
# HANDOFF_DIR are set -- see studio_helpers.R for the contract.
source("studio_helpers.R", local = TRUE)

app   <- readRDS("app_data.rds")
meta  <- app$meta; expr <- app$expr; genes <- app$genes
cmm   <- app$cm$meta; RES <- setdiff(app$cm$res, c("0.1", "0.3"))  # only expose res 0.2 in the app
heat  <- app$heat; tabs <- app$tables; CONF <- app$confound
ctDE  <- tabs$ct_DE; subDE <- tabs$sub_DE; subSum <- tabs$sub_summary; subType <- tabs$sub_subtype
figs  <- app$figs
GI    <- app$geneinfo          # per-gene info table (NULL on an un-enriched build)
ENR   <- app$enrich            # list(gsea, go, tf) of precomputed enrichment (may be NULL)
EXPR  <- app$deg_expr          # broad-gene log-norm matrix (genes x downsampled cells) for the DEG explorer
DMETA <- if (!is.null(app$deg_meta)) app$deg_meta else meta   # metadata aligned to EXPR cols (fallback so UI builds)
GENES_FULL <- app$deg_genes
SCOREMETA  <- app$score_meta   # per-cell module-score definitions/coverage (build_signature_scores.R; may be NULL)
COMMUN     <- app$commun       # curated ligand-receptor scores (build_communication.R; may be NULL)
REFMAP     <- app$refmap       # reference-marker annotation check (build_refmap.R; may be NULL)
FG         <- app$fourgroup    # four-group CM analysis (build_fourgroup.R; may be NULL)
# Present WT before KO on every plot (genotype is otherwise alphabetical -> KO first).
GENO_LEVELS <- c("WT","KO")
relevel_geno <- function(df) {
  if (!is.null(df) && "genotype" %in% names(df))
    df$genotype <- factor(df$genotype, levels = intersect(GENO_LEVELS, unique(as.character(df$genotype))))
  df
}
meta <- relevel_geno(meta); cmm <- relevel_geno(cmm); DMETA <- relevel_geno(DMETA)
# searchable union: curated panel (full 30k cells) + broad genes (deg_expr, ~8k cells).
# Lets the UMAP/Gene-detail views show ANY gene that shows up in a volcano, falling
# back to the broad matrix when a gene is outside the curated panel (expr_vec below).
ALL_GENES <- sort(unique(c(genes, GENES_FULL)))
in_panel  <- function(g) !is.null(g) && nzchar(g) && g %in% rownames(expr)

# Curated gene sets for the "Gene set" quick-pick (mouse symbols, copied from the
# pipeline's _common.R). Each is intersected with the available genes; empties drop.
# Note: every gene here lives in the curated panel, so picking via a set => full cells.
.gene_sets_raw <- list(
  "E2F targets" = c("Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7","Pcna","Cdc6","Cdt1",
                    "Ccne1","Ccne2","Ccna2","Ccnb1","Ccnb2","Cdk1","Cdc20","E2f1",
                    "Mki67","Top2a","Birc5","Bub1","Rrm2","Tk1","Cdc25a","Foxm1",
                    "Aurkb","Plk1","Cenpa","Cenpe"),
  "E2F repressors" = c("E2f7","E2f8"),
  "Cell cycle (S)" = c("Mcm2","Mcm4","Mcm5","Mcm6","Mcm7","Pcna","Rrm1","Rrm2","Tyms",
                       "Slbp","Gins2","Cdc6","Cdc45","Cdca7","Dtl","Uhrf1","Usp1","Hells"),
  "Cell cycle (G2/M)" = c("Mki67","Top2a","Ccnb1","Ccnb2","Cdk1","Cdc20","Aurka","Aurkb",
                          "Bub1","Birc5","Cenpa","Cenpe","Cenpf","Cks2","Kif11","Kif23",
                          "Tpx2","Nusap1","Anln","Cdca2","Cdca3","Cdca8","Ube2c","Tubb4b"),
  "Cardiomyocyte" = c("Tnnt2","Myh6","Actc1","Ttn","Tnni3","Nppa"),
  "Fibroblast"    = c("Col1a1","Col1a2","Dcn","Pdgfra","Gsn"),
  "Endothelial"   = c("Pecam1","Cdh5","Kdr","Fabp4","Egfl7"),
  "Immune / myeloid" = c("Ptprc","Cd68","Lyz2","C1qa","Csf1r"),
  "Mural / pericyte" = c("Rgs5","Pdgfrb","Myh11","Acta2","Tagln"),
  "CM subtypes"   = c("Myl2","Myh7","Myl7","Sln","Nppa","Bmp10","Hey2","Irx3","Tbx20",
                      "Mki67","Top2a","Ccnb1","Aurkb","Cdca8"),
  "CM maturation" = c("Myh6","Tnni3","Pln","Atp2a2","Ckm","Myl2","Cox6a2",          # mature
                      "Myh7","Tnni1","Nppa","Nppb","Ccnd1","Mki67","Top2a"))         # immature
GENE_SETS <- Filter(length, lapply(.gene_sets_raw, function(g) intersect(g, ALL_GENES)))
# choices for the "Gene set" dropdowns: sentinel "__all__" maps to the full list
GENE_SET_CHOICES <- c("All genes" = "__all__", setNames(names(GENE_SETS), names(GENE_SETS)))
genes_for_set <- function(set) if (is.null(set) || set == "__all__" || !set %in% names(GENE_SETS))
  ALL_GENES else GENE_SETS[[set]]

has <- function(col, df = meta) col %in% names(df)
CAT_COLS  <- Filter(has, c("celltype","genotype","timepoint","Phase","cycling","cm_subtype",
                           "cm_subcluster","seurat_clusters"))
# per-cell module scores added by build_signature_scores.R (headline scores only;
# absent columns drop out so an un-rebuilt bundle still loads).
SCORE_COLS <- Filter(has, c("sig_prolif","sig_cytokinesis","sig_ccexit","sig_ploidy",
                            "sig_maturation","sig_maturation_nocc","sig_metabolic"))
CONT_COLS <- Filter(has, c("pseudotime","S.Score","G2M.Score", SCORE_COLS))

nice <- c(gene = "Gene expression", celltype = "Cell type", genotype = "Genotype (KO/WT)",
          timepoint = "Timepoint (P0/P7)", Phase = "Cell-cycle phase", cycling = "Cycling (S/G2M)",
          cm_subtype = "CM subtype", seurat_clusters = "Cluster", pseudotime = "Pseudotime",
          S.Score = "S-phase score", G2M.Score = "G2/M score", subcluster = "Subcluster",
          splitgrp = "Genotype × timepoint", cm_subcluster = "CM subcluster (res 0.2)",
          sig_maturation_nocc = "CM maturation (cycle-free)",
          sig_mat_immature_nocc = "Immature-CM program (cycle-free)",
          sig_prolif = "Proliferation score", sig_cytokinesis = "Cytokinesis score",
          sig_ccexit = "Cell-cycle exit score", sig_ploidy = "Polyploidization proxy",
          sig_maturation = "CM maturation score", sig_metabolic = "Metabolic maturation (FAO−glyc)",
          sig_glycolysis = "Glycolysis score", sig_faox = "Fatty-acid oxidation score",
          sig_mat_mature = "Mature-CM program", sig_mat_immature = "Immature-CM program")
labof <- function(x) ifelse(x %in% names(nice), nice[x], gsub("_", " ", x))

theme_umap <- theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "right", strip.text = element_text(face = "bold"))
div_scale <- scale_fill_gradient2(low = "#1565c0", mid = "white", high = "#c62828", midpoint = 0, na.value = "grey92")

# expression for `cells`, drawn from the curated panel when possible, else the
# broad deg_expr matrix (only ~8k cells overlap, the rest stay NA = grey). This is
# what lets a gene that appears in a volcano but not the curated panel be coloured.
expr_vec <- function(gene, cells) {
  if (in_panel(gene)) { v <- as.numeric(expr[gene, ]); names(v) <- colnames(expr); return(v[cells]) }
  if (!is.null(gene) && nzchar(gene) && !is.null(EXPR) && gene %in% rownames(EXPR)) {
    v <- as.numeric(EXPR[gene, ]); names(v) <- colnames(EXPR); return(v[cells])
  }
  setNames(rep(NA_real_, length(cells)), cells)
}
cm_subcol <- function(res) paste0("SCT_snn_res.", res)
cm_subs   <- function(res) { v <- unique(paste0("CM", cmm[[cm_subcol(res)]]))
                             v[order(as.integer(sub("CM", "", v)))] }
# label a subcluster with its nearest CM subtype + size, for dropdowns/titles
sub_label <- function(res, sub) {
  st <- subType[[paste0("res", res)]]; sm <- subSum[[paste0("res", res)]]
  lab <- sub
  if (!is.null(st) && sub %in% st$subcluster) lab <- paste0(lab, " · ", st$nearest_CM_subtype[match(sub, st$subcluster)])
  if (!is.null(sm) && sub %in% sm$subcluster) lab <- paste0(lab, " (", sm$n_cells[match(sub, sm$subcluster)], " cells)")
  lab
}
VOLC_PAL <- c("up in KO" = "#c62828", "up in WT" = "#1565c0", "n.s." = "#cccccc", "sex/construct" = "#9e9e9e")
volc_pal <- function(pos, neg) setNames(c("#c62828","#1565c0","#cccccc","#9e9e9e"),
                                        c(pos, neg, "n.s.", "sex/construct"))
# add the derived columns (-log10 p, up/down class) the volcano + hover need.
# pos/neg name the two directions (default KO/WT for the precomputed DE tabs).
# lfc_cut is the colouring threshold, NOT a filter -- every gene stays on the plot and in
# the table. It defaults to 1 for continuity, but 1 is a poor default on this data: in the
# P7 KO-vs-WT contrasts almost nothing reaches |log2FC| = 1, so the volcano renders almost
# entirely grey and reads as "nothing changed". Each volcano tab exposes a slider.
DE_LFC_CUT <- 1
de_annot <- function(d, pos = "up in KO", neg = "up in WT", lfc_cut = DE_LFC_CUT) {
  d$neglogp <- -log10(pmax(d$pvalue, 1e-300))
  d$class <- ifelse(d$confounder, "sex/construct",
              ifelse(abs(d$log2FoldChange) >= lfc_cut, ifelse(d$log2FoldChange > 0, pos, neg), "n.s."))
  d$class <- factor(d$class, levels = c(pos, neg, "n.s.", "sex/construct"))
  d
}
# DE volcano, static ggplot twin of de_volcano_ly: same annotation, colours and
# threshold lines, no hover/click. The interactive volcano stays on screen; this
# is what the vector downloads and the Figure Studio handoff export.
de_volcano <- function(d, ttl, pos = "up in KO", neg = "up in WT",
                       xlab = "log2 fold change (KO / WT)", lfc_cut = DE_LFC_CUT,
                       highlight = NULL) {
  validate(need(!is.null(d) && nrow(d), "No DE results for this selection (cluster too small / unbalanced)."))
  d <- de_annot(d, pos, neg, lfc_cut)
  p <- ggplot(d, aes(log2FoldChange, neglogp, color = class)) +
    geom_point(size = 1.1, alpha = .6, shape = 16) +
    scale_color_manual(values = volc_pal(pos, neg), drop = FALSE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey60") +
    theme_minimal(base_size = 13) +
    labs(x = xlab, y = "-log10 p (ranking only, n=1)", color = NULL, title = ttl)
  if (!is.null(highlight) && nzchar(highlight) && highlight %in% d$gene) {
    hd <- d[match(highlight, d$gene), , drop = FALSE]
    p <- p + geom_point(data = hd, shape = 21, size = 4.5, stroke = 1.1, fill = NA, color = "#111") +
      geom_text(data = hd, aes(label = gene), color = "#111", fontface = "bold",
                size = 3.2, vjust = -1.2)
  }
  p
}
# interactive plotly volcano: hover shows gene/stats, click emits the gene via
# customdata (captured by event_data(source = source_id)) to drive the DE table.
de_volcano_ly <- function(d, ttl, source_id, pos = "up in KO", neg = "up in WT",
                           xlab = "log2 fold change (KO / WT)", highlight = NULL,
                           lfc_cut = DE_LFC_CUT) {
  validate(need(!is.null(d) && nrow(d), "No DE results for this selection (cluster too small / unbalanced)."))
  d <- de_annot(d, pos, neg, lfc_cut)
  d$hover <- sprintf(
    "<b>%s</b><br>logFC: %.2f<br>-log10 p: %.2f<br>padj: %.2g<br>%s",
    d$gene, d$log2FoldChange, d$neglogp, d$padj, as.character(d$class))
  p <- plot_ly(d, x = ~log2FoldChange, y = ~neglogp, color = ~class, colors = volc_pal(pos, neg),
          customdata = ~gene, text = ~hover, hovertemplate = "%{text}<extra></extra>",
          type = "scattergl", mode = "markers",
          marker = list(size = 6, opacity = 0.6, line = list(width = 0)),
          source = source_id) |>
    layout(
      title = list(text = ttl, font = list(size = 13)),
      xaxis = list(title = xlab, zeroline = FALSE),
      yaxis = list(title = "-log10 p (ranking only, n=1)", zeroline = FALSE),
      legend = list(title = list(text = ""), itemsizing = "constant"),
      shapes = lapply(c(-lfc_cut, lfc_cut), function(v) list(type = "line", x0 = v, x1 = v,
        yref = "paper", y0 = 0, y1 = 1, line = list(color = "grey60", width = 1, dash = "dot"))),
      margin = list(t = 34))
  # ring the selected gene's point (set by clicking a point OR a table row)
  if (!is.null(highlight) && nzchar(highlight) && highlight %in% d$gene) {
    hd <- d[match(highlight, d$gene), , drop = FALSE]
    p <- add_trace(p, x = hd$log2FoldChange, y = hd$neglogp, type = "scattergl", mode = "markers",
                   marker = list(size = 15, color = "rgba(0,0,0,0)", line = list(color = "#111", width = 3)),
                   name = "selected", showlegend = FALSE, hoverinfo = "skip", inherit = FALSE)
  }
  event_register(p, "plotly_click")
}
# ordered/filtered DE table; returns a data frame (rendered by DT). No row cap so
# every volcano point has a corresponding table row for click->highlight.
DE_PAGELEN <- 25
DE_DISP <- c(gene = "gene", log2FoldChange = "log2FC", neglog10p = "-log10(p)",
             baseMean = "baseMean", pvalue = "p", padj = "padj")   # data col -> header label
de_table <- function(d, search = "") {
  validate(need(!is.null(d) && nrow(d), "No DE table for this selection."))
  d <- d[order(-abs(d$log2FoldChange)), ]
  d$neglog10p <- -log10(pmax(d$pvalue, 1e-300))
  cols <- intersect(c("gene","log2FoldChange","neglog10p","baseMean","pvalue","padj"), names(d))
  d <- d[, cols]
  if (nzchar(search)) d <- d[grepl(search, d$gene, ignore.case = TRUE), ]
  d
}
# shared DT renderer for the DE tables (single-row select, scroll body)
de_datatable <- function(df, scroll = "380px") {
  disp <- ifelse(names(df) %in% names(DE_DISP), DE_DISP[names(df)], names(df))
  opts <- list(pageLength = DE_PAGELEN, dom = "ftip", order = list())
  if (!is.null(scroll)) { opts$scrollY <- scroll; opts$scrollCollapse <- TRUE }  # NULL = natural height, paginated
  DT::datatable(df, rownames = FALSE, selection = "single", colnames = unname(disp),
    options = opts, class = "compact stripe hover") |>
    DT::formatSignif(intersect(c("log2FoldChange","neglog10p","baseMean","pvalue","padj"), names(df)), 3)
}
# One slider definition for all four volcanoes. 0 colours every significant gene.
volc_lfc_ui <- function(id)
  sliderInput(id, "Colour genes at |log2FC| ≥", 0, 3, DE_LFC_CUT, 0.05)
# optionally drop the sex/construct confounder genes (Xist, Y-genes, ROSA26) from a
# DE frame before it reaches a volcano/table — n=1 makes these dominate the contrast.
drop_conf <- function(d, hide) {
  if (isTRUE(hide) && !is.null(d) && "confounder" %in% names(d)) d <- d[!d$confounder, , drop = FALSE]
  d
}
# small "Selected <gene>" banner above the table + a link to clear the selection
pick_banner <- function(gene, clear_id) {
  if (is.null(gene) || !nzchar(gene)) return(NULL)
  div(style = "margin-bottom:6px; font-size:13px",
      span(HTML(paste0("Selected <b>", gene, "</b> &middot; "))),
      actionLink(clear_id, "clear"))
}
# info card for a picked gene, from the bundled app$geneinfo table (GI). All data
# is precomputed (no runtime network calls); external links go to the full records.
gene_info_card <- function(gene, close_id = NULL) {
  if (is.null(gene) || !nzchar(gene)) return(NULL)
  info <- if (!is.null(GI) && gene %in% rownames(GI)) as.list(GI[gene, ]) else NULL
  has  <- function(x) !is.null(x) && length(x) && !is.na(x) && nzchar(x)
  lnk  <- function(href, label) tags$a(label, href = href, target = "_blank", style = "margin-right:14px")
  links <- list()
  if (!is.null(info)) {
    if (has(info$entrez))  links <- c(links, list(lnk(paste0("https://www.ncbi.nlm.nih.gov/gene/", info$entrez), "NCBI Gene")))
    if (has(info$ensembl)) links <- c(links, list(lnk(paste0("https://www.ensembl.org/Mus_musculus/Gene/Summary?g=", info$ensembl), "Ensembl")))
    if (has(info$mgi))     links <- c(links, list(lnk(paste0("https://www.informatics.jax.org/marker/", info$mgi), "MGI")))
  }
  links <- c(links, list(lnk(paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", toupper(gene)), "GeneCards")))
  name  <- if (!is.null(info) && has(info$name)) info$name else "(name not found)"
  bits  <- character(0)
  if (!is.null(info)) {
    if (has(info$type)) bits <- c(bits, info$type)
    if (has(info$chr))  bits <- c(bits, paste0("chr ", info$chr,
        if (has(info$start)) paste0(":", info$start, "-", info$end) else ""))
  }
  summ <- if (!is.null(info) && has(info$summary)) info$summary else
          "No functional summary available for this gene (see links)."
  card(class = "mt-2",
    card_header(class = "d-flex justify-content-between align-items-center",
      HTML(paste0("<b>", gene, "</b> &mdash; ", name)),
      if (!is.null(close_id)) actionLink(close_id, HTML("&times;"), title = "close",
        style = "font-size:20px;line-height:1;color:#888;text-decoration:none")),
    card_body(
      if (length(bits)) tags$p(tags$small(paste(bits, collapse = " · ")), style = "color:#666;margin:0 0 4px"),
      if (!is.null(info) && has(info$alias)) tags$p(tags$small(paste0("Aliases: ", info$alias)), style = "color:#666;margin:0 0 6px"),
      tags$p(summ, style = "font-size:13px;margin-bottom:8px"),
      tags$div(links)))
}
# ggplot heatmap -> plotly with hover (gene / group / value)
ggheat <- function(p) {
  ggplotly(p, tooltip = c("x", "y", "fill")) |> layout(margin = list(t = 40))
}

# ---- enrichment helpers (precomputed tables in ENR; reuse VOLC_PAL for KO/WT) --
# cell-type values are matched EXACTLY: "Cardiomyocyte" and "Cardiomyocyte(cardiac-subset)"
# are distinct analyses and must not be merged.
enr_celltypes <- function() sort(unique(c(ENR$gsea$celltype, ENR$go$celltype, ENR$tf$celltype)))
enr_dt <- function(df, scroll = "320px") {
  validate(need(!is.null(df) && nrow(df), "No results for this selection."))
  DT::datatable(df, rownames = FALSE, selection = "none",
    options = list(pageLength = 15, scrollY = scroll, scrollX = TRUE,
                   scrollCollapse = TRUE, dom = "ftip"),
    class = "compact stripe hover")
}
# ---- df-based enrichment plot cores (shared by the cell-type tab and the
# per-subcluster tab so both render identically) --------------------------------
# up_lab/down_lab default to the KO-vs-WT wording the existing call sites use.
# They are parameters because the four-group contrasts include WT: P0 vs P7,
# where a legend reading "up in KO" is not a style question -- it is wrong.
gsea_barplot_gg <- function(d, ttl, topn = 20, up_lab = "up in KO", down_lab = "up in WT",
                            base_size = 12, pal_choice = NULL, label_chars = 46L,
                            axis_scale = 1) {
  validate(need(!is.null(d) && nrow(d), "No GSEA results for this selection."))
  d <- head(d[order(-abs(d$NES)), ], topn)
  d$dir <- ifelse(d$NES > 0, up_lab, down_lab)
  # hover keeps the full name; only the axis label is shortened
  d$full <- d$pathway
  lab <- shorten_lab(d$pathway, label_chars)         # once: values and levels must agree
  d$pathway <- factor(lab, levels = lab[order(d$NES)])
  pal <- updown_pal(pal_choice)
  ggplot(d, aes(NES, pathway, fill = dir,
        text = paste0(full, "<br>NES: ", NES, "<br>padj: ", padj, "<br>size: ", size))) +
    geom_col() + geom_vline(xintercept = 0, color = "grey60") +
    scale_fill_manual(values = setNames(pal, c(up_lab, down_lab))) +
    theme_minimal(base_size = base_size) +
    theme(axis.text = element_text(size = rel(axis_scale))) +
    labs(x = sprintf("NES (>0 enriched toward %s)", up_lab), y = NULL, fill = NULL, title = ttl)
}
gsea_barplot_df <- function(d, ttl, topn = 20, up_lab = "up in KO", down_lab = "up in WT", ...)
  ggplotly(gsea_barplot_gg(d, ttl, topn, up_lab, down_lab, ...), tooltip = "text") |>
    layout(margin = list(l = 0, t = 40))
# colour_by = "padj" keeps the significance ramp; "none" drops the colour aesthetic and
# its legend entirely and draws every point in one colour. Worth having as an option
# rather than a default: the terms are already ordered and filtered by significance, so on
# a projector the gradient often buys nothing and the legend costs width the panel needs.
go_dotplot_gg <- function(d, ttl, topn = 20, base_size = 11, pal_choice = NULL,
                          label_chars = 46L, axis_scale = 1, colour_by = c("padj", "none")) {
  colour_by <- match.arg(colour_by)
  validate(need(!is.null(d) && nrow(d), "No GO BP results for this selection."))
  d <- head(d[order(d$p.adjust), ], topn)
  d$full <- d$Description
  lab <- shorten_lab(d$Description, label_chars)     # once: values and levels must agree
  d$Description <- factor(lab, levels = rev(lab))
  tip <- paste0(d$full, "<br>fold: ", d$FoldEnrichment, "<br>padj: ", d$p.adjust,
                "<br>genes: ", d$Count)             # padj stays in the hover either way
  p <- if (identical(colour_by, "none")) {
    # the darkest end of the chosen ramp, so the palette selector still means something
    # and the points stay the most visible thing on the panel
    ggplot(d, aes(FoldEnrichment, Description, size = Count, text = tip)) +
      geom_point(colour = cont_pal(pal_choice)[1])
  } else {
    ggplot(d, aes(FoldEnrichment, Description, size = Count, color = p.adjust, text = tip)) +
      geom_point() +
      # colours[1] is the low-p.adjust end, so the most significant terms are the darkest.
      # This replaces scale_color_viridis_c(option="magma", direction=-1), which did the
      # opposite and made the significant points nearly invisible.
      scale_color_gradientn(colours = cont_pal(pal_choice), na.value = "grey85")
  }
  # labs(color=) only when a colour aesthetic exists, else ggplot warns
  # "Ignoring unknown labels" on every single render.
  p <- p + theme_minimal(base_size = base_size) +
    theme(axis.text = element_text(size = rel(axis_scale))) +
    labs(x = "fold enrichment", y = NULL, size = "genes", title = ttl)
  if (identical(colour_by, "padj")) p <- p + labs(color = "padj")
  p
}
go_dotplot_df <- function(d, ttl, topn = 20, ...)
  ggplotly(go_dotplot_gg(d, ttl, topn, ...), tooltip = "text") |> layout(margin = list(l = 0, t = 40))
# static "all subclusters at once" overviews (faceted; renderPlot, low-memory) --
# "All clusters" enrichment view: one full-size plot per subcluster, two per row
# (server registers a renderPlot per subcluster id "cm_<kind>_all_<CMn>").
# `per_row` is a real control, not a constant: a GO dot plot two-up gets ~450 px of width
# inside this sidebar layout, and 40-character term names then take most of it. One-up is
# the readable default; two-up stays available because comparing subclusters side by side
# is the point of this view, and the labels are truncated hard to make it viable.
cm_enr_grid <- function(kind, per_row = 1L, height = "420px") {
  subs <- cm_subs("0.2")
  outs <- lapply(subs, function(cl) plotOutput(paste0("cm_", kind, "_all_", cl), height = height))
  w <- if (per_row >= 2L) 6L else 12L
  rows <- lapply(seq(1, length(outs), by = per_row), function(i)
    fluidRow(lapply(seq(i, min(i + per_row - 1L, length(outs))), function(j) column(w, outs[[j]]))))
  do.call(tagList, rows)
}
# per-subcluster enrichment slice (precomputed into ENR$sub by build_subcluster_enrichment.R)
enr_sub_df <- function(kind, sub) {
  d <- ENR$sub[[kind]]
  validate(need(!is.null(d), "Per-subcluster enrichment isn't in this data build — run build_subcluster_enrichment.R and redeploy."))
  d[d$subcluster == sub, , drop = FALSE]
}

enr_gsea <- function(ct, tp) { d <- ENR$gsea
  validate(need(!is.null(d), "GSEA results are not in this data build."))
  d[d$celltype == ct & d$timepoint == tp, , drop = FALSE] }
enr_gsea_plot <- function(ct, tp, topn = 20, ...)
  gsea_barplot_df(enr_gsea(ct, tp), paste0("GSEA — ", ct, " ", tp), topn, ...)
# Split df-from-widget so downloads can reuse the frame (same for the three below).
enr_gsea_table_df <- function(ct, tp) {
  d <- enr_gsea(ct, tp)
  d[order(d$padj), intersect(c("pathway","NES","padj","size","leadingEdge"), names(d))]
}
enr_gsea_table <- function(ct, tp) enr_dt(enr_gsea_table_df(ct, tp))
enr_go <- function(ct, tp) { d <- ENR$go
  validate(need(!is.null(d), "GO results are not in this data build."))
  d[d$celltype == ct & d$timepoint == tp, , drop = FALSE] }
enr_go_plot <- function(ct, tp, topn = 20, ...)
  go_dotplot_df(enr_go(ct, tp), paste0("GO BP enriched in KO-up genes — ", ct, " ", tp), topn, ...)
enr_go_table_df <- function(ct, tp) {
  d <- enr_go(ct, tp)
  d[order(d$p.adjust), intersect(c("ID","Description","FoldEnrichment","p.adjust","Count","geneID"), names(d))]
}
enr_go_table <- function(ct, tp) enr_dt(enr_go_table_df(ct, tp))
# E2F-family regulon activity (KO - WT) across cell type x timepoint
enr_e2f_heat_gg <- function() {
  d <- tabs$e2f_regulon
  validate(need(!is.null(d) && nrow(d), "No E2F regulon activity table."))
  d$col <- paste0(d$celltype, " ", d$timepoint)
  p <- ggplot(d, aes(col, source, fill = KO_minus_WT,
        text = paste0(source, "<br>", col, "<br>KO-WT: ", round(KO_minus_WT, 3)))) +
    geom_tile(color = "grey92") +
    scale_fill_gradient2(low = "#1565c0", mid = "white", high = "#c62828", midpoint = 0) +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "KO - WT", title = "E2F-family regulon activity (KO - WT)")
  p
}
enr_e2f_heat <- function() ggplotly(enr_e2f_heat_gg(), tooltip = "text") |> layout(margin = list(t = 40))
# top TFs by |KO - WT| activity for a cell type (from ENR$tf)
enr_tf_top_gg <- function(ct, topn = 20, base_size = 11, pal_choice = NULL, axis_scale = 1) {
  d <- ENR$tf
  validate(need(!is.null(d), "TF activity is not in this data build."))
  d <- d[d$celltype == ct, , drop = FALSE]
  validate(need(nrow(d), "No TF activity for this cell type."))
  ko <- setNames(d$mean_activity[d$genotype == "KO"], d$source[d$genotype == "KO"])
  wt <- setNames(d$mean_activity[d$genotype == "WT"], d$source[d$genotype == "WT"])
  src <- intersect(names(ko), names(wt))
  validate(need(length(src), "TF activity needs both KO and WT."))
  w <- data.frame(source = src, KO = ko[src], WT = wt[src]); w$diff <- w$KO - w$WT
  w <- head(w[order(-abs(w$diff)), ], topn)
  w$source <- factor(w$source, levels = w$source[order(w$diff)])
  p <- ggplot(w, aes(diff, source, fill = diff > 0,
        text = paste0(source, "<br>KO: ", round(KO,3), "<br>WT: ", round(WT,3), "<br>KO-WT: ", round(diff,3)))) +
    geom_col() + geom_vline(xintercept = 0, color = "grey60") +
    scale_fill_manual(values = setNames(updown_pal(pal_choice), c("TRUE", "FALSE")), guide = "none") +
    theme_minimal(base_size = base_size) +
    theme(axis.text = element_text(size = rel(axis_scale))) +
    labs(x = "KO - WT activity", y = NULL, title = paste0("Top TFs by |KO-WT| — ", ct))
  p
}
enr_tf_top <- function(ct, topn = 20, ...)
  ggplotly(enr_tf_top_gg(ct, topn, ...), tooltip = "text") |> layout(margin = list(l = 0, t = 40))
# log2FC heatmap: top genes (by max |LFC| across groups) x groups, fill = KO/WT log2FC
lfc_heat <- function(de_list, topn = 22, ttl = NULL, fill_lab = "log2FC\n(KO/WT)") {
  de_list <- de_list[!vapply(de_list, is.null, logical(1))]
  validate(need(length(de_list) >= 1, "No DE tables to compare."))
  allg <- unique(unlist(lapply(de_list, `[[`, "gene")))
  M <- vapply(de_list, function(d) d$log2FoldChange[match(allg, d$gene)], numeric(length(allg)))
  rownames(M) <- allg
  score <- apply(abs(M), 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
  score[allg %in% CONF] <- NA
  top <- names(sort(score, decreasing = TRUE))[seq_len(min(topn, sum(!is.na(score))))]
  long <- expand.grid(gene = top, grp = colnames(M), stringsAsFactors = FALSE)
  long$lfc <- M[cbind(match(long$gene, rownames(M)), match(long$grp, colnames(M)))]
  long$gene <- factor(long$gene, levels = rev(top))
  ggplot(long, aes(grp, gene, fill = lfc)) + geom_tile(color = "grey92") + div_scale +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = fill_lab, title = ttl)
}

# ---- plotly UMAP helpers: distinct colours, hover tooltips, hover-to-highlight ----
PAL <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFD92F","#A65628","#F781BF",
         "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666",
         "#1F78B4","#33A02C","#FB9A99","#FDBF6F","#CAB2D6","#B15928","#6A3D9A","#B2DF8A")
pal_for <- function(levs) setNames(PAL[((seq_along(levs) - 1) %% length(PAL)) + 1], levs)
# ---- publication aesthetics: palettes + label-rename helpers --------------------
# colourblind-safe discrete palettes (hardcoded hex -> no extra package deps).
PALETTES_DISCRETE <- list(
  "Default"   = PAL,
  "Okabe-Ito" = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000"),
  "Set2"      = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
  "Dark2"     = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
  "Viridis"   = NULL)   # NULL -> generated per-n via grDevices::hcl.colors below
# ---- continuous palettes for significance-coloured plots --------------------
# Written as explicit low->high ramps rather than viridis option+direction, because the
# direction argument is what went wrong here: go_dotplot_gg had
# scale_color_viridis_c(option="magma", direction=-1) on p.adjust, which put the PALE
# YELLOW end on the SMALLEST p.adjust -- the most significant terms were drawn in the
# least visible colour. Reported from a projector, but it was wrong on any display.
#
# Invariant for everything in this list: element 1 is the colour for the LOW end of the
# mapped variable. Since these scales carry p.adjust, element 1 is what "most significant"
# gets, so element 1 is always the darkest. Encoding that in the data makes it impossible
# to reintroduce the bug by flipping a direction flag.
PALETTES_CONTINUOUS <- list(
  "Blue (significant = dark)"     = c("#08306b", "#2171b5", "#6baed6", "#c6dbef"),
  "Red (significant = dark)"      = c("#67000d", "#cb181d", "#fb6a4a", "#fcbba1"),
  "Magma (significant = dark)"    = c("#000004", "#51127c", "#b63679", "#fb8861"),
  "Viridis (significant = dark)"  = c("#440154", "#31688e", "#35b779", "#8fd744"),
  "Grey (print-safe)"             = c("#111111", "#555555", "#999999", "#cccccc"))
cont_pal <- function(choice = NULL) {
  if (is.null(choice) || !choice %in% names(PALETTES_CONTINUOUS))
    choice <- names(PALETTES_CONTINUOUS)[1]
  PALETTES_CONTINUOUS[[choice]]
}
# Up/down pairs for the diverging bar charts, high-contrast on a projector. As above, no
# pale end: both directions have to be readable, not just the one you happen to care about.
PALETTES_UPDOWN <- list(
  "Red / blue"      = c("#b2182b", "#2166ac"),
  "Orange / purple" = c("#b35806", "#542788"),
  "Okabe-Ito"       = c("#D55E00", "#0072B2"),
  "Grey (print-safe)" = c("#252525", "#969696"))
updown_pal <- function(choice = NULL) {
  if (is.null(choice) || !choice %in% names(PALETTES_UPDOWN))
    choice <- names(PALETTES_UPDOWN)[1]
  PALETTES_UPDOWN[[choice]]
}
# Long GO/pathway names are the other half of the "axis labels too large" complaint: at
# 440 px a 90-character term eats the plot. Truncate on a word boundary where possible.
shorten_lab <- function(x, chars = 46L) {
  x <- as.character(x)
  out <- ifelse(nchar(x) <= chars, x,
                paste0(sub("\\s+\\S*$", "", substr(x, 1, chars)), "\u2026"))
  # Two long terms sharing a prefix truncate to the same string, and a factor cannot carry
  # duplicate levels -- "GO BP regulation of ..." twice is a real case, not a corner one.
  # sep = " #" so the disambiguator reads as one, rather than looking like part of the term
  # name. The full name is preserved in the hover text and in the downloaded table.
  make.unique(out, sep = " #")
}

# named colour vector over `levs` for the chosen palette (recycles like pal_for)
disc_pal <- function(levs, choice = "Default") {
  if (is.null(choice) || !choice %in% names(PALETTES_DISCRETE)) choice <- "Default"
  pal <- if (identical(choice, "Viridis")) grDevices::hcl.colors(max(length(levs), 1), "viridis")
         else PALETTES_DISCRETE[[choice]]
  setNames(pal[((seq_along(levs) - 1) %% length(pal)) + 1], levs)
}
# map a character/level vector through a rename list; unmapped values pass through.
relab <- function(x, map = list()) {
  x <- as.character(x)
  if (!length(map)) return(x)
  hit <- x %in% names(map)
  x[hit] <- unlist(map[x[hit]]); x
}
.ax       <- list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "")
.exprsc   <- list(c(0,"#eeeeee"), c(0.45,"#fec44f"), c(0.75,"#fc4e2a"), c(1,"#800026"))
.umap_mar <- list(l = 0, r = 0, b = 0, t = 26)
centroids <- function(df, size = 12, map = list()) { a <- aggregate(cbind(UMAP1, UMAP2) ~ val, df, median)
  lapply(seq_len(nrow(a)), function(i) list(x = a$UMAP1[i], y = a$UMAP2[i], text = relab(a$val[i], map),
    showarrow = FALSE, font = list(size = size, color = "#111"), bgcolor = "rgba(255,255,255,0.55)")) }

# single-panel categorical: hover a cell -> highlight ALL cells of that cluster (crosstalk),
# double-click a legend entry to isolate; centroid labels make clusters easy to read.
# labels/labsize/legend/pal/map are publication-figure controls (defaults = original look).
umap_cat <- function(df, colvar, ttl = NULL, psize = 4.5, labels = TRUE,
                     labsize = 12, legend = TRUE, pal_choice = "Default", map = list()) {
  df$val <- factor(df[[colvar]])
  if (length(map)) levels(df$val) <- relab(levels(df$val), map)   # rename legend + centroid labels
  cols <- disc_pal(levels(df$val), pal_choice)
  grp <- paste0("umaphl_", make.names(colvar))                 # named crosstalk group to clear from JS
  sd <- crosstalk::SharedData$new(df, key = ~val, group = grp)
  p <- plot_ly(sd, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
               color = ~val, colors = cols,
               marker = list(size = psize, line = list(width = 0)), text = ~val, hoverinfo = "text") |>
    layout(xaxis = .ax, yaxis = .ax, margin = .umap_mar, title = list(text = ttl, font = list(size = 13)),
           showlegend = legend, legend = list(itemsizing = "constant"),
           annotations = if (labels) centroids(df, size = labsize) else NULL)
  # hover a cell -> highlight its whole subtype (crosstalk). crosstalk's own `off`
  # can't fire on un-hover, so we clear the selection ourselves on unhover and on a
  # click that hits no point -> all cells return to normal.
  p <- highlight(p, on = "plotly_hover", off = "plotly_deselect",
                 opacityDim = 0.10, persistent = FALSE)
  htmlwidgets::onRender(p, sprintf("
    function(el) {
      var clear = function() {
        try { crosstalk.group('%s').var('selection').set(null); } catch(e) {}
      };
      el.on('plotly_unhover', clear);
      el.on('plotly_click', function(d) {
        if (!d || !d.points || d.points.length === 0) clear();
      });
    }", grp))
}
# single-panel continuous (gene expression / scores)
umap_cont <- function(df, val, ttl = NULL, psize = 4.5, legend = TRUE) {
  df$val <- as.numeric(val); df <- df[order(df$val, na.last = FALSE), ]
  plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
          marker = list(size = psize, color = ~val, colorscale = .exprsc, showscale = legend,
                        colorbar = list(title = ""), line = list(width = 0)),
          text = ~round(val, 2), hoverinfo = "text") |>
    layout(xaxis = .ax, yaxis = .ax, margin = .umap_mar, title = list(text = ttl, font = list(size = 13)))
}
# side-by-side split panels (legend / colourbar shown once). pal_choice/map rename + recolour.
# nrows > 1 lays the panels out on a grid (row-major), e.g. genotype x timepoint as a 2x2.
umap_split <- function(df, colvar, splitvar, gene = NULL, continuous = FALSE, psize = 4,
                       legend = TRUE, pal_choice = "Default", map = list(), nrows = 1) {
  levs_sp  <- if (is.factor(df[[splitvar]])) levels(droplevels(factor(df[[splitvar]]))) else sort(unique(as.character(df[[splitvar]])))
  val_levs <- if (!continuous) relab(levels(factor(df[[colvar]])), map) else NULL
  cols     <- if (!continuous) disc_pal(val_levs, pal_choice) else NULL
  plts <- lapply(seq_along(levs_sp), function(j) {
    d <- df[as.character(df[[splitvar]]) == levs_sp[j], ]
    if (continuous) {
      d$val <- if (!is.null(gene)) expr_vec(gene, d$cell) else d[[colvar]]; d <- d[order(d$val, na.last = FALSE), ]
      plot_ly(d, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
              marker = list(size = psize, color = ~val, colorscale = .exprsc, showscale = (j == 1 && legend),
                            colorbar = list(title = ""), line = list(width = 0)),
              text = ~round(val, 2), hoverinfo = "text") |> layout(xaxis = .ax, yaxis = .ax)
    } else {
      d$val <- factor(relab(d[[colvar]], map), levels = val_levs); p <- plot_ly()
      for (i in seq_along(val_levs)) { di <- d[d$val == val_levs[i], ]
        p <- add_trace(p, data = di, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
                       name = val_levs[i], legendgroup = val_levs[i], showlegend = (j == 1 && legend),
                       marker = list(size = psize, color = cols[i], line = list(width = 0)),
                       text = val_levs[i], hoverinfo = "text") }
      p |> layout(xaxis = .ax, yaxis = .ax)
    }
  })
  ncol <- ceiling(length(levs_sp) / nrows)
  anns <- lapply(seq_along(levs_sp), function(j) {
    col <- (j - 1) %% ncol; row <- (j - 1) %/% ncol           # row-major (matches subplot fill order)
    list(text = paste0(labof(splitvar), ": ", levs_sp[j]),
         x = (col + 0.5) / ncol, y = 1 - row / nrows,
         xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 13)) })
  subplot(plts, nrows = nrows, shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE) |>
    layout(margin = .umap_mar, showlegend = legend, legend = list(itemsizing = "constant"), annotations = anns)
}
# Static ggplot twin of the umap_cat / umap_cont / umap_split plotly views: same
# palettes, centroid labels and split layout, drawn with geom_point so the UMAP
# has a vector (PDF/SVG) export and a Figure Studio handoff. No hover/highlight.
umap_gg <- function(df, colvar, splitvar = NULL, gene = NULL, continuous = FALSE,
                    ttl = NULL, psize = 4.5, labels = TRUE, labsize = 12, legend = TRUE,
                    pal_choice = "Default", map = list(), nrows = 1) {
  pt <- psize / 3.4                       # plotly px -> ggplot mm, matched by eye
  if (continuous) {
    df$val <- if (!is.null(gene)) as.numeric(expr_vec(gene, df$cell)) else as.numeric(df[[colvar]])
    df <- df[order(df$val, na.last = FALSE), ]
  } else {
    df$val <- factor(df[[colvar]])
    if (length(map)) levels(df$val) <- relab(levels(df$val), map)
  }
  if (!is.null(splitvar)) {
    levs_sp <- if (is.factor(df[[splitvar]])) levels(droplevels(factor(df[[splitvar]])))
               else sort(unique(as.character(df[[splitvar]])))
    df$panel <- factor(paste0(labof(splitvar), ": ", as.character(df[[splitvar]])),
                       levels = paste0(labof(splitvar), ": ", levs_sp))
  }
  p <- ggplot(df, aes(UMAP1, UMAP2, color = val)) + geom_point(size = pt, shape = 16) +
    theme_umap + labs(title = ttl, color = NULL)
  p <- if (continuous)
    p + scale_color_gradientn(colours = vapply(.exprsc, `[[`, "", 2),
                              values = as.numeric(vapply(.exprsc, `[[`, "", 1)),
                              na.value = "grey92", name = NULL)
  else
    p + scale_color_manual(values = disc_pal(levels(df$val), pal_choice), drop = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 3)))
  if (!continuous && labels && is.null(splitvar)) {
    a <- aggregate(cbind(UMAP1, UMAP2) ~ val, df, median)
    p <- p + geom_label(data = a, aes(label = val), color = "#111", size = labsize / 2.845,
                        fill = scales::alpha("white", 0.55), linewidth = 0, label.padding = unit(0.1, "lines"))
  }
  if (!is.null(splitvar)) p <- p + facet_wrap(~panel, nrow = nrows)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# Three-panel UMAP for the object-mode diagnostic. Not umap_gg(): its splitvar path
# prefixes every facet with labof(splitvar), which reads badly when the facet IS the
# variant name. theme_umap and disc_pal are reused so it matches the rest of the app.
objtest_gg <- function(d, rc, pal_choice = "Default", psize = 0.35) {
  lv <- sort(unique(suppressWarnings(as.integer(as.character(d[[rc]])))))
  d$val   <- factor(paste0("CM", d[[rc]]), levels = paste0("CM", lv))
  d$panel <- factor(unname(CMTEST$labels[d$variant]), levels = unname(CMTEST$labels))
  ggplot(d, aes(UMAP1, UMAP2, color = val)) +
    geom_point(size = psize, shape = 16) + theme_umap +
    facet_wrap(~panel, nrow = 1) + labs(color = NULL) +
    scale_color_manual(values = disc_pal(levels(d$val), pal_choice), drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3)))
}

# Faceted UMAP for the PC-dimension sweep. The facet strip carries the cumulative PCA
# variance for each cut, because "dims 1:10" on its own does not tell you whether the
# discarded components held anything.
pcdims_gg <- function(g, colvar, pal_choice = "Default", psize = 0.3) {
  d <- g$percell
  d$val <- if (colvar == "cluster") {
    lv <- sort(unique(suppressWarnings(as.integer(as.character(d$cluster)))))
    factor(paste0(g$prefix, d$cluster), levels = paste0(g$prefix, lv))
  } else factor(d[[colvar]])
  # "of the top-50 PC variance", never "of PCA variance": Stdev() only returns the PCs
  # that were computed, so the 50-PC figure is 100% by construction. Labelling it as a
  # share of total variance would invite reading that as "50 PCs capture everything".
  strip <- function(x) sprintf("dims 1:%d \u2014 %.1f%% of the top-50 PC variance",
                               x, g$varpct[as.character(x)])
  d$panel <- factor(strip(d$dims), levels = strip(g$dims))
  ggplot(d, aes(UMAP1, UMAP2, color = val)) +
    geom_point(size = psize, shape = 16) + theme_umap +
    facet_wrap(~panel, nrow = 1) + labs(color = NULL) +
    scale_color_manual(values = disc_pal(levels(d$val), pal_choice), drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3)))
}

# ---- publication "Figure options": reusable control block + export wiring -------
# emits a namespaced (prefix_*) control set; `export` picks the download UI (vector
# ggsave for ggplot panels, camera-button PNG for the WebGL UMAP).
# `palette` picks WHICH registry the selector offers: TRUE/"discrete" for categorical
# fills, "continuous" for a significance ramp, "updown" for a two-direction bar chart,
# FALSE for none. `export = "none"` emits options WITHOUT download buttons, for a plot that
# already carries dl_fig_ui() above it -- emitting both would duplicate the input ids.
figure_controls <- function(prefix, export = c("ggplot","umap","none"), base = TRUE,
                            palette = TRUE, rename = TRUE, labels = FALSE,
                            default_base = 13, axis = FALSE, labelchars = FALSE,
                            colourby = FALSE) {
  export <- match.arg(export); p <- function(s) paste0(prefix, "_", s)
  pal_choices <- if (isFALSE(palette)) NULL
    else switch(if (isTRUE(palette)) "discrete" else palette,
                discrete = names(PALETTES_DISCRETE),
                continuous = names(PALETTES_CONTINUOUS),
                updown = names(PALETTES_UPDOWN), NULL)
  tagList(
    textInput(p("title"), "Title (blank = default)", ""),
    if (labels) checkboxInput(p("labels"), "Show cluster labels", TRUE),
    if (labels) sliderInput(p("labelsize"), "Label font size", 8, 28, 12, 1),
    if (base)   sliderInput(p("basesize"), "Base font size", 8, 24, default_base, 1),
    # Axis text separately from base size: the projector complaint was specifically that
    # the AXIS labels were too large, and scaling everything down to fix them shrinks the
    # title and legend too.
    if (axis)   sliderInput(p("axisscale"), "Axis label size (relative)", 0.5, 1.6, 1, 0.05),
    # The other half of "labels too large": a 90-character GO term eats the panel. Hover
    # and the download table keep the full name.
    if (labelchars) sliderInput(p("labelchars"), "Truncate long term names at", 20, 90, 46, 2),
    if (colourby) radioButtons(p("colourby"), "Colour points by",
                               c("Adjusted p-value" = "padj", "Nothing (single colour)" = "none"),
                               selected = "padj"),
    checkboxInput(p("legend"), "Show legend", TRUE),
    if (!is.null(pal_choices)) selectInput(p("palette"), "Palette", pal_choices),
    if (rename) tagList(tags$small("Rename categories (double-click a label cell):"),
                        DTOutput(p("renametab"))),
    hr(),
    # The download buttons are always visible. They used to sit inside the
    # "Custom export options" conditionalPanel, which meant that unless you
    # thought to tick a box labelled as being about options, the app looked like
    # it could not export figures at all. Only the size/DPI inputs are optional --
    # they have sensible defaults and most people never touch them.
    if (export == "ggplot")
      div(style = "display:flex;gap:6px;margin-bottom:6px;flex-wrap:wrap",
        downloadButton(p("dl_pdf"), "PDF", class = "btn-sm btn-outline-secondary"),
        downloadButton(p("dl_svg"), "SVG", class = "btn-sm btn-outline-secondary"),
        downloadButton(p("dl_png"), "PNG", class = "btn-sm btn-outline-secondary"),
        studio_btn(prefix)),
    if (export != "none") checkboxInput(p("export_on"), "Custom export size", FALSE),
    if (export != "none") conditionalPanel(sprintf("input.%s_export_on", prefix),
      if (export == "ggplot") tagList(
        div(style = "display:flex;gap:6px",
          numericInput(p("w"), "W (in)", 7, 1, 20, 0.5),
          numericInput(p("h"), "H (in)", 5, 1, 20, 0.5),
          numericInput(p("dpi"), "DPI", 300, 72, 600, 1)))
      else tagList(
        div(style = "display:flex;gap:6px",
          numericInput(p("w"), "W (px)", 1200, 300, 4000, 50),
          numericInput(p("h"), "H (px)", 900, 300, 4000, 50)),
        sliderInput(p("scale"), "Resolution scale", 1, 4, 2, 0.5),
        helpText("Camera icon on the plot saves a PNG at this size/scale. ",
                 "Leave unchecked to use the camera's default download."))))
}
# editable 2-column rename table (category read-only, label editable)
rename_table <- function(levs, map = list()) {
  df <- data.frame(category = levs, label = relab(levs, map), stringsAsFactors = FALSE)
  DT::datatable(df, rownames = FALSE, selection = "none",
    editable = list(target = "cell", disable = list(columns = 0)),
    options = list(dom = "t", paging = FALSE, scrollY = "180px", ordering = FALSE),
    class = "compact stripe")
}
# apply the generic options a built ggplot can take without knowing its data:
# title override (blank = keep the plot's own) + legend toggle. Palette/rename/base
# size are applied inside each plot's build reactive (they need the data levels).
apply_fig_opts <- function(p, prefix, input) {
  ttl <- input[[paste0(prefix, "_title")]]
  if (!is.null(ttl) && nzchar(ttl)) p <- p + labs(title = ttl)
  if (isFALSE(input[[paste0(prefix, "_legend")]])) p <- p + theme(legend.position = "none")
  p
}
# downloadHandler factory: regenerates the SAME plot object for export via ggsave.
# opts_prefix lets a plot have its own download id while sharing another plot's
# title/legend/size controls -- several tabs show two or three figures under one
# "Figure options" accordion, and duplicating the accordion per panel would be
# worse than sharing it.
dl_ggplot <- function(prefix, plot_reactive, input, fmt, opts_prefix = prefix) {
  downloadHandler(
    filename = function() paste0(prefix, "_", Sys.Date(), ".", fmt),
    content  = function(file) {
      g <- function(s) input[[paste0(opts_prefix, "_", s)]]
      ggsave(file, apply_fig_opts(plot_reactive(), opts_prefix, input), device = fmt,
             width = g("w") %||% 7, height = g("h") %||% 5, dpi = g("dpi") %||% 300, units = "in")
    })
}
# Register all three formats for one figure, so a plot is one line to export.
# Also wires the Figure Studio handoff for the same prefix (a no-op when the
# studio isn't configured), so every exportable figure gets both for free.
register_fig <- function(output, prefix, plot_reactive, input, opts_prefix = prefix) {
  for (.f in c("pdf", "svg", "png")) local({
    f <- .f
    output[[paste0(prefix, "_dl_", f)]] <- dl_ggplot(prefix, plot_reactive, input, f, opts_prefix)
  })
  register_studio(prefix, plot_reactive, input, opts_prefix)
  invisible(NULL)
}
# Figure download buttons for a plot that has no sidebar accordion of its own
# (i.e. one that shares another prefix's options). Sits above the plot.
dl_fig_ui <- function(prefix, label = "Download figure") {
  btn <- function(f) downloadButton(paste0(prefix, "_dl_", f), toupper(f),
                                    class = "btn-sm btn-outline-secondary")
  div(class = "d-flex align-items-center gap-2 mb-2 flex-wrap",
      tags$small(class = "text-muted", label), lapply(c("pdf", "svg", "png"), btn),
      studio_btn(prefix))
}

# ---- data downloads: the table-side counterpart to dl_ggplot ----------------
# Every table in this app should be obtainable as a file. There are ~20 of them,
# so these are factories driven from a registry (TABLE_DL, in the server) rather
# than twenty hand-written handlers -- the same shape the figure exports use.
#
# Rule: the button exports the FULL underlying frame, not the display subset the
# DT shows. The display frame drops columns to fit on screen, which is exactly
# the wrong thing to do to data someone is taking away to analyse. Where the two
# differ only in columns that is a harmless superset; where the display frame is
# a genuinely different SHAPE (a wide pivot, say), pass that one instead.

# One data frame, one format. `basename` may be a string or a function, so a
# filename can carry the current selection (which is what makes a downloaded file
# still identifiable a week later, in a folder full of them).
dl_name <- function(basename) if (is.function(basename)) basename() else basename
dl_data <- function(basename, df_fn, fmt, title = NULL, notes = NULL, sheet = "data") {
  downloadHandler(
    filename = function() paste0(dl_name(basename), "_", dl_stamp(), ".", fmt),
    content  = function(file) {
      d <- tryCatch(df_fn(), error = function(e) NULL)
      if (is.null(d) || !nrow(d)) d <- data.frame(note = "No rows for the current selection.")
      if (identical(fmt, "csv")) dl_write_csv(d, file)
      else dl_write_xlsx(setNames(list(d), sheet), file, title = title, notes = notes)
    })
}

# Many data frames -> one workbook. sheets_fn returns a named list.
dl_book <- function(basename, sheets_fn, title = NULL, notes = NULL) {
  downloadHandler(
    filename = function() paste0(dl_name(basename), "_", dl_stamp(), ".xlsx"),
    content  = function(file) {
      sh <- tryCatch(sheets_fn(), error = function(e) NULL)
      if (is.null(sh) || !length(sh))
        sh <- list(data = data.frame(note = "Nothing to export for the current selection."))
      # title/notes may be functions so they can name the current selection.
      dl_write_xlsx(sh, file,
                    title = if (is.function(title)) title() else title,
                    notes = if (is.function(notes)) notes() else notes)
    })
}

# A "how this was made" block, placed directly under its figure. A <details> element, so
# it costs one line of screen until someone opens it, and it sits beside the plot rather
# than in a document that drifts out of date. `code` names FUNCTIONS AND FILES, never line
# numbers -- those rot the first time anything above them is edited.
method_note <- function(..., code = NULL, title = "How this plot was made") {
  tags$details(
    style = paste0("font-size:12.5px;margin-top:12px;border:1px solid #e6e6e6;",
                   "border-radius:6px;padding:6px 12px;background:#fafafa"),
    tags$summary(style = "cursor:pointer;font-weight:600;color:#2c3e50", title),
    div(style = "padding-top:8px;line-height:1.55", ...),
    if (!is.null(code)) div(style = "margin-top:8px;color:#555",
      HTML(paste0("<b>Code:</b> ",
                  paste(sprintf("<code>%s</code>", code), collapse = " &nbsp;·&nbsp; ")))))
}
# The button pair. Placed in the content pane directly above its table, which is
# where a reader looks for it and which works inside navset_card_tab panels that
# the sidebar cannot cleanly address.
dl_data_ui <- function(id, label = "Download table", formats = c("csv", "xlsx")) {
  btn <- function(f) downloadButton(paste0(id, "_dl_", f), toupper(f),
                                    class = "btn-sm btn-outline-secondary")
  div(class = "d-flex align-items-center gap-2 mb-2",
      tags$small(class = "text-muted", label), lapply(formats, btn))
}

# Register <id>_dl_csv and <id>_dl_xlsx for one table.
register_dl <- function(output, id, df_fn, basename, title = NULL, notes = NULL,
                        formats = c("csv", "xlsx")) {
  for (f in formats) {
    local({
      ff <- f
      output[[paste0(id, "_dl_", ff)]] <- dl_data(basename, df_fn, ff, title, notes)
    })
  }
  invisible(NULL)
}

# ---- interactive subset DEG (descriptive Wilcoxon via presto on log-norm) ----
# filters: named list colname -> selected levels (NULL/empty = all). Returns a
# logical mask over meta rows (= columns of EXPR, same order).
deg_mask <- function(filters) {
  keep <- rep(TRUE, nrow(DMETA))
  for (col in names(filters)) {
    sel <- filters[[col]]
    if (!is.null(sel) && length(sel) && col %in% names(DMETA))
      keep <- keep & as.character(DMETA[[col]]) %in% sel
  }
  keep
}
# compute descriptive DE between two groups within a cell subset.
# grpvar = meta column to split on; a/b = the two level-sets (group A vs B).
deg_compute <- function(mask, grpvar, a_levels, b_levels) {
  validate(need(!is.null(EXPR), "Expanded expression matrix not in this data build."))
  gv <- as.character(DMETA[[grpvar]])
  inA <- mask & gv %in% a_levels
  inB <- mask & gv %in% b_levels
  nA <- sum(inA); nB <- sum(inB)
  validate(need(nA >= 10 && nB >= 10,
    sprintf("Need >= 10 cells per group (A = %d, B = %d). Loosen the filters.", nA, nB)))
  cols <- which(inA | inB)
  grp  <- ifelse(inA[cols], "A", "B")
  X <- EXPR[, cols, drop = FALSE]
  res <- presto::wilcoxauc(X, grp)
  res <- res[res$group == "A", ]                       # logFC > 0 => up in group A
  data.frame(gene = res$feature, log2FoldChange = res$logFC,
             pvalue = res$pval, padj = res$padj,
             pct_A = round(res$pct_in, 1), pct_B = round(res$pct_out, 1),
             confounder = res$feature %in% CONF, n_A = nA, n_B = nB,
             stringsAsFactors = FALSE)
}

# ---- new-analysis helpers: module scores / communication / annotation check --
# All read precomputed slots (build_signature_scores.R / build_communication.R /
# build_refmap.R). Each validate()s a clear "run the builder" message so an
# un-rebuilt bundle still loads.
CELLTYPE_CHOICES <- c("All cells" = "All",
  if ("celltype" %in% names(meta)) setNames(sort(unique(as.character(meta$celltype))),
    gsub("_", " ", sort(unique(as.character(meta$celltype))))) else NULL)
CM_DEFAULT_CT <- if ("celltype" %in% names(meta) && "Cardiomyocyte" %in% meta$celltype) "Cardiomyocyte" else "All"

score_df <- function(scol, ct) {
  validate(need(scol %in% names(meta), paste0(
    "Score not in this data build — run build_signature_scores.R and redeploy.")))
  df <- meta
  if (!is.null(ct) && ct != "All" && "celltype" %in% names(df)) df <- df[as.character(df$celltype) == ct, ]
  df <- df[!is.na(df[[scol]]), ]
  validate(need(nrow(df), "No scored cells for this selection."))
  df
}
# violin+box of a score across the four groups. One axis rather than genotype +
# timepoint facets, so WT-P0/WT-P7/KO-P0/KO-P7 can be compared directly; n and the
# mean are drawn on, and the G1 stratum holds cycling composition fixed (P7 was
# cycling-enriched 4.5-5.2x relative to P0, which shifts any score that tracks cycling).
score_violin <- function(scol, ct, bs = 13, pal = "Default", stratum = "all") {
  df <- score_df(scol, ct)
  if (stratum == "G1" && "Phase" %in% names(df)) df <- df[as.character(df$Phase) == "G1", , drop = FALSE]
  validate(need(nrow(df), "No cells for this selection — the G1 stratum may be empty here."))
  df$y <- df[[scol]]
  has_tp <- "timepoint" %in% names(df)
  df$grp <- if (has_tp) factor(paste(df$genotype, df$timepoint, sep = "-"), levels = FG_GROUPS)
            else factor(df$genotype)
  df <- df[!is.na(df$grp), , drop = FALSE]
  validate(need(nrow(df), "No cells for this selection."))
  fill <- if (pal != "Default") disc_pal(levels(df$grp), pal)
          else if (has_tp) FG_PAL else setNames(c("#1565c0","#c62828"), levels(df$grp))
  lo <- min(df$y, na.rm = TRUE)
  n_lab <- as.data.frame(table(df$grp)); names(n_lab) <- c("grp", "n")
  ggplot(df, aes(grp, y, fill = grp)) +
    geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2) +
    geom_boxplot(width = .12, outlier.size = .3, alpha = .5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4,
                 fill = "white", colour = "grey20") +
    geom_text(data = n_lab, inherit.aes = FALSE, aes(x = grp, label = paste0("n=", n)),
              y = lo, vjust = 1.4, size = 3, colour = "#666") +
    scale_fill_manual(values = fill, na.value = "grey85") +
    guides(fill = "none") + theme_minimal(base_size = bs) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(x = NULL, y = labof(scol),
         title = paste0(labof(scol), " — ", ct,
                        if (stratum == "G1") " (G1 cells only)" else ""),
         caption = "White diamond = mean. Descriptive only: n = 1 animal per group.")
}
# proliferation vs cytokinesis scatter — the polyploidization quadrant
ploidy_scatter <- function(ct, bs = 13) {
  need_cols <- c("sig_prolif","sig_cytokinesis")
  validate(need(all(need_cols %in% names(meta)),
    "Run build_signature_scores.R for the proliferation/cytokinesis scores."))
  df <- meta; if (ct != "All" && "celltype" %in% names(df)) df <- df[as.character(df$celltype) == ct, ]
  df <- df[stats::complete.cases(df[need_cols]), ]
  validate(need(nrow(df), "No scored cells for this selection."))
  p <- ggplot(df, aes(sig_prolif, sig_cytokinesis, color = genotype)) +
    geom_point(size = .5, alpha = .35) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    theme_minimal(base_size = bs) +
    labs(x = "Proliferation score", y = "Cytokinesis score", color = NULL,
         title = paste0("Cycling vs cytokinesis — ", ct),
         caption = "High proliferation + low cytokinesis (lower-right) = karyokinesis without cytokinesis → polyploidization.")
  if ("timepoint" %in% names(df)) p <- p + facet_wrap(~ timepoint)
  p
}
# Transcriptional vs metabolic maturation, per CELL. 30k points plotted raw is an
# unreadable blob, so this draws density contours over a thinned point layer and puts
# each group's centroid on top; the WT->KO segment in each panel is the shift the KO
# question actually asks about, annotated with its size.
mat_scatter <- function(ct, bs = 13, stratum = "all", show_cells = TRUE) {
  need_cols <- c("sig_maturation","sig_metabolic")
  validate(need(all(need_cols %in% names(meta)),
    "Run build_signature_scores.R for the maturation/metabolic scores."))
  df <- meta; if (ct != "All" && has("celltype")) df <- df[as.character(df$celltype) == ct, ]
  df <- df[stats::complete.cases(df[need_cols]), ]
  if (stratum == "G1" && "Phase" %in% names(df)) df <- df[as.character(df$Phase) == "G1", , drop = FALSE]
  validate(need(nrow(df) > 10, "No cells with both scores for this selection."))
  has_tp <- "timepoint" %in% names(df)
  df$grp <- if (has_tp) factor(paste(df$genotype, df$timepoint, sep = "-"), levels = FG_GROUPS)
            else factor(df$genotype)
  df <- df[!is.na(df$grp), , drop = FALSE]
  pal <- if (has_tp) FG_PAL else setNames(c("#1565c0","#c62828"), levels(df$grp))
  # centroids, and the KO-WT displacement within each timepoint
  cen <- do.call(rbind, lapply(split(df, df$grp), function(x) if (!nrow(x)) NULL else data.frame(
    grp = x$grp[1], genotype = x$genotype[1],
    timepoint = if (has_tp) x$timepoint[1] else NA,
    x = mean(x$sig_maturation), y = mean(x$sig_metabolic), n = nrow(x), stringsAsFactors = FALSE)))
  seg <- NULL; shift_txt <- ""
  if (has_tp && !is.null(cen)) {
    parts <- lapply(split(cen, cen$timepoint), function(cc) {
      w <- cc[cc$genotype == "WT", ]; k <- cc[cc$genotype == "KO", ]
      if (!nrow(w) || !nrow(k)) return(NULL)
      data.frame(timepoint = cc$timepoint[1], x = w$x, y = w$y, xend = k$x, yend = k$y,
                 d = sqrt((k$x - w$x)^2 + (k$y - w$y)^2), dmat = k$x - w$x, stringsAsFactors = FALSE)
    })
    seg <- do.call(rbind, Filter(Negate(is.null), parts))
    if (!is.null(seg)) shift_txt <- paste0("KO−WT centroid shift — ",
      paste(sprintf("%s: %.3f (maturation %+.3f)", seg$timepoint, seg$d, seg$dmat), collapse = "; "))
  }
  # quadrant occupancy per group, relative to 0 on both scores
  qd <- do.call(rbind, lapply(split(df, df$grp), function(x) if (!nrow(x)) NULL else data.frame(
    grp = x$grp[1], pct = round(100 * mean(x$sig_maturation > 0 & x$sig_metabolic > 0), 1),
    stringsAsFactors = FALSE)))
  quad_txt <- if (is.null(qd)) "" else paste0("Mature+oxidative quadrant — ",
    paste(sprintf("%s %.1f%%", qd$grp, qd$pct), collapse = ", "))
  thin <- if (nrow(df) > 6000) df[sample.int(nrow(df), 6000), ] else df
  p <- ggplot(df, aes(sig_maturation, sig_metabolic, colour = grp))
  if (isTRUE(show_cells)) p <- p + geom_point(data = thin, size = .35, alpha = .16)
  p <- p +
    stat_density_2d(bins = 5, linewidth = .45, contour = TRUE) +
    geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70")
  if (!is.null(seg)) p <- p + geom_segment(data = seg, inherit.aes = FALSE,
    aes(x = x, y = y, xend = xend, yend = yend), colour = "grey25", linewidth = .6,
    arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed"))
  p <- p +
    geom_point(data = cen, inherit.aes = FALSE, aes(x, y, fill = grp),
               shape = 21, size = 4.2, stroke = 1.1, colour = "grey15") +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal, guide = "none") +
    theme_minimal(base_size = bs) +
    labs(x = "CM maturation score", y = "Metabolic maturation (FAO−glyc)", colour = NULL,
         title = paste0("Transcriptional vs metabolic maturation — ", ct,
                        if (stratum == "G1") " (G1 cells only)" else ""),
         subtitle = if (nzchar(shift_txt)) shift_txt else NULL,
         caption = paste0(quad_txt,
           "\nContours = per-group density; large circles = centroids; arrow = WT→KO.",
           "\nDescriptive only — n = 1 animal per group."))
  if (has_tp) p <- p + facet_wrap(~ timepoint)
  p
}
# ---- curated cell-cell communication (COMMUN$scores) --------------------------
commun_pathways <- function() if (is.null(COMMUN$scores)) character(0) else sort(unique(COMMUN$scores$pathway))
commun_heat_gg <- function(pathway, tp, metric) {
  validate(need(!is.null(COMMUN$scores),
    "Cell-cell communication isn't in this data build — run build_communication.R and redeploy."))
  d <- COMMUN$scores[COMMUN$scores$timepoint == tp, , drop = FALSE]
  if (!is.null(pathway) && pathway != "All") d <- d[d$pathway == pathway, , drop = FALSE]
  validate(need(nrow(d), "No interactions for this selection."))
  d$val <- d[[metric]]
  agg <- aggregate(val ~ sender + receiver, d, function(x) mean(x, na.rm = TRUE))
  fillscale <- if (metric == "delta") div_scale
               else scale_fill_gradient(low = "#eeeeee", high = "#800026", na.value = "grey92")
  ttl <- paste0(if (pathway == "All") "All pathways" else pathway, " — ", tp,
                " (", if (metric == "delta") "KO − WT" else metric, ")")
  p <- ggplot(agg, aes(receiver, sender, fill = val,
        text = paste0(sender, " → ", receiver, "<br>", metric, ": ", round(val, 3)))) +
    geom_tile(color = "grey92") + fillscale +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = "receiver", y = "sender", fill = metric, title = ttl)
  p
}
commun_heat <- function(pathway, tp, metric)
  ggplotly(commun_heat_gg(pathway, tp, metric), tooltip = "text") |> layout(margin = list(t = 40))
commun_table_df <- function(pathway, tp) {
  validate(need(!is.null(COMMUN$scores), "No communication table in this build."))
  d <- COMMUN$scores[COMMUN$scores$timepoint == tp, , drop = FALSE]
  if (pathway != "All") d <- d[d$pathway == pathway, , drop = FALSE]
  d <- d[order(-abs(d$delta)), ]
  d[, intersect(c("pathway","ligand","receptor","sender","receiver","WT","KO","delta"), names(d))]
}
commun_table <- function(pathway, tp) enr_dt(commun_table_df(pathway, tp))
# ---- reference-marker annotation check (REFMAP$confusion) --------------------
refmap_heat_gg <- function() {
  validate(need(!is.null(REFMAP$confusion),
    "Annotation check isn't in this data build — run build_refmap.R and redeploy."))
  d <- REFMAP$confusion
  p <- ggplot(d, aes(predicted, celltype, fill = prop,
        text = paste0(celltype, " → ", predicted, "<br>row prop: ", round(prop, 2), "<br>n: ", n))) +
    geom_tile(color = "grey92") +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey92") +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = "predicted (marker argmax)", y = "existing celltype", fill = "row prop",
         title = "Annotation concordance (row-normalised)")
  p
}
refmap_heat <- function() ggplotly(refmap_heat_gg(), tooltip = "text") |> layout(margin = list(t = 40))
refmap_table_df <- function() {
  validate(need(!is.null(REFMAP$confusion), "No annotation-check table in this build."))
  d <- REFMAP$confusion
  d[order(d$celltype, -d$prop), ]
}
refmap_table <- function() enr_dt(refmap_table_df())

# ---- four-group CM analysis (FG; build_fourgroup.R) ---------------------------
# WT-P0 / WT-P7 / KO-P0 / KO-P7 within each res-0.2 CM subcluster. Every contrast
# exists in two phase strata: "G1" (phase-matched, the default) and "all" (raw).
# Phase-matching matters because P7 was FACS cycling-enriched 4.5-5.2x and P0 was
# not, so a raw P0-vs-P7 contrast largely reads out the sort.
FG_MSG      <- "Four-group analysis isn't in this data build — run build_fourgroup.R and redeploy."
FG_GROUPS   <- if (!is.null(FG)) FG$built$groups else c("WT-P0","WT-P7","KO-P0","KO-P7")
FG_CLUSTERS <- if (!is.null(FG)) unique(FG$counts$cluster) else character(0)
FG_CTAB     <- if (!is.null(FG)) FG$built$contrasts else NULL
# light = P0, dark = P7; blue = WT, red = KO — same colour language as VOLC_PAL.
FG_PAL <- setNames(c("#90caf9","#1565c0","#ef9a9a","#c62828"),
                   c("WT-P0","WT-P7","KO-P0","KO-P7"))
FG_SORT_NOTE <- paste(
  "P7 was FACS cycling-enriched 4.5–5.2× and P0 essentially unenriched, so a raw",
  "P0-vs-P7 contrast reads out the sort as much as development. Contrasts default to",
  "the phase-matched (G1-only) stratum for this reason.")

fg_ok <- function() validate(need(!is.null(FG), FG_MSG))

# ---- four-group enrichment (build_fourgroup_enrichment.R) --------------------
# GO + GSEA for the four-group contrasts specifically. Distinct from ENR$sub,
# which is the same KO-vs-WT question POOLED across P0 and P7 -- that pooled
# version cannot answer "what changes in the KO at P7", which is what was asked.
FGE     <- app$enrich$fourgroup     # may be NULL on a bundle built before this
# ---- linked precomputed results (build_linked_results.R) ---------------------
# Results the upstream pipeline already computed, carried into the bundle rather than
# re-approximated in the app. Two of these supersede work the app does by hand: the
# CellChat run has permutation tests where the curated ligand score has none, and
# propeller tests abundance where the Composition tab only draws proportions.
LK      <- app$linked
LKM     <- app$linked_manifest
LK_MSG  <- paste("Linked upstream results aren't in this data build —",
                 "run build_linked_results.R and redeploy.")

# ---- CM object-mode diagnostic (build_cm_objectmode.R) -----------------------
# Review asked whether subsetting with comb[, celltype == "Cardiomyocyte"] leaves the
# embedding contaminated by the whole-heart object because the subset was never made
# standalone. cm_objectmode_check.R answers it by building the same cells three ways
# and diffing the results; this carries the answer, not an argument about it.
CMTEST     <- app$cmtest
CMTEST_MSG <- paste("The object-mode diagnostic isn't in this data build —",
                    "run our_analysis/04_integrate_annotate/cm_objectmode_check.R,",
                    "then build_cm_objectmode.R, and redeploy.")
cmtest_ok  <- function() validate(need(!is.null(CMTEST), CMTEST_MSG))

# ---- PC-dimension sweep (build_pcdims.R) ------------------------------------
# Both production embeddings carry dims 1:30 into FindNeighbors and RunUMAP with no
# recorded justification -- docs/01-upstream.qmd lists the PC count as a parameter the
# book cannot cite. This holds SCT/PCA/Harmony fixed and varies only that cut.
PCD     <- app$pcdims
PCD_MSG <- paste("The PC-dimension sweep isn't in this data build —",
                 "run our_analysis/04_integrate_annotate/pcdims_sweep.R for each object,",
                 "then build_pcdims.R, and redeploy.")
pcd_ok  <- function() validate(need(!is.null(PCD), PCD_MSG))

# ---- gene-set provenance (build_gene_provenance.R) --------------------------
# Every gene set this app scores or filters on, and where it came from. Most were typed
# into a script here with no citation; that is worth stating on the page rather than
# leaving a reader to assume the panels are pulled from somewhere.
GSP     <- app$genesets
GSP_MSG <- paste("The gene-set provenance audit isn't in this data build —",
                 "run our_analysis/05_analyses/gene_set_provenance.R,",
                 "then build_gene_provenance.R, and redeploy.")
gsp_ok  <- function() validate(need(!is.null(GSP), GSP_MSG))

# ---- clustering variants (build_clusterings.R) ------------------------------
# The CM labelling is a choice, not a fact: the PC-dimension sweep showed the subcluster
# count moving 7 -> 11 -> 12 with the PC cut. This makes that choice selectable and lets
# the downstream numbers follow it. Markers, composition and phase are recomputed live
# (presto is in the runtime image, ~1-3 s); pseudobulk DE and enrichment are precomputed
# offline per variant because DESeq2 and clusterProfiler are not.
# Immune cells that annotate.R labelled Cardiomyocyte (docs/05-cm-deepdive.qmd#cm12).
# Surfaced rather than filtered: measured, excluding them moves the CM cycling fraction
# 21.4% -> 21.3%, so exclusion machinery in every statistic would not be earning its
# complexity. What IS wrong is a per-subcluster row, and being able to see which cells
# they are is what makes that legible.
ICN     <- app$immune_contam
ICN_MSG <- paste("The immune-contamination flag isn't in this data build —",
                 "run our_analysis/05_analyses/cm_immune_contamination.R,",
                 "then build_immune_flag.R, and redeploy.")
CLU     <- app$clusterings
CLU_MSG <- paste("Clustering variants aren't in this data build —",
                 "run pcdims_sweep.R, cm_subcluster_analyze.R --variant=...,",
                 "then build_clusterings.R, and redeploy.")
clu_ok  <- function() validate(need(!is.null(CLU), CLU_MSG))
clu_choices <- function() {
  if (is.null(CLU)) return(character(0))
  ids <- names(CLU$variants)
  setNames(ids, vapply(ids, function(i) {
    v <- CLU$variants[[i]]
    sprintf("%s — %d subclusters%s", v$label, v$n_clusters,
            if (isTRUE(v$is_production)) "  (production)"
            else if (!isTRUE(v$has_downstream)) "  (labels only)" else "") }, ""))
}

FGE_MSG <- paste("Four-group enrichment isn't in this data build —",
                 "run build_fourgroup_enrichment.R and redeploy.")
fg_enr_ok <- function() validate(need(!is.null(FGE) && !is.null(FGE$go), FGE_MSG))
# One cluster x contrast x stratum slice. `kind` is "go" or "gsea".
fg_enr_df <- function(kind, cluster, contrast, stratum, ont = "BP", direction = NULL) {
  fg_enr_ok()
  d <- FGE[[kind]]
  validate(need(!is.null(d), FGE_MSG))
  d <- d[d$cluster == cluster & d$contrast == contrast & d$stratum == stratum, , drop = FALSE]
  if (kind == "go" && "ontology" %in% names(d)) d <- d[d$ontology == ont, , drop = FALSE]
  if (!is.null(direction) && "direction" %in% names(d)) d <- d[d$direction == direction, , drop = FALSE]
  d
}
# The audit row behind an empty result: an empty GO table means either "tested,
# nothing passed" or "list too small to test", and those are different answers.
fg_enr_why <- function(cluster, contrast, stratum, ont, direction) {
  if (is.null(FGE) || is.null(FGE$audit)) return("")
  a <- FGE$audit
  a <- a[a$cluster == cluster & a$contrast == contrast & a$stratum == stratum &
         a$ontology == ont & a$direction == direction, , drop = FALSE]
  if (!nrow(a)) return("")
  sprintf("%d genes in, %d-gene universe, selected by [%s].",
          a$n_input[1], a$n_universe[1], a$input_rule[1])
}
# Two DE grids over the same contrasts, trading gene coverage against cell coverage.
# Neither dominates: the broad matrix has ~24k genes but 8k cells, so CM2's KO-P0 arm
# falls below the floor and CM4/CM9 lose their G1 strata; the curated matrix keeps all
# 30k cells so every contrast runs, over ~2.2k genes. The user picks per question.
fg_grid_choices <- function() {
  if (is.null(FG)) return(c("Gene coverage" = "de"))
  g <- c(setNames("de", sprintf("Gene coverage — %s genes, %s CM cells",
                   format(FG$built$n_genes, big.mark = ","),
                   format(FG$built$n_cells_de, big.mark = ","))))
  if (!is.null(FG$de2)) g <- c(g, setNames("de2", sprintf("Cell coverage — %s genes, all %s CM cells",
                   format(FG$built$n_genes2 %||% 0, big.mark = ","),
                   format(FG$built$n_cells_total %||% 0, big.mark = ","))))
  g
}
fg_grid <- function(grid) if (identical(grid, "de2") && !is.null(FG$de2)) FG$de2 else FG$de
fg_skipset <- function(grid) if (identical(grid, "de2")) FG$skipped2 else FG$skipped
# is this contrast available in the OTHER grid? worth saying so rather than just "no"
fg_other_has <- function(cluster, contrast, stratum, grid) {
  other <- if (identical(grid, "de2")) FG$de else FG$de2
  if (is.null(other)) return(FALSE)
  !is.null(other[[cluster]][[paste0(contrast, "__", stratum)]])
}
# cluster dropdown: "All cardiomyocytes" + the usual "CM2 · Ventricular (2397 cells)"
fg_cluster_choices <- function() {
  if (is.null(FG)) return(character(0))
  setNames(FG_CLUSTERS, vapply(FG_CLUSTERS, function(x)
    if (x == "AllCM") "All cardiomyocytes" else sub_label("0.2", x), ""))
}
fg_contrast_choices <- function() {
  if (is.null(FG_CTAB)) return(character(0))
  setNames(FG_CTAB$key, FG_CTAB$label)
}
fg_ct <- function(key) {
  if (is.null(FG_CTAB)) return(NULL)
  r <- FG_CTAB[FG_CTAB$key == key, , drop = FALSE]
  if (nrow(r)) as.list(r[1, ]) else NULL
}
# when a table is missing, say WHY (and with what cell counts) rather than "no data"
fg_skip_msg <- function(cluster, contrast, stratum, grid = "de") {
  hint <- if (fg_other_has(cluster, contrast, stratum, grid))
    paste0(" This contrast IS available in the “",
           if (identical(grid, "de2")) "Gene coverage" else "Cell coverage",
           "” matrix — switch it in the sidebar.") else ""
  s <- fg_skipset(grid)
  if (!is.null(s)) {
    r <- s[s$cluster == cluster & s$contrast == contrast & s$stratum == stratum, , drop = FALSE]
    if (nrow(r)) return(sprintf(
      "Not computed — %s. Group A: %d cells, group B: %d cells (floor is %d).%s%s",
      r$reason[1], r$n_A[1], r$n_B[1], FG$built$min_cells,
      if (stratum == "G1") " Try the “All cells” stratum." else "", hint))
  }
  paste0("No DE table for this selection.", hint)
}
fg_de <- function(cluster, contrast, stratum, grid = "de") {
  fg_ok(); req(cluster, contrast, stratum)
  d <- fg_grid(grid)[[cluster]][[paste0(contrast, "__", stratum)]]
  validate(need(!is.null(d) && nrow(d), fg_skip_msg(cluster, contrast, stratum, grid)))
  d
}
# The three caveats a four-group contrast has to carry wherever it is shown: the cells
# per arm, a hard warning when one arm is too thin to trust, and the sort-confound flag
# on a raw P0-vs-P7 comparison. Shared by the Four-group tab and the CM deep-dive so the
# two cannot drift apart -- a caveat that appears on one tab and not the other is worse
# than no caveat, because it makes the quiet tab look safe.
de_context_note <- function(cluster, ct, stratum, grid, d) {
  ok <- is.data.frame(d) && nrow(d) > 0 && all(c("n_A", "n_B") %in% names(d))
  nA <- if (ok) d$n_A[1] else NA_integer_
  nB <- if (ok) d$n_B[1] else NA_integer_
  ok <- ok && !is.na(nA) && !is.na(nB)
  thin <- FG$built$thin_cells %||% 50
  div(style = "font-size:13px;margin-bottom:6px",
      HTML(paste0("<b>", ct$label, "</b> in ",
                  if (identical(cluster, "AllCM")) "all cardiomyocytes" else cluster,
                  if (identical(stratum, "G1")) ", G1 cells only" else ", all cells",
                  if (ok) sprintf(" &mdash; %d vs %d cells, %d genes past the gate", nA, nB, nrow(d)) else "",
                  if (identical(grid, "de2")) " &middot; curated panel" else "")),
      if (ok && min(nA, nB) < thin)
        div(style = "color:#c62828;font-weight:600",
            HTML(sprintf("&#9888; Smallest group is %d cells — treat this contrast as unreliable%s.",
              min(nA, nB),
              if (max(nA, nB) / max(min(nA, nB), 1) >= 10)
                sprintf(" (%.0f× imbalance)", max(nA, nB) / max(min(nA, nB), 1)) else ""))),
      if (identical(stratum, "all") && grepl("P0_vs_P7", ct$key))
        span(style = "color:#c62828",
             HTML(" &middot; <b>sort-confounded</b> — P7 was FACS cycling-enriched 4.5–5.2× relative to P0; read the G1 stratum instead.")))
}
# One enrichment question, two precomputed sources. ENR$sub is KO-vs-WT POOLED over P0
# and P7 (build_subcluster_enrichment.R); FGE is one timepoint-specific contrast at a
# time (build_fourgroup_enrichment.R). They render identically, so the deep-dive picks
# between them on the same Comparison dropdown its DE tab uses. Only the direction keys
# differ -- KO_up/KO_down vs A_up/B_up -- so callers ask for "up"/"down" and this maps.
enr_src_df <- function(kind, cluster, ct, stratum = "all", ont = "BP", dir = NULL) {
  if (is.null(ct)) {
    d <- enr_sub_df(kind, cluster)
    if (!is.null(dir) && "direction" %in% names(d))
      d <- d[d$direction == if (identical(dir, "up")) "KO_up" else "KO_down", , drop = FALSE]
    return(d)
  }
  fg_enr_df(kind, cluster, ct$key, stratum, ont,
            if (is.null(dir)) NULL else if (identical(dir, "up")) "A_up" else "B_up")
}
# direction wording + a title fragment for whichever source is selected
enr_src_labs <- function(ct) {
  if (is.null(ct)) list(up = "up in KO", dn = "up in WT", label = "KO vs WT (P0 + P7 pooled)")
  else list(up = ct$pos, dn = ct$neg, label = ct$label)
}
# one row per cluster: counts + percentages for the four groups, with the
# under-powered arms named explicitly rather than left to be discovered.
fg_counts_wide <- function() {
  fg_ok(); d <- FG$counts; cl <- FG_CLUSTERS
  st <- subType[[paste0("res", FG$built$res)]]
  out <- data.frame(cluster = cl, stringsAsFactors = FALSE)
  if (!is.null(st)) out$subtype <- st$nearest_CM_subtype[match(cl, st$subcluster)]
  key <- paste(d$cluster, d$group)
  for (g in FG_GROUPS) {
    i <- match(paste(cl, g), key)
    out[[g]] <- d$n[i]
    out[[paste0(g, " %")]] <- d$pct_of_cluster[i]
    if ("n_G1" %in% names(d)) out[[paste0(g, " G1")]] <- d$n_G1[i]
  }
  out$total <- rowSums(out[, FG_GROUPS, drop = FALSE], na.rm = TRUE)
  # name the arms that can't carry a contrast, and say which way they fail —
  # "too few cells" kills the contrast, "thin in G1" only kills the default stratum.
  out$underpowered <- vapply(cl, function(c) {
    r <- d[d$cluster == c & d$status != "ok", , drop = FALSE]
    if (!nrow(r)) return("—")
    paste(sprintf("%s (%s)", r$group, sub(" \\(.*", "", r$status)), collapse = ", ") }, "")
  out
}
fg_counts_plot <- function(mode = "prop", bs = 13) {
  fg_ok(); d <- FG$counts
  d$group   <- factor(d$group, levels = FG_GROUPS)
  d$cluster <- factor(d$cluster, levels = FG_CLUSTERS)
  d$y <- if (mode == "prop") d$pct_of_cluster else d$n
  ggplot(d, aes(cluster, y, fill = group)) +
    geom_col(position = if (mode == "prop") "stack" else position_dodge(width = .85),
             width = .8) +
    scale_fill_manual(values = FG_PAL, na.value = "grey85") +
    theme_minimal(base_size = bs) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = "CM subcluster", fill = NULL,
         y = if (mode == "prop") "% of subcluster" else "cell count",
         title = paste0("Four-group composition per CM subcluster — res ", FG$built$res))
}
# cell-cycle phase composition per cluster x group. Answers the G1-proportion
# question and makes the P0->P7 S-phase jump (the sort) visible at the same time.
fg_phase_plot <- function(clusters, bs = 13) {
  fg_ok(); validate(need(!is.null(FG$phase), "No phase table in this data build."))
  d <- FG$phase[FG$phase$cluster %in% clusters, , drop = FALSE]
  validate(need(nrow(d), "Pick at least one subcluster."))
  d$group   <- factor(d$group, levels = FG_GROUPS)
  d$Phase   <- factor(d$Phase, levels = c("G1","S","G2M"))
  d$cluster <- factor(d$cluster, levels = intersect(FG_CLUSTERS, clusters))
  lab <- d[d$Phase == "G1", ]
  ggplot(d, aes(group, pct, fill = Phase)) +
    geom_col(width = .8) + facet_wrap(~ cluster) +
    geom_text(data = lab, inherit.aes = FALSE, y = 4, size = 3,
              colour = "white", fontface = "bold",
              aes(x = group, label = ifelse(is.na(pct), "", paste0(round(pct), "%")))) +
    scale_fill_manual(values = c(G1 = "#bdbdbd", S = "#1565c0", G2M = "#c62828")) +
    theme_minimal(base_size = bs) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = "% of cells in group", fill = NULL,
         title = "Cell-cycle phase composition — four groups",
         caption = "G1 % is printed on each bar. The P0→P7 S-phase jump is largely the FACS enrichment.")
}
# per-cell score violins, four groups, faceted by cluster. Built from cmm (real
# per-cell values) rather than FG$scores (summaries) so the distribution shows.
fg_score_plot <- function(scol, clusters, stratum, bs = 13) {
  validate(need(scol %in% names(cmm),
    "Score not in this data build — run build_signature_scores.R and redeploy."))
  df <- cmm
  df$cluster <- if ("cm_subcluster" %in% names(df)) df$cm_subcluster else
                paste0("CM", df[[cm_subcol("0.2")]])
  parts <- list()
  if ("AllCM" %in% clusters) { a <- df; a$cluster <- "AllCM"; parts[[1]] <- a }
  rest <- setdiff(clusters, "AllCM")
  if (length(rest)) parts[[length(parts) + 1]] <- df[df$cluster %in% rest, , drop = FALSE]
  validate(need(length(parts), "Pick at least one subcluster."))
  df <- do.call(rbind, parts)
  if (stratum == "G1") df <- df[as.character(df$Phase) == "G1", , drop = FALSE]
  df$group   <- factor(paste(df$genotype, df$timepoint, sep = "-"), levels = FG_GROUPS)
  df$cluster <- factor(df$cluster, levels = intersect(FG_CLUSTERS, clusters))
  df <- df[!is.na(df[[scol]]) & !is.na(df$group), , drop = FALSE]
  validate(need(nrow(df) > 0, paste0(
    "No scored cells for this selection",
    if (stratum == "G1") " — this subcluster may have no G1 cells." else ".")))
  df$y <- df[[scol]]
  ggplot(df, aes(group, y, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2) +
    geom_boxplot(width = .12, outlier.size = .3, alpha = .5) +
    facet_wrap(~ cluster) + scale_fill_manual(values = FG_PAL) + guides(fill = "none") +
    theme_minimal(base_size = bs) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = labof(scol),
         title = paste0(labof(scol), " — four groups",
                        if (stratum == "G1") " (G1 cells only)" else ""))
}
# ---- does the KO leave P7 cardiomyocytes less mature, and more cycling? -------
# The email's actual question. Everything needed was already stored per group; what
# was missing was the subtraction. Per subcluster this contrasts KO against WT at each
# timepoint on both readouts, so "the P7 gap is bigger than the P0 gap" is a number
# rather than an impression formed from two violins.
fg_summary_df <- function(score = "sig_maturation_nocc", stratum = "G1") {
  fg_ok()
  validate(need(!is.null(FG$scores), "Score summaries aren't in this build — run build_fourgroup.R."))
  sc <- FG$scores[FG$scores$score == score & FG$scores$stratum == stratum, , drop = FALSE]
  validate(need(nrow(sc), "No score summary for this selection."))
  ph <- if (!is.null(FG$phase)) FG$phase[FG$phase$Phase == "G1", , drop = FALSE] else NULL
  pick <- function(df, cl, gr, col) {
    if (is.null(df)) return(NA_real_)
    v <- df[[col]][df$cluster == cl & df$group == gr]
    if (length(v)) v[1] else NA_real_
  }
  rows <- lapply(FG_CLUSTERS, function(cl) {
    mn <- function(g) pick(sc, cl, g, "mean")
    sd <- function(g) pick(sc, cl, g, "sd")
    nn <- function(g) pick(sc, cl, g, "n")
    # KO - WT, plus Cohen's d so a shift is readable against the spread it sits in
    gap <- function(ko, wt) {
      if (is.na(mn(ko)) || is.na(mn(wt))) return(c(NA_real_, NA_real_))
      diff <- mn(ko) - mn(wt)
      sp <- suppressWarnings(sqrt(((nn(ko) - 1) * sd(ko)^2 + (nn(wt) - 1) * sd(wt)^2) /
                                  max(nn(ko) + nn(wt) - 2, 1)))
      c(diff, if (is.finite(sp) && sp > 0) diff / sp else NA_real_)
    }
    p0 <- gap("KO-P0", "WT-P0"); p7 <- gap("KO-P7", "WT-P7")
    g1 <- function(g) pick(ph, cl, g, "pct")
    # how much of WT's P0->P7 maturation gain does the KO achieve?
    wt_gain <- mn("WT-P7") - mn("WT-P0"); ko_gain <- mn("KO-P7") - mn("KO-P0")
    data.frame(
      cluster = cl,
      G1_pct_WT_P7 = g1("WT-P7"), G1_pct_KO_P7 = g1("KO-P7"),
      G1_gap_P7 = round(g1("KO-P7") - g1("WT-P7"), 1),
      G1_gap_P0 = round(g1("KO-P0") - g1("WT-P0"), 1),
      mat_WT_P7 = round(mn("WT-P7"), 3), mat_KO_P7 = round(mn("KO-P7"), 3),
      mat_gap_P7 = round(p7[1], 3), d_P7 = round(p7[2], 2),
      mat_gap_P0 = round(p0[1], 3), d_P0 = round(p0[2], 2),
      KO_gain_pct_of_WT = if (is.na(wt_gain) || is.na(ko_gain) || wt_gain == 0) NA_real_
                          else round(100 * ko_gain / wt_gain, 1),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  # the verdict column: the hypothesis is less mature AND less G1 at P7, and more so
  # at P7 than at P0. Anything else is spelled out rather than left to interpretation.
  d$verdict <- with(d, ifelse(
    is.na(mat_gap_P7) | is.na(G1_gap_P7), "insufficient data",
    ifelse(mat_gap_P7 < 0 & G1_gap_P7 < 0 &
             (is.na(mat_gap_P0) | abs(mat_gap_P7) > abs(mat_gap_P0)),
           "less mature + more cycling at P7",
    ifelse(mat_gap_P7 < 0 & G1_gap_P7 < 0, "less mature + more cycling, but not P7-specific",
    ifelse(mat_gap_P7 < 0, "less mature, but not more cycling",
    ifelse(G1_gap_P7 < 0, "more cycling, but not less mature", "neither"))))))
  d
}
# the same thing as a picture: where each group's score sits, with the KO-WT gap drawn
fg_summary_plot <- function(score = "sig_maturation_nocc", stratum = "G1",
                            clusters = NULL, bs = 13) {
  fg_ok()
  validate(need(!is.null(FG$scores), "Score summaries aren't in this build."))
  sc <- FG$scores[FG$scores$score == score & FG$scores$stratum == stratum, , drop = FALSE]
  if (!is.null(clusters) && length(clusters)) sc <- sc[sc$cluster %in% clusters, , drop = FALSE]
  validate(need(nrow(sc), "Pick at least one subcluster."))
  sc$group <- factor(sc$group, levels = FG_GROUPS)
  sc$cluster <- factor(sc$cluster, levels = intersect(FG_CLUSTERS, unique(sc$cluster)))
  sc$timepoint <- sub("^.*-", "", as.character(sc$group))
  sc$genotype  <- sub("-.*$", "", as.character(sc$group))
  # segment joining WT to KO within each timepoint = the gap the question is about
  seg <- do.call(rbind, lapply(split(sc, list(sc$cluster, sc$timepoint), drop = TRUE), function(x) {
    w <- x[x$genotype == "WT", ]; k <- x[x$genotype == "KO", ]
    if (!nrow(w) || !nrow(k)) return(NULL)
    data.frame(cluster = x$cluster[1], timepoint = x$timepoint[1],
               y = w$mean, yend = k$mean, gap = k$mean - w$mean, stringsAsFactors = FALSE)
  }))
  p <- ggplot(sc, aes(timepoint, mean, colour = group)) +
    { if (!is.null(seg)) geom_segment(data = seg, inherit.aes = FALSE,
        aes(x = timepoint, xend = timepoint, y = y, yend = yend),
        colour = "grey45", linewidth = .8,
        arrow = grid::arrow(length = grid::unit(6, "pt"), type = "closed")) } +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = .12, linewidth = .5) +
    geom_point(size = 3) +
    { if (!is.null(seg)) geom_text(data = seg, inherit.aes = FALSE,
        aes(x = timepoint, y = pmin(y, yend), label = sprintf("%+.2f", gap)),
        vjust = 1.9, size = 3, colour = "grey25") } +
    scale_colour_manual(values = FG_PAL) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_minimal(base_size = bs) +
    labs(x = NULL, y = paste0(labof(score), " (mean ± 95% CI)"), colour = NULL,
         title = paste0(labof(score), " — KO vs WT at each age",
                        if (stratum == "G1") " (G1 cells only)" else ""),
         caption = paste("Arrow runs WT → KO; the number is the gap.",
                         "\nG1-only holds cycling composition fixed, so the gap is not the FACS sort.",
                         "\nDescriptive only — n = 1 animal per group."))
  p
}

# ---- maturation axis x P7 KO-vs-WT -------------------------------------------
FG_QUAD <- c(immature_up_in_KO = "#c62828", mature_down_in_KO = "#1565c0",
             immature_down_in_KO = "#cccccc", mature_up_in_KO = "#cccccc", ns = "#e8e8e8")
fg_int_msg <- paste("Needs the maturation scores — run build_signature_scores.R,",
                    "then build_fourgroup.R, then redeploy.")
fg_intersect_df <- function(cluster, quadrants = NULL, hide_conf = TRUE) {
  fg_ok(); validate(need(!is.null(FG$intersect), fg_int_msg))
  d <- FG$intersect[FG$intersect$cluster == cluster, , drop = FALSE]
  validate(need(nrow(d), "No intersection rows for this subcluster."))
  if (isTRUE(hide_conf)) d <- d[!d$confounder, , drop = FALSE]
  if (!is.null(quadrants) && length(quadrants)) d <- d[d$quadrant %in% quadrants, , drop = FALSE]
  d[order(-abs(d$p7ko_log2FC)), , drop = FALSE]
}
fg_quadrant_plot <- function(cluster, hide_conf = TRUE, label_n = 20, bs = 13) {
  d <- fg_intersect_df(cluster, NULL, hide_conf)
  d$quadrant <- factor(d$quadrant, levels = names(FG_QUAD))
  hit <- d[d$quadrant %in% c("immature_up_in_KO","mature_down_in_KO"), , drop = FALSE]
  lab <- hit[order(-abs(hit$mat_auc - 0.5) - abs(hit$p7ko_log2FC)), , drop = FALSE]
  lab <- head(lab, label_n)
  # x is AUC, not log2FC: immature markers are far higher dynamic-range genes than
  # mature ones here (log2FC spans -1.86 to +0.43), so a log2FC axis would put every
  # classified gene on the left. AUC is rank-based, so the axis is symmetric and the
  # point position agrees with the colour.
  mid <- 0.5
  # an all-NA or empty frame makes max() return -Inf, which silently yields a
  # degenerate plot rather than an error — fall back to a sane range instead
  span <- function(v, default) { m <- suppressWarnings(max(abs(v), na.rm = TRUE))
                                 if (!is.finite(m) || m <= 0) default else m }
  xr <- span(d$mat_auc - mid, 0.1); yr <- span(d$p7ko_log2FC, 1)
  ggplot(d, aes(mat_auc, p7ko_log2FC)) +
    annotate("rect", xmin = mid - xr, xmax = mid, ymin = 0, ymax = yr, fill = "#c62828", alpha = .06) +
    annotate("rect", xmin = mid, xmax = mid + xr, ymin = -yr, ymax = 0, fill = "#1565c0", alpha = .06) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = mid, colour = "grey70") +
    geom_point(aes(colour = quadrant), size = 1.2, alpha = .7) +
    geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.7, check_overlap = TRUE) +
    scale_colour_manual(values = FG_QUAD, drop = FALSE) +
    theme_minimal(base_size = bs) +
    labs(x = "maturation association, AUC  (← immature | 0.5 | mature →)",
         y = "P7 KO vs WT  log2FC  (↑ up in KO)", colour = NULL,
         title = paste0("Maturation axis × P7 KO response — ",
                        if (cluster == "AllCM") "all cardiomyocytes" else cluster),
         caption = paste("Shaded: immature genes up in P7 KO (red) and mature genes down in P7 KO (blue)",
                         "— the two quadrants consistent with delayed maturation.",
                         "\nMaturation axis uses the cycle-free score, so it is not circular with cycling."))
}
# ---- candidate genes: computed live, so any gene can be asked about ----------
FG_SHORTLIST <- intersect(
  c("Birc5","Foxm1","Rrm2","Aurkb","Prc1","Gabbr2","Tcf4","Adamts9"), ALL_GENES)
# mean expression + % expressing per cluster x four-group, for a set of genes
fg_candidate_df <- function(genes, clusters, stratum) {
  validate(need(length(genes), "Pick at least one gene."))
  df <- cmm
  df$cluster <- if ("cm_subcluster" %in% names(df)) df$cm_subcluster else
                paste0("CM", df[[cm_subcol("0.2")]])
  parts <- list()
  if ("AllCM" %in% clusters) { a <- df; a$cluster <- "AllCM"; parts[[1]] <- a }
  rest <- setdiff(clusters, "AllCM")
  if (length(rest)) parts[[length(parts) + 1]] <- df[df$cluster %in% rest, , drop = FALSE]
  validate(need(length(parts), "Pick at least one subcluster."))
  df <- do.call(rbind, parts)
  if (stratum == "G1") df <- df[as.character(df$Phase) == "G1", , drop = FALSE]
  df$group <- factor(paste(df$genotype, df$timepoint, sep = "-"), levels = FG_GROUPS)
  df <- df[!is.na(df$group), , drop = FALSE]
  validate(need(nrow(df), "No cells for this selection."))
  out <- list()
  for (g in genes) {
    v <- expr_vec(g, df$cell)
    if (all(is.na(v))) next
    for (cl in unique(df$cluster)) for (gr in FG_GROUPS) {
      s <- df$cluster == cl & df$group == gr & !is.na(v)
      if (!any(s)) next
      out[[length(out) + 1]] <- data.frame(
        gene = g, cluster = cl, group = gr, n = sum(s),
        pct_expressing = round(100 * mean(v[s] > 0), 1),
        mean_expr = round(mean(v[s]), 3), stringsAsFactors = FALSE)
    }
  }
  validate(need(length(out), "None of the selected genes are in this data build."))
  d <- do.call(rbind, out)
  d$cluster <- factor(d$cluster, levels = intersect(FG_CLUSTERS, unique(d$cluster)))
  d$group   <- factor(d$group, levels = FG_GROUPS)
  d$gene    <- factor(d$gene, levels = intersect(genes, unique(d$gene)))
  d
}
fg_candidate_plot <- function(genes, clusters, stratum, bs = 13) {
  d <- fg_candidate_df(genes, clusters, stratum)
  ggplot(d, aes(group, gene, size = pct_expressing, colour = mean_expr)) +
    geom_point() + facet_wrap(~ cluster) +
    scale_size_area(max_size = 9) +
    scale_color_viridis_c(option = "magma", direction = -1) +
    theme_minimal(base_size = bs) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = NULL, size = "% expressing", colour = "mean expr",
         title = paste0("Candidate genes across CM subclusters × four groups",
                        if (stratum == "G1") " (G1 cells only)" else ""))
}
# Is a gene's KO effect P7-specific, and is it concentrated in CM2/CM4/CM5?
# Read straight off the precomputed contrasts, so it matches the DEG tables.
FG_PRIORITY <- c("CM2","CM4","CM5")
fg_specificity_df <- function(genes, stratum, grid = "de") {
  fg_ok(); validate(need(length(genes), "Pick at least one gene."))
  pull <- function(cl, key, g) {
    d <- fg_grid(grid)[[cl]][[paste0(key, "__", stratum)]]
    if (is.null(d)) return(NA_real_)
    d$log2FoldChange[match(g, d$gene)]
  }
  rows <- list()
  for (g in genes) for (cl in FG_CLUSTERS) {
    p7 <- pull(cl, "P7_KO_vs_WT", g); p0 <- pull(cl, "P0_KO_vs_WT", g)
    if (is.na(p7) && is.na(p0)) next
    rows[[length(rows) + 1]] <- data.frame(
      gene = g, cluster = cl,
      P7_KO_vs_WT = p7, P0_KO_vs_WT = p0,
      P7_specificity = round(abs(p7) - abs(p0), 3), stringsAsFactors = FALSE)
  }
  validate(need(length(rows), paste0(
    "No contrasts cover these genes in the “",
    if (stratum == "G1") "G1" else "All cells", "” stratum.")))
  d <- do.call(rbind, rows)
  # per gene: is the effect bigger inside CM2/CM4/CM5 than outside?
  agg <- do.call(rbind, lapply(split(d, d$gene), function(x) {
    pri <- mean(abs(x$P7_KO_vs_WT[x$cluster %in% FG_PRIORITY]), na.rm = TRUE)
    oth <- mean(abs(x$P7_KO_vs_WT[!x$cluster %in% c(FG_PRIORITY, "AllCM")]), na.rm = TRUE)
    data.frame(gene = x$gene[1],
               CM2_4_5_mean_absLFC = round(pri, 3),
               other_clusters_mean_absLFC = round(oth, 3),
               priority_concentration = round(pri - oth, 3), stringsAsFactors = FALSE)
  }))
  d <- merge(d, agg, by = "gene", all.x = TRUE)
  # the email's third question: is the gene associated with a less mature or more
  # cycling state? Read off the same axes the Gene map and intersection tabs use, so
  # a candidate's state and its KO behaviour sit on one row.
  if (!is.null(GM)) {
    i <- match(d$gene, GM$gene)
    d$mat_auc <- GM$mat_auc[i]
    d$state_maturation <- GM$mat_class[i]
    if ("cyc_auc" %in% names(GM)) {
      d$cyc_auc <- GM$cyc_auc[i]
      d$state_cycling <- GM$cyc_class[i]
    }
  } else if (!is.null(FG$maturation)) {
    i <- match(d$gene, FG$maturation$gene)
    d$mat_auc <- FG$maturation$mat_auc[i]; d$state_maturation <- FG$maturation$mat_class[i]
  }
  d[order(d$gene, factor(d$cluster, levels = FG_CLUSTERS)), , drop = FALSE]
}
# --- "any additional top candidates identified from the analyses above" ----------
# The email asks for the shortlist PLUS whatever the other analyses threw up, so the
# gene box can be populated from each of them rather than retyped by hand.
FG_CAND_SOURCES <- c("Email shortlist" = "shortlist",
                     "Top four-group DE — P7 KO vs WT" = "de",
                     "Intersection hits — immature-up / mature-down" = "intersect",
                     "Gene map extremes — by distance" = "genemap")
fg_candidate_pool <- function(src, n = 20, grid = "de") {
  drop_conf_genes <- function(g) setdiff(g, CONF)
  if (is.null(src) || src == "shortlist") return(FG_SHORTLIST)
  if (src == "de") {
    cls <- intersect(c("AllCM", FG_PRIORITY), FG_CLUSTERS)
    g <- unlist(lapply(cls, function(cl) {
      d <- fg_grid(grid)[[cl]][["P7_KO_vs_WT__G1"]]
      if (is.null(d)) d <- fg_grid(grid)[[cl]][["P7_KO_vs_WT__all"]]
      if (is.null(d)) return(character(0))
      d <- d[!d$confounder, , drop = FALSE]
      head(d$gene[order(-abs(d$log2FoldChange))], n)
    }))
    return(drop_conf_genes(unique(g)))
  }
  if (src == "intersect" && !is.null(FG$intersect)) {
    d <- FG$intersect[!FG$intersect$confounder &
          FG$intersect$cluster %in% intersect(c("AllCM", FG_PRIORITY), FG_CLUSTERS) &
          FG$intersect$quadrant %in% c("immature_up_in_KO", "mature_down_in_KO"), , drop = FALSE]
    if (!nrow(d)) return(character(0))
    return(drop_conf_genes(unique(d$gene[order(-abs(d$p7ko_log2FC))])))
  }
  if (src == "genemap" && !is.null(GM)) {
    d <- GM[is.na(GM$in_score_set), , drop = FALSE]   # exclude the circular ones
    return(drop_conf_genes(head(d$gene[order(-d$distance)], n * 2)))
  }
  FG_SHORTLIST
}

# ---- the curated gene lists behind every sig_* score --------------------------
# They live as a literal in build_signature_scores.R and are NOT written into
# app_data.rds, so rather than keep a second copy here -- which would drift the first
# time one is edited -- the app parses that one expression out of the source. Same trick
# analysis/2026-08-21_email/01_de.R already uses, for the same reason. The builder ships
# next to app.R in every deploy path (rsconnect, both Dockerfiles); if it is ever missing,
# the score panels say so rather than showing a stale copy.
SCORE_SETS <- local({
  path <- Filter(file.exists, c("build_signature_scores.R",
                                file.path("shiny_app", "build_signature_scores.R")))
  s <- if (!length(path)) NULL else tryCatch({
    got <- NULL
    for (x in parse(path[1]))
      if (is.null(got) && is.call(x) && identical(as.character(x[[1]]), "<-") &&
          identical(as.character(x[[2]]), "SETS")) got <- eval(x[[3]])
    got
  }, error = function(e) NULL)
  if (is.null(s)) NULL else lapply(s, function(g) setdiff(g, CONF))
})
# What each score is FOR -- the question it exists to answer. The maths is shared; the
# purpose is not, and it is the part that is impossible to reconstruct from the code.
SCORE_PURPOSE <- c(
  sig_maturation = paste(
    "How far a cardiomyocyte has moved along the fetal-to-adult transition. Positive means the",
    "mature program dominates. This is the axis the P0-vs-P7 comparison is about, and the one",
    "the headline KO question asks about: do P7 knockout cardiomyocytes sit at a lower",
    "maturation score than P7 wild-type?"),
  sig_maturation_nocc = paste(
    "The same maturation axis with the three cell-cycle genes dropped from the immature pole.",
    "Use this one whenever the argument involves cycling. sig_maturation's immature program",
    "contains Mki67, Top2a and Ccnd1, so arguing 'less mature, therefore more cycling-competent'",
    "from it would be partly circular. The gene map and the intersection tab default to this."),
  sig_metabolic = paste(
    "Which fuel the cell is set up to burn. Positive means fatty-acid oxidation and oxidative",
    "phosphorylation; negative means glycolysis. Cardiomyocytes switch from glycolysis to fatty",
    "acids over the first postnatal week, so this is the metabolic half of maturation and should",
    "track sig_maturation if both are measuring the same transition. Note the OXPHOS genes",
    "(Ppargc1a, Cox6a2, Ndufa4, Sdha, Etfa) sit inside the FAO set, so this score cannot",
    "separate fatty-acid oxidation from oxidative phosphorylation."),
  sig_mat_mature   = "The mature pole of the maturation axis on its own, before the subtraction.",
  sig_mat_immature = "The immature pole on its own, cell-cycle genes included.",
  sig_mat_immature_nocc = "The immature pole with the three cell-cycle genes removed.",
  sig_faox         = "The fatty-acid-oxidation pole of the metabolic axis on its own.",
  sig_glycolysis   = "The glycolytic pole of the metabolic axis on its own.",
  sig_prolif       = "The G2/M proliferation program: is this cell cycling at all?",
  sig_cytokinesis  = "The cytokinesis machinery specifically: can this cell finish a division?",
  sig_ploidy = paste(
    "Proliferation minus cytokinesis. Cycling WITHOUT the machinery to complete division is the",
    "route to binucleation and endoreduplication, so a high score flags cells likely to become",
    "polyploid. A heuristic built from two module scores, not a ploidy measurement."),
  sig_ccexit = "The cell-cycle exit / arrest program (CDK inhibitors, Rb1, Meis1, Btg2, Gadd45a).")

# What a sig_* score is for, how it is built, and the actual genes -- the gene lists being
# the one thing app_data.rds does not carry and the thing people actually ask for.
score_def_ui <- function(scol) {
  r <- if (!is.null(SCOREMETA)) SCOREMETA[SCOREMETA$score == scol, , drop = FALSE] else NULL
  purpose <- if (scol %in% names(SCORE_PURPOSE)) SCORE_PURPOSE[[scol]] else NULL
  parts <- if (!is.null(r) && nrow(r))
    trimws(strsplit(r$sets[1], " - ", fixed = TRUE)[[1]]) else character(0)
  sgn <- if (length(parts) == 2) c("+", "−") else rep("+", length(parts))
  mat <- if (!is.null(r) && nrow(r)) r$matrix[1] else NA_character_
  # a set gene absent from the matrix the score was computed on simply did not contribute
  present <- if (identical(mat, "curated")) genes
             else if (identical(mat, "broad")) (GENES_FULL %||% genes) else ALL_GENES
  gene_line <- function(nm, s) {
    g <- if (is.null(SCORE_SETS)) NULL else SCORE_SETS[[nm]]
    if (is.null(g)) return(tags$li(HTML(sprintf("<b>%s %s</b> &mdash; gene list unavailable", s, nm))))
    miss <- setdiff(g, present)
    tags$li(HTML(sprintf("<b>%s %s</b> (%d genes)<br><code>%s</code>%s", s, nm, length(g),
      paste(g, collapse = ", "),
      if (length(miss)) sprintf(paste0("<br><span style='color:#c62828'>not present in the %s ",
        "matrix, so it did not contribute: %s</span>"), mat, paste(miss, collapse = ", ")) else "")))
  }
  tagList(
    if (!is.null(purpose)) tags$p(HTML(sprintf("<b>What it is for.</b> %s", purpose))),
    tags$p(HTML(paste0("<b>How it is built.</b> ",
      if (length(parts) == 2) sprintf(paste0("A difference of two module scores, ",
        "<code>%s</code> minus <code>%s</code> &mdash; a difference, not a ratio. Each is "),
        parts[1], parts[2]) else "A single module score: ",
      "an <code>AddModuleScore</code> equivalent, the mean log-normalised expression of the set ",
      "minus the mean of 100 control genes drawn from the same expression bin (24 quantile bins, ",
      "<code>seed = 1</code>). Zero therefore means &ldquo;no higher than a random set of ",
      "equally-expressed genes&rdquo;. Raw values &mdash; no z-scoring or rescaling. ",
      if (!is.null(r) && nrow(r)) sprintf(paste0("Scored on the <b>%s</b> matrix; %s of %s set ",
        "genes found; %s cells scored."), mat, r$n_genes_used[1], r$n_genes_set[1],
        format(r$n_cells_scored[1], big.mark = ",")) else ""))),
    if (length(parts)) tags$div(tags$b("The genes."), tags$ul(style = "margin-top:4px",
      lapply(seq_along(parts), function(i) gene_line(parts[i], sgn[i])))),
    tags$p(style = "color:#555", HTML(paste0("These are <b>hand-curated canonical markers</b>, ",
      "not a database term &mdash; no MSigDB, no GO. There is no recorded source for the ",
      "specific choices; see the README section &ldquo;Module scores&rdquo; before citing them."))),
    if (is.null(SCORE_SETS)) div(style = "color:#c62828", paste(
      "Gene lists unavailable in this deploy — build_signature_scores.R is not next to app.R.")))
}
# every set, once, for the reference table under Help
score_sets_df <- function() {
  validate(need(!is.null(SCORE_SETS), paste(
    "Gene lists unavailable — build_signature_scores.R is not next to app.R in this deploy.")))
  used <- vapply(names(SCORE_SETS), function(nm) {
    if (is.null(SCOREMETA)) return("")
    paste(SCOREMETA$score[vapply(strsplit(SCOREMETA$sets, " - ", fixed = TRUE),
                                 function(p) nm %in% trimws(p), TRUE)], collapse = ", ")
  }, "")
  data.frame(set = names(SCORE_SETS), n_genes = unname(lengths(SCORE_SETS)),
             used_by = unname(used),
             genes = unname(vapply(SCORE_SETS, paste, "", collapse = ", ")),
             row.names = NULL, stringsAsFactors = FALSE)
}

# ---- gene-set Venn -----------------------------------------------------------
# Crosses any two or three of the sets the other tabs already produce. A Venn hides
# the three things that decide whether an overlap means anything -- the threshold that
# built each set, the direction of change, and the overlap expected by chance -- so
# every region here is reported alongside its null rather than left to look definitive.
#
# The cycling circle defaults to the CURATED canonical genes, not the data-driven axis.
# The axis calls 531 genes cycling-associated but only 46 are canonical; the rest include
# Ran, Nap1l1, Calm1, Ppia and other housekeeping genes, because cycling cells are
# globally more transcriptionally active and the AUC partly measures output. Labelled
# "cycling genes" on a Venn that would be read as far more than it says.
VN_CANONICAL <- unique(unlist(GENE_SETS[intersect(
  c("Cell cycle (S)", "Cell cycle (G2/M)", "E2F targets"), names(GENE_SETS))]))

vn_choices <- function() {
  de <- if (is.null(FG_CTAB)) character(0) else setNames(
    unlist(lapply(FG_CTAB$key, function(k) paste0("de:", k, c(":both", ":up", ":down")))),
    unlist(lapply(seq_len(nrow(FG_CTAB)), function(i)
      paste0(FG_CTAB$label[i], c("", " — up only", " — down only")))))
  ax <- c("Maturation: mature-associated"   = "ax:mat:mature",
          "Maturation: immature-associated" = "ax:mat:immature",
          "Metabolic: oxidative"            = "ax:met:oxidative",
          "Metabolic: glycolytic"           = "ax:met:glycolytic",
          "Cycling (data-driven axis)"      = "ax:cyc:cycling")
  iq <- c("Intersection: immature UP in P7 KO" = "iq:immature_up_in_KO",
          "Intersection: mature DOWN in P7 KO" = "iq:mature_down_in_KO")
  cur <- c(setNames("cur:__canonical__", "Canonical cell cycle (S + G2/M + E2F targets)"),
           setNames(paste0("cur:", names(GENE_SETS)), paste0("Curated: ", names(GENE_SETS))))
  list("Differential expression" = as.list(de), "Gene axes" = as.list(ax),
       "Intersection quadrants"  = as.list(iq), "Curated panels" = as.list(cur))
}
# Resolve one set id to its genes AND its testable space. The universe matters: a
# 52-gene curated panel crossed against a set drawn from 24k tested genes gives a
# meaningless hypergeometric unless the universe is the space both could have come from.
# What "up" and "down" mean, in one place. Both the Venn tab and the WT-programs tab
# build gene lists from the same DE tables; if they each spelled this out the two could
# drift and nobody would notice, because both would still return plausible lists.
#
# `measure` picks the effect size, and the choice is not cosmetic. presto's logFC is a
# difference of mean log-normalised expression, so it scales with how highly expressed a
# gene is. Cell-cycle genes in cardiomyocytes are sparse: between WT P0 and P7, Mcm3 goes
# from 6.6% to 27.2% detection and that is a log2FC of 0.16, while Myh7 moves 2.5. A
# symmetric |log2FC| cut therefore classifies maturation genes and can never classify a
# cell-cycle one -- it would report "no cell-cycle gene changes", which is false. AUC is
# rank-based and scale-free, so the same number means the same thing on both. This is the
# identical failure mode build_fourgroup.R hit when it classified its maturation axis, and
# it resolved it the same way (MAT_AUC = 0.60).
de_pass <- function(d, dir, padj_cut, eff_cut, measure = "lfc") {
  keep <- !d$confounder & !is.na(d$padj) & d$padj < padj_cut
  if (identical(measure, "auc")) {
    if (is.null(d$auc)) return(character(0))
    keep <- keep & !is.na(d$auc)
    if (identical(dir, "up"))   keep <- keep & d$auc >= eff_cut
    if (identical(dir, "down")) keep <- keep & d$auc <= 1 - eff_cut
    if (identical(dir, "both")) keep <- keep & (d$auc >= eff_cut | d$auc <= 1 - eff_cut)
  } else {
    keep <- keep & abs(d$log2FoldChange) >= eff_cut
    if (identical(dir, "up"))   keep <- keep & d$log2FoldChange > 0
    if (identical(dir, "down")) keep <- keep & d$log2FoldChange < 0
  }
  d$gene[which(keep)]
}
vn_set <- function(id, cluster, stratum, grid, padj_cut, lfc_cut, measure = "lfc") {
  if (is.null(id) || !nzchar(id) || identical(id, "none")) return(NULL)
  kind <- sub(":.*$", "", id); rest <- sub("^[^:]*:", "", id)
  ax_univ <- if (!is.null(GM)) GM$gene else ALL_GENES
  out <- function(g, lab, univ) list(genes = setdiff(unique(g), CONF),
                                     label = lab, universe = setdiff(unique(univ), CONF))
  if (kind == "de") {
    key <- sub(":.*$", "", rest); dir <- sub("^[^:]*:", "", rest)
    d <- fg_grid(grid)[[cluster]][[paste0(key, "__", stratum)]]
    ct <- fg_ct(key); lab <- paste0(if (is.null(ct)) key else ct$label,
                                    switch(dir, up = " (up)", down = " (down)", ""))
    # The testable space is every gene the matrix carries, NOT the genes in this table.
    # The DE tables are gated on expression AND a significance-or-effect condition, so an
    # expressed-but-unchanging gene is absent -- and those are exactly the genes a
    # hypergeometric universe must contain. Using the table as the universe would shrink
    # it to roughly the changing genes and inflate every expected overlap.
    univ <- if (identical(grid, "de2")) genes else (GENES_FULL %||% genes)
    if (is.null(d)) return(out(character(0), lab, univ))
    return(out(de_pass(d, dir, padj_cut, lfc_cut, measure), lab, univ))
  }
  if (kind == "ax") {
    validate(need(!is.null(GM), GM_MSG))
    axis <- sub(":.*$", "", rest); side <- sub("^[^:]*:", "", rest)
    col <- paste0(axis, "_class")
    validate(need(col %in% names(GM), paste0("The ", axis, " axis isn't in this data build.")))
    return(out(GM$gene[!is.na(GM[[col]]) & GM[[col]] == paste0(side, "-associated")],
               paste0(side, "-associated"), ax_univ))
  }
  if (kind == "iq") {
    validate(need(!is.null(FG$intersect), fg_int_msg))
    d <- FG$intersect[FG$intersect$cluster == cluster & !FG$intersect$confounder, , drop = FALSE]
    if (!nrow(d)) d <- FG$intersect[FG$intersect$cluster == "AllCM" & !FG$intersect$confounder, , drop = FALSE]
    return(out(d$gene[d$quadrant == rest], gsub("_", " ", rest), unique(d$gene)))
  }
  if (kind == "cur") {
    g <- if (identical(rest, "__canonical__")) VN_CANONICAL else GENE_SETS[[rest]]
    lab <- if (identical(rest, "__canonical__")) "canonical cell cycle" else rest
    return(out(intersect(g, ALL_GENES), lab, ALL_GENES))
  }
  NULL
}
# Collect the selected slots, and build the shared universe they must be judged in.
vn_sets <- function(ids, cluster, stratum, grid, padj_cut, lfc_cut, measure = "lfc") {
  ss <- Filter(Negate(is.null), lapply(ids, vn_set, cluster = cluster, stratum = stratum,
                                       grid = grid, padj_cut = padj_cut, lfc_cut = lfc_cut,
                                       measure = measure))
  validate(need(length(ss) >= 2, "Pick at least two gene sets."))
  univ <- Reduce(intersect, lapply(ss, `[[`, "universe"))
  validate(need(length(univ) > 0, "These sets have no shared testable space — the overlap statistics would be meaningless."))
  for (i in seq_along(ss)) ss[[i]]$genes <- intersect(ss[[i]]$genes, univ)
  list(sets = ss, universe = univ)
}
# Every disjoint region, keyed by which circles it belongs to ("A", "AB", "ABC", ...)
vn_regions <- function(sets) {
  n <- length(sets); nm <- LETTERS[seq_len(n)]
  g <- lapply(sets, `[[`, "genes")
  allg <- unique(unlist(g))
  memb <- vapply(g, function(x) allg %in% x, logical(length(allg)))
  if (is.null(dim(memb))) memb <- matrix(memb, ncol = n)
  key <- apply(memb, 1, function(r) paste(nm[r], collapse = ""))
  split(allg, key)
}
# Fixed 2- and 3-circle layouts, drawn with geom_polygon so no new package is needed.
vn_layout <- function(n) {
  if (n == 2) list(
    circles = data.frame(set = c("A","B"), x = c(-0.45, 0.45), y = c(0, 0), r = c(1, 1)),
    labels  = data.frame(key = c("A","B","AB"), x = c(-0.95, 0.95, 0), y = c(0, 0, 0)),
    titles  = data.frame(set = c("A","B"), x = c(-1.05, 1.05), y = c(1.15, 1.15)))
  else list(
    circles = data.frame(set = c("A","B","C"), x = c(0, -0.5, 0.5), y = c(0.55, -0.35, -0.35), r = 1),
    labels  = data.frame(
      key = c("A","B","C","AB","AC","BC","ABC"),
      x = c(0, -0.85, 0.85, -0.55, 0.55, 0, 0),
      y = c(1.05, -0.7, -0.7, 0.15, 0.15, -0.7, -0.2)),
    titles  = data.frame(set = c("A","B","C"), x = c(0, -1.25, 1.25), y = c(1.85, -1.35, -1.35)))
}
VN_PAL <- c(A = "#1565c0", B = "#c62828", C = "#f9a825")
vn_plot <- function(vs, bs = 13, ttl = NULL) {
  sets <- vs$sets; n <- length(sets); nm <- LETTERS[seq_len(n)]
  lay <- vn_layout(n); reg <- vn_regions(sets)
  th <- seq(0, 2 * pi, length.out = 181)
  poly <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(
    set = lay$circles$set[i],
    x = lay$circles$x[i] + lay$circles$r[i] * cos(th),
    y = lay$circles$y[i] + lay$circles$r[i] * sin(th))))
  lab <- lay$labels
  # reg[[k]] errors on a missing name for lists, so index defensively — an empty
  # region is normal (two sets can simply not overlap)
  lab$n <- vapply(lab$key, function(k) if (k %in% names(reg)) length(reg[[k]]) else 0L, 0L)
  tit <- lay$titles
  tit$txt <- vapply(seq_len(n), function(i) sprintf("%s\n(%d)", sets[[i]]$label,
                                                    length(sets[[i]]$genes)), "")
  ggplot() +
    geom_polygon(data = poly, aes(x, y, fill = set, group = set), alpha = .32, colour = "grey30") +
    geom_text(data = lab, aes(x, y, label = n), size = bs / 2.4, fontface = "bold") +
    geom_text(data = tit, aes(x, y, label = txt), size = bs / 3.4, fontface = "bold",
              lineheight = .95) +
    scale_fill_manual(values = VN_PAL, guide = "none") +
    coord_equal(clip = "off") + theme_void(base_size = bs) +
    theme(plot.margin = margin(24, 24, 24, 24)) +
    labs(title = ttl %||% "Gene-set overlap",
         caption = sprintf(paste("Universe: %s genes tested in common. Region sizes move with the",
                                 "padj and log2FC cuts in the sidebar.\nDescriptive only — n = 1 animal per group."),
                           format(length(vs$universe), big.mark = ",")))
}
# Observed vs expected for every pair — the null the picture cannot show on its own.
vn_stats <- function(vs) {
  sets <- vs$sets; N <- length(vs$universe); n <- length(sets)
  cmb <- utils::combn(seq_len(n), 2, simplify = FALSE)
  do.call(rbind, lapply(cmb, function(ij) {
    a <- sets[[ij[1]]]$genes; b <- sets[[ij[2]]]$genes
    o <- length(intersect(a, b)); e <- length(a) * length(b) / N
    data.frame(set_A = sets[[ij[1]]]$label, set_B = sets[[ij[2]]]$label,
               n_A = length(a), n_B = length(b), overlap = o,
               expected = round(e, 1), fold = round(o / max(e, 1e-9), 2),
               p_hypergeom = signif(stats::phyper(o - 1, length(b), N - length(b),
                                                  length(a), lower.tail = FALSE), 3),
               stringsAsFactors = FALSE)
  }))
}
# Region -> gene list, for the table and the CSV
vn_region_df <- function(vs) {
  sets <- vs$sets; nm <- LETTERS[seq_along(sets)]
  labs <- vapply(sets, `[[`, "", "label")
  reg <- vn_regions(sets)
  do.call(rbind, lapply(names(reg), function(k) {
    inn <- strsplit(k, "")[[1]]
    data.frame(region = k,
               sets = paste(labs[match(inn, nm)], collapse = " & "),
               n = length(reg[[k]]),
               genes = paste(sort(reg[[k]]), collapse = ", "),
               stringsAsFactors = FALSE)
  }))[order(-vapply(reg, length, 0L)), , drop = FALSE]
}


# ---- WT programs x P7 KO clusters --------------------------------------------
# The collaborator's four crossings, in one place. Two steps, then a crossing:
#   1. WT P0->P7, restricted to a curated gene category (maturation / cell cycle)
#   2. P7 KO-vs-WT, unioned over a named group of CM subclusters
#   3. each WT list crossed against the OPPOSITE category's cluster group
#
# The KO side is defined by CLUSTERS, not by a gene category -- that is the email's
# "8 clusters including G1/Maturation related genes (CM1, CM2, CM3, CM7, CM8), and
# Cell cycle related genes (CM2, 4, 5)". Filtering both sides by category instead
# would make crossings 1 and 2 intersect the maturation set with the cell-cycle set,
# and those two curated sets share only Mki67 and Top2a -- the Venns would read empty
# by construction rather than by biology.
XC_MAT_CLUSTERS <- c("CM1","CM2","CM3","CM7","CM8")   # "G1/Maturation related"
XC_CYC_CLUSTERS <- c("CM2","CM4","CM5")               # "Cell cycle related". CM2 is in BOTH.
XC_WT_KEY  <- "WT_P0_vs_P7"     # A = WT-P7, so log2FC > 0 means up at P7
XC_KO_KEY  <- "P7_KO_vs_WT"     # A = KO-P7, so log2FC > 0 means up in KO
XC_CANON   <- "__canonical__"
XC_COMPARISONS <- list(
  list(key = "matUP_cycDN", wt = "mat", wt_dir = "up",   ko = "cyc", ko_dir = "down"),
  list(key = "matDN_cycUP", wt = "mat", wt_dir = "down", ko = "cyc", ko_dir = "up"),
  list(key = "cycUP_matDN", wt = "cyc", wt_dir = "up",   ko = "mat", ko_dir = "down"),
  list(key = "cycDN_matUP", wt = "cyc", wt_dir = "down", ko = "mat", ko_dir = "up"))

# mt- genes are KO-up in all seven subclusters and KO-down in none, which is a library
# read-fraction difference rather than biology (analysis/2026-08-21_email). Left in,
# they dominate every KO-up union and carry the oxidative-phosphorylation terms with them.
xc_is_mt <- function(g) grepl("^mt-", g, ignore.case = TRUE)
# The testable space is the matrix, never the gene category: a hypergeometric against a
# 14-gene universe is not a null, it is a tautology.
xc_universe <- function(grid) setdiff(if (identical(grid, "de2")) genes else (GENES_FULL %||% genes), CONF)
xc_category <- function(set) {
  g <- if (identical(set, XC_CANON)) VN_CANONICAL else GENE_SETS[[set]]
  intersect(g %||% character(0), ALL_GENES)
}
xc_cat_label <- function(set) if (identical(set, XC_CANON)) "canonical cell cycle" else set
xc_eff_label <- function(p) if (identical(p$measure, "auc"))
  sprintf("AUC >= %s (or <= %s)", p$eff, signif(1 - p$eff, 3)) else
  sprintf("|log2FC| >= %s", p$eff)
# Both effect-size measures side by side for the four WT lists. Reported ALWAYS, not
# only on disagreement: "no cell-cycle gene changes between WT P0 and P7" is a conclusion
# someone would carry away from this tab, and on this data it is a property of the log2FC
# scale rather than a result.
xc_measure_audit <- function(p) {
  d <- fg_grid(p$grid)[[p$wt_cluster]][[paste0(XC_WT_KEY, "__", p$stratum)]]
  if (is.null(d)) return(NULL)
  cats <- list(maturation = p$mat_set, "cell cycle" = p$cyc_set)
  rows <- list()
  for (cn in names(cats)) {
    cg <- xc_category(cats[[cn]])
    for (dir in c("up", "down")) {
      na <- length(intersect(de_pass(d, dir, p$padj, p$eff_auc %||% 0.60, "auc"), cg))
      nl <- length(intersect(de_pass(d, dir, p$padj, p$eff_lfc %||% 0.25, "lfc"), cg))
      rows[[length(rows) + 1L]] <- data.frame(category = cn,
        direction = if (identical(dir, "up")) "up at P7" else "up at P0",
        n_AUC = na, n_log2FC = nl, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
# The four comparisons carry different per-cluster lfc_ columns (CM2/4/5 vs CM1/2/3/7/8),
# so a plain rbind fails. Union the columns and pad rather than dropping the extras.
xc_rbind <- function(lst) {
  lst <- Filter(function(d) !is.null(d) && nrow(d), lst)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (cn in setdiff(cols, names(d))) d[[cn]] <- NA
    d[, cols, drop = FALSE] }))
}
xc_cat_choices <- function()
  c(setNames(XC_CANON, "Canonical cell cycle (S + G2/M + E2F targets)"),
    setNames(names(GENE_SETS), names(GENE_SETS)))

# Step 1: WT P0->P7 within one curated category.
xc_wt_set <- function(set, dir, cluster, stratum, grid, padj_cut, eff_cut,
                      measure = "auc", hide_mt = TRUE) {
  cat_g <- xc_category(set)
  lab <- sprintf("WT P0→P7\n%s %s", xc_cat_label(set),
                 if (identical(dir, "up")) "up at P7" else "up at P0")
  d <- fg_grid(grid)[[cluster]][[paste0(XC_WT_KEY, "__", stratum)]]
  g <- if (is.null(d)) character(0) else intersect(de_pass(d, dir, padj_cut, eff_cut, measure), cat_g)
  if (hide_mt) g <- g[!xc_is_mt(g)]
  list(genes = setdiff(unique(g), CONF), label = lab, universe = xc_universe(grid),
       n_category = length(cat_g), cluster = cluster, dir = dir, set = set)
}

# Step 2: P7 KO-vs-WT per cluster, tidied. One row per gene x cluster that passes.
xc_ko_long <- function(clusters, dir, stratum, grid, padj_cut, eff_cut,
                       measure = "auc", hide_mt = TRUE) {
  rows <- lapply(clusters, function(cl) {
    d <- fg_grid(grid)[[cl]][[paste0(XC_KO_KEY, "__", stratum)]]
    if (is.null(d) || !nrow(d)) return(NULL)
    i <- which(d$gene %in% de_pass(d, dir, padj_cut, eff_cut, measure))
    if (!length(i)) return(NULL)
    data.frame(cluster = cl, gene = d$gene[i], log2FoldChange = d$log2FoldChange[i],
               auc = d$auc[i], padj = d$padj[i], pct_KO = d$pct_A[i], pct_WT = d$pct_B[i],
               n_KO = d$n_A[i], n_WT = d$n_B[i], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) return(NULL)
  if (hide_mt) out <- out[!xc_is_mt(out$gene), , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out[order(out$gene, out$cluster), , drop = FALSE]
}
# Which of the named clusters actually carry a table for this stratum. CM4 has no G1
# stratum for ANY contrast, so the G1 option silently shrinks the cycling group unless
# the tab says so -- hence this is reported, not inferred by the reader.
xc_ko_present <- function(clusters, stratum, grid)
  clusters[vapply(clusters, function(cl)
    !is.null(fg_grid(grid)[[cl]][[paste0(XC_KO_KEY, "__", stratum)]]), TRUE)]
# Arms that clear the 10-cell DE floor but are still too thin to trust.
xc_ko_thin <- function(clusters, stratum, grid) {
  thin <- FG$built$thin_cells %||% 50
  keep <- vapply(clusters, function(cl) {
    d <- fg_grid(grid)[[cl]][[paste0(XC_KO_KEY, "__", stratum)]]
    !is.null(d) && nrow(d) && min(d$n_A[1], d$n_B[1]) < thin }, TRUE)
  clusters[keep]
}

xc_ko_set <- function(clusters, dir, stratum, grid, padj_cut, eff_cut,
                      measure = "auc", minc = 1, hide_mt = TRUE) {
  lab <- sprintf("P7 KO-vs-WT\n%s in %s", if (identical(dir, "up")) "up in KO" else "up in WT",
                 paste(clusters, collapse = "/"))
  L <- xc_ko_long(clusters, dir, stratum, grid, padj_cut, eff_cut, measure, hide_mt)
  g <- if (is.null(L)) character(0) else { tb <- table(L$gene); names(tb)[tb >= minc] }
  list(genes = setdiff(unique(g), CONF), label = lab, universe = xc_universe(grid),
       long = L, clusters = clusters, present = xc_ko_present(clusters, stratum, grid),
       dir = dir, minc = minc)
}

# One comparison -> exactly the list(sets, universe) shape vn_plot/vn_stats/vn_regions
# already take, so the drawing and the hypergeometric are reused unchanged.
xc_comparison <- function(cmp, p) {
  wt_set <- if (identical(cmp$wt, "mat")) p$mat_set else p$cyc_set
  ko_cls <- if (identical(cmp$ko, "mat")) p$mat_clusters else p$cyc_clusters
  A <- xc_wt_set(wt_set, cmp$wt_dir, p$wt_cluster, p$stratum, p$grid, p$padj, p$eff,
                 p$measure, p$hide_mt)
  B <- xc_ko_set(ko_cls, cmp$ko_dir, p$stratum, p$grid, p$padj, p$eff, p$measure,
                 p$minc, p$hide_mt)
  univ <- intersect(A$universe, B$universe)
  A$genes <- intersect(A$genes, univ); B$genes <- intersect(B$genes, univ)
  list(sets = list(A, B), universe = univ, cmp = cmp)
}
xc_label <- function(cmp, p) {
  wt_set <- if (identical(cmp$wt, "mat")) p$mat_set else p$cyc_set
  ko_cls <- if (identical(cmp$ko, "mat")) p$mat_clusters else p$cyc_clusters
  sprintf("WT %s %s  ×  P7 KO %s in %s", xc_cat_label(wt_set),
          if (identical(cmp$wt_dir, "up")) "up" else "down",
          if (identical(cmp$ko_dir, "up")) "up" else "down", paste(ko_cls, collapse = "/"))
}

# Step 1 table: both categories x both directions, one frame.
xc_wt_df <- function(p) {
  d <- fg_grid(p$grid)[[p$wt_cluster]][[paste0(XC_WT_KEY, "__", p$stratum)]]
  validate(need(!is.null(d) && nrow(d), fg_skip_msg(p$wt_cluster, XC_WT_KEY, p$stratum, p$grid)))
  cats <- list(maturation = p$mat_set, "cell cycle" = p$cyc_set)
  rows <- list()
  for (cn in names(cats)) for (dir in c("up", "down")) {
    g <- intersect(de_pass(d, dir, p$padj, p$eff, p$measure), xc_category(cats[[cn]]))
    if (p$hide_mt) g <- g[!xc_is_mt(g)]
    g <- setdiff(g, CONF)
    if (!length(g)) next
    i <- match(g, d$gene)
    rows[[length(rows) + 1L]] <- data.frame(
      category = cn, gene_set = xc_cat_label(cats[[cn]]),
      direction = if (identical(dir, "up")) "up at P7" else "up at P0",
      gene = g, log2FoldChange = d$log2FoldChange[i], auc = d$auc[i], padj = d$padj[i],
      pct_P7 = d$pct_A[i], pct_P0 = d$pct_B[i], stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  validate(need(!is.null(out), paste0(
    "No gene in either category passes padj < ", p$padj, " and ", xc_eff_label(p),
    " in this cluster and stratum. The categories are small by design (a few dozen genes) ",
    "-- try the other effect-size measure, loosen the cut, or pick a subcluster instead of ",
    "all cardiomyocytes.")))
  out[order(out$category, out$direction, -abs(out$log2FoldChange)), , drop = FALSE]
}
# Counts behind Step 1, so an empty circle later is traceable to its list.
xc_wt_counts <- function(p) {
  d <- try(xc_wt_df(p), silent = TRUE)
  if (inherits(d, "try-error")) return(NULL)
  as.data.frame(table(category = d$category, direction = d$direction), stringsAsFactors = FALSE)
}

# Step 2 table: every named cluster, both groups, both directions.
xc_ko_df <- function(p) {
  grp <- list(maturation = p$mat_clusters, "cell cycle" = p$cyc_clusters)
  rows <- list()
  for (gn in names(grp)) for (dir in c("up", "down")) {
    L <- xc_ko_long(grp[[gn]], dir, p$stratum, p$grid, p$padj, p$eff, p$measure, p$hide_mt)
    if (is.null(L)) next
    L$cluster_group <- gn
    L$direction <- if (identical(dir, "up")) "up in KO" else "up in WT"
    rows[[length(rows) + 1L]] <- L
  }
  out <- do.call(rbind, rows)
  validate(need(!is.null(out), "No gene passes the current cuts in these clusters."))
  out <- out[, c("cluster_group","cluster","direction","gene","log2FoldChange","auc","padj",
                 "pct_KO","pct_WT","n_KO","n_WT")]
  out[order(out$cluster_group, out$direction, out$cluster, -abs(out$log2FoldChange)), , drop = FALSE]
}
# gene x cluster log2FC -- "which genes in each" in one glance rather than a filter.
xc_ko_pivot_df <- function(p) {
  d <- xc_ko_df(p)
  cls <- unique(c(p$mat_clusters, p$cyc_clusters))
  cls <- cls[order(suppressWarnings(as.integer(sub("^CM", "", cls))), cls)]
  g <- sort(unique(d$gene))
  out <- data.frame(gene = g, stringsAsFactors = FALSE)
  for (cl in cls) { v <- d[d$cluster == cl, , drop = FALSE]
    out[[cl]] <- if (nrow(v)) v$log2FoldChange[match(g, v$gene)] else NA_real_ }
  out$n_clusters <- rowSums(!is.na(out[, cls, drop = FALSE]))
  out[order(-out$n_clusters, out$gene), , drop = FALSE]
}

# The deliverable table: one row per gene per comparison, both sides' evidence on the
# same row so a shared gene can be judged without a second lookup.
xc_gene_df <- function(vs, p) {
  A <- vs$sets[[1]]; B <- vs$sets[[2]]
  allg <- sort(union(A$genes, B$genes))
  if (!length(allg)) return(NULL)
  wt <- fg_grid(p$grid)[[p$wt_cluster]][[paste0(XC_WT_KEY, "__", p$stratum)]]
  iw <- if (is.null(wt)) rep(NA_integer_, length(allg)) else match(allg, wt$gene)
  L  <- B$long
  # Summarise the KO side once, vectorised. A per-gene subset inside a loop is quadratic
  # and this table is rebuilt on every keystroke in the threshold boxes.
  ks <- if (is.null(L)) NULL else {
    o   <- L[order(L$gene, -abs(L$log2FoldChange)), , drop = FALSE]
    top <- o[!duplicated(o$gene), , drop = FALSE]
    cls <- vapply(split(L$cluster, L$gene), paste, "", collapse = ",")
    nc  <- lengths(split(L$cluster, L$gene))
    data.frame(gene = top$gene, ko_log2FC_top = top$log2FoldChange, ko_padj_top = top$padj,
               ko_cluster_top = top$cluster, ko_clusters = unname(cls[top$gene]),
               ko_n_clusters = as.integer(unname(nc[top$gene])), stringsAsFactors = FALSE)
  }
  ik <- if (is.null(ks)) rep(NA_integer_, length(allg)) else match(allg, ks$gene)
  na_if_null <- function(x, i, default) if (is.null(x)) rep(default, length(i)) else x[i]
  out <- data.frame(
    comparison = xc_label(vs$cmp, p),
    region = ifelse(allg %in% A$genes & allg %in% B$genes, "shared",
             ifelse(allg %in% A$genes, "WT only", "KO only")),
    gene = allg,
    in_maturation_set = allg %in% xc_category(p$mat_set),
    in_cellcycle_set  = allg %in% xc_category(p$cyc_set),
    wt_log2FC = na_if_null(wt$log2FoldChange, iw, NA_real_),
    wt_padj   = na_if_null(wt$padj,           iw, NA_real_),
    ko_clusters    = na_if_null(ks$ko_clusters,    ik, NA_character_),
    ko_n_clusters  = na_if_null(ks$ko_n_clusters,  ik, NA_integer_),
    ko_cluster_top = na_if_null(ks$ko_cluster_top, ik, NA_character_),
    ko_log2FC_top  = na_if_null(ks$ko_log2FC_top,  ik, NA_real_),
    ko_padj_top    = na_if_null(ks$ko_padj_top,    ik, NA_real_),
    stringsAsFactors = FALSE)
  # per-cluster log2FC, so the union never hides which cluster carried a gene
  for (cl in B$clusters) {
    v <- if (is.null(L)) NULL else L[L$cluster == cl, , drop = FALSE]
    out[[paste0("lfc_", cl)]] <-
      if (is.null(v) || !nrow(v)) NA_real_ else v$log2FoldChange[match(allg, v$gene)]
  }
  out$ko_n_clusters[is.na(out$ko_n_clusters)] <- 0L
  out[order(factor(out$region, c("shared","WT only","KO only")),
            -abs(out$wt_log2FC), -abs(out$ko_log2FC_top)), , drop = FALSE]
}


# ---- gene map: maturation axis x metabolic axis (FG$geneaxes) ----------------
# Each point is a GENE, placed by how strongly it marks mature-vs-immature CMs (x)
# and oxidative-vs-glycolytic metabolism (y). Both coordinates are AUCs from a
# tertile split computed within each timepoint and averaged, so neither axis is a
# restatement of P0-vs-P7.
GM        <- if (!is.null(FG)) FG$geneaxes else NULL
GM_MSG    <- "Gene map isn't in this data build — run build_fourgroup.R and redeploy."
# The split point is each axis's own MEDIAN, not 0.5: wilcoxauc's AUC carries a small
# global offset (the tertile groups differ in detection rate), and since most genes sit
# within ~0.02 of the median, splitting at 0.5 would pile 65% of them into one corner.
# One centre per panel. The AUC offset is a property of the cells being compared, so
# P0 and P7 do not share the pooled median -- here P0 sits at 0.498 and P7 at 0.519.
GM_CENTRES <- local({
  a <- if (!is.null(GM)) attr(GM, "centre") else NULL
  # older builds stored a bare vector rather than one centre per panel
  if (is.null(a)) list(avg = c(mat = 0.5, met = 0.5)) else if (is.list(a)) a else list(avg = a)
})
GM_PANELS <- names(GM_CENTRES)
# column pair backing each panel; "avg" uses the averaged AUCs
gm_cols <- function(panel) {
  if (identical(panel, "avg")) c("mat_auc", "met_auc")
  else paste0(c("mat_auc_", "met_auc_"), panel)
}
gm_centre <- function(panel) if (!is.null(GM_CENTRES[[panel]])) GM_CENTRES[[panel]] else GM_CENTRES[[1]]
GM_NCELL <- if (!is.null(GM)) attr(GM, "n_cells_per_tp") else NULL
# diagonal = the expected coupling (mature<->oxidative, immature<->glycolytic);
# off-diagonal = genes that uncouple the two programs, which is the interesting part.
GM_QUAD_PAL <- c("mature+oxidative"    = "#c62828", "immature+glycolytic" = "#1565c0",
                 "mature+glycolytic"   = "#ef6c00", "immature+oxidative"  = "#00838f")
GM_QUADS <- names(GM_QUAD_PAL)
# Below this distance a gene is close enough to the centre that which side it falls on
# is jitter — the P0 and P7 centres alone differ by 0.021 on the maturation axis.
GM_MOVE_MIN <- 0.05
# A gene counts as "confidently labelled" on an axis when it clears this margin from
# that panel's centre AND is significant in that panel. Same margin the intersection
# tab uses (its mat_class threshold is 0.60 against a 0.5 split), but measured from the
# panel's own centre, so it means the same thing on the P0 and P7 panels too. Most
# genes never clear it -- that is the point: it marks the ones worth trusting.
GM_CONF_MARGIN <- 0.10

# `panel` selects which pair of coordinates to plot: the timepoint-averaged axes, or
# P0 / P7 on their own. Quadrant and distance are recomputed against that panel's own
# centre rather than reused from the averaged one.
gm_df <- function(panel = "avg", quadrants = NULL, min_dist = 0,
                  hide_sets = TRUE, geneset = "__all__", conf_only = FALSE) {
  validate(need(!is.null(GM), GM_MSG))
  panel <- panel %||% "avg"
  cols <- gm_cols(panel)
  validate(need(all(cols %in% names(GM)), paste0(
    "The ", panel, " panel isn't in this data build — re-run build_fourgroup.R and redeploy.")))
  d <- GM
  d$x <- d[[cols[1]]]; d$y <- d[[cols[2]]]
  d <- d[!is.na(d$x) & !is.na(d$y), , drop = FALSE]
  validate(need(nrow(d), paste0("No genes have coordinates in the ", panel, " panel.")))
  ctr <- gm_centre(panel)
  d$quadrant <- paste0(ifelse(d$x >= ctr[["mat"]], "mature", "immature"), "+",
                       ifelse(d$y >= ctr[["met"]], "oxidative", "glycolytic"))
  d$distance <- round(sqrt((d$x - ctr[["mat"]])^2 + (d$y - ctr[["met"]])^2), 4)
  # which points are the confident ones — effect clear of the centre AND significant,
  # judged inside this panel rather than borrowed from the averaged axis
  pcol <- function(pre) if (identical(panel, "avg")) paste0(pre, "_padj")
                        else paste0(pre, "_padj_", panel)
  conf1 <- function(v, c0, pc) {
    if (!pc %in% names(d)) return(rep(NA, nrow(d)))
    pv <- d[[pc]]; !is.na(pv) & pv < 0.05 & abs(v - c0) >= GM_CONF_MARGIN
  }
  d$mat_confident <- conf1(d$x, ctr[["mat"]], pcol("mat"))
  d$met_confident <- conf1(d$y, ctr[["met"]], pcol("met"))
  d$confidence <- ifelse(d$mat_confident & d$met_confident, "both axes",
                  ifelse(d$mat_confident, "maturation only",
                  ifelse(d$met_confident, "metabolic only", "neither")))
  # genes inside the scoring sets sit at the extremes of their own axis by
  # construction; hidden by default so the map isn't just recovering its own inputs
  if (isTRUE(hide_sets)) d <- d[is.na(d$in_score_set), , drop = FALSE]
  if (!is.null(geneset) && geneset != "__all__")
    d <- d[d$gene %in% genes_for_set(geneset), , drop = FALSE]
  if (!is.null(quadrants) && length(quadrants)) d <- d[d$quadrant %in% quadrants, , drop = FALSE]
  if (!is.null(min_dist)) d <- d[d$distance >= min_dist, , drop = FALSE]
  if (isTRUE(conf_only)) d <- d[d$confidence != "neither", , drop = FALSE]
  validate(need(nrow(d), "No genes pass these filters — lower the distance cut, re-enable a quadrant, or untick “confidently labelled only”."))
  d[order(-d$distance), , drop = FALSE]
}
# table view: the columns worth reading, plus the per-timepoint AUCs so a gene's
# movement between P0 and P7 is visible without switching panels
gm_table <- function(d) {
  tps <- setdiff(GM_PANELS, "avg")
  cols <- intersect(c("gene","quadrant","confidence","distance","x","y",
                      paste0("mat_auc_", tps), paste0("met_auc_", tps),
                      "p7ko_log2FC","mat_class","met_class","in_score_set"), names(d))
  out <- d[, cols, drop = FALSE]
  names(out)[match(c("x","y"), names(out), nomatch = 0)] <-
    c("mat_auc","met_auc")[c("x","y") %in% names(out)]
  out
}
gm_plot_ly <- function(d, label_n = 20, highlight = NULL, panel = "avg") {
  validate(need(nrow(d), "No genes to plot."))
  ctr <- gm_centre(panel)
  cx <- unname(ctr[["mat"]]); cy <- unname(ctr[["met"]])
  d$quadrant <- factor(d$quadrant, levels = GM_QUADS)
  d$hover <- sprintf(
    "<b>%s</b><br>maturation AUC: %.3f<br>metabolic AUC: %.3f<br>%s<br>distance: %.3f%s%s",
    d$gene, d$x, d$y, as.character(d$quadrant), d$distance,
    ifelse(is.na(d$p7ko_log2FC), "", sprintf("<br>P7 KO vs WT log2FC: %+.2f", d$p7ko_log2FC)),
    ifelse(is.na(d$in_score_set), "", paste0("<br><i>in the ", d$in_score_set, " scoring set</i>")))
  rng <- function(v, c0) { m <- max(abs(v - c0), na.rm = TRUE) * 1.08; c(c0 - m, c0 + m) }
  xr <- rng(d$x, cx); yr <- rng(d$y, cy)
  ln <- function(x0, x1, y0, y1) list(type = "line", x0 = x0, x1 = x1, y0 = y0, y1 = y1,
                                      line = list(color = "grey65", width = 1))
  corner <- function(x, y, txt, col) list(x = x, y = y, text = txt, showarrow = FALSE,
    font = list(size = 10, color = col), xanchor = "center", yanchor = "middle", opacity = 0.75)
  p <- plot_ly(d, x = ~x, y = ~y, color = ~quadrant, colors = GM_QUAD_PAL,
        customdata = ~gene, text = ~hover, hovertemplate = "%{text}<extra></extra>",
        type = "scattergl", mode = "markers",
        marker = list(size = 6, opacity = 0.5, line = list(width = 0)),
        source = "gm_scatter") |>
    layout(
      title = list(text = paste0("Gene map — ",
        if (identical(panel, "avg")) "timepoints averaged" else paste0(panel, " only")),
        font = list(size = 13)),
      xaxis = list(title = "maturation association, AUC  (← immature | mature →)",
                   range = xr, zeroline = FALSE),
      yaxis = list(title = "metabolic association, AUC  (← glycolytic | oxidative →)",
                   range = yr, zeroline = FALSE),
      legend = list(title = list(text = ""), itemsizing = "constant"),
      shapes = list(ln(xr[1], xr[2], cy, cy), ln(cx, cx, yr[1], yr[2])),
      annotations = list(
        corner(mean(c(cx, xr[2])), yr[2], "mature + oxidative",    "#c62828"),
        corner(mean(c(xr[1], cx)), yr[1], "immature + glycolytic", "#1565c0"),
        corner(mean(c(cx, xr[2])), yr[1], "mature + glycolytic",   "#ef6c00"),
        corner(mean(c(xr[1], cx)), yr[2], "immature + oxidative",  "#00838f")),
      margin = list(t = 30))
  # ring the confidently-labelled points so it is visible at a glance which of the
  # 11k are the ones that actually clear significance — most do not
  cf <- d[!is.na(d$confidence) & d$confidence != "neither", , drop = FALSE]
  if (nrow(cf)) p <- add_trace(p, data = cf, x = ~x, y = ~y, type = "scattergl", mode = "markers",
    marker = list(size = 9, color = "rgba(0,0,0,0)", line = list(color = "#222", width = 1.2)),
    name = sprintf("confidently labelled (%d)", nrow(cf)),
    text = ~hover, hovertemplate = "%{text}<extra></extra>", inherit = FALSE)
  # name the genes furthest from the centre — that ranking IS the "distance" question
  if (!is.null(label_n) && label_n > 0) {
    lab <- head(d[order(-d$distance), , drop = FALSE], label_n)
    if (nrow(lab)) p <- add_annotations(p, data = lab, x = ~x, y = ~y, text = ~gene,
      showarrow = FALSE, yshift = 9, font = list(size = 9, color = "#333"), inherit = FALSE)
  }
  if (!is.null(highlight) && nzchar(highlight) && highlight %in% d$gene) {
    hd <- d[match(highlight, d$gene), , drop = FALSE]
    p <- add_trace(p, x = hd$x, y = hd$y, type = "scattergl", mode = "markers",
                   marker = list(size = 16, color = "rgba(0,0,0,0)", line = list(color = "#111", width = 3)),
                   name = "selected", showlegend = FALSE, hoverinfo = "skip", inherit = FALSE)
  }
  event_register(p, "plotly_click")
}
# Static ggplot twin of gm_plot_ly: same quadrant colours, centre lines, corner
# labels, confidence rings and top-N gene labels, no hover/click. This is what
# the vector downloads and the Figure Studio handoff export for the gene map.
gm_plot_gg <- function(d, label_n = 20, highlight = NULL, panel = "avg") {
  validate(need(nrow(d), "No genes to plot."))
  ctr <- gm_centre(panel); cx <- unname(ctr[["mat"]]); cy <- unname(ctr[["met"]])
  d$quadrant <- factor(d$quadrant, levels = GM_QUADS)
  rng <- function(v, c0) { m <- max(abs(v - c0), na.rm = TRUE) * 1.08; c(c0 - m, c0 + m) }
  xr <- rng(d$x, cx); yr <- rng(d$y, cy)
  corner <- data.frame(
    x = c(mean(c(cx, xr[2])), mean(c(xr[1], cx)), mean(c(cx, xr[2])), mean(c(xr[1], cx))),
    y = c(yr[2], yr[1], yr[1], yr[2]),
    lab = c("mature + oxidative", "immature + glycolytic",
            "mature + glycolytic", "immature + oxidative"),
    col = c("#c62828", "#1565c0", "#ef6c00", "#00838f"))
  p <- ggplot(d, aes(x, y, color = quadrant)) +
    geom_hline(yintercept = cy, color = "grey65", linewidth = .4) +
    geom_vline(xintercept = cx, color = "grey65", linewidth = .4) +
    geom_point(size = 1.4, alpha = .5, shape = 16) +
    scale_color_manual(values = GM_QUAD_PAL, drop = FALSE) +
    annotate("text", x = corner$x, y = corner$y, label = corner$lab,
             color = corner$col, size = 3, alpha = .75) +
    coord_cartesian(xlim = xr, ylim = yr) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Gene map — ",
           if (identical(panel, "avg")) "timepoints averaged" else paste0(panel, " only")),
         x = "maturation association, AUC  (← immature | mature →)",
         y = "metabolic association, AUC  (← glycolytic | oxidative →)", color = NULL)
  cf <- d[!is.na(d$confidence) & d$confidence != "neither", , drop = FALSE]
  if (nrow(cf)) p <- p + geom_point(data = cf, shape = 21, size = 2.2, stroke = .4,
                                    fill = NA, color = "#222")
  if (!is.null(label_n) && label_n > 0) {
    lab <- head(d[order(-d$distance), , drop = FALSE], label_n)
    if (nrow(lab)) p <- p + geom_text(data = lab, aes(label = gene), color = "#333",
                                      size = 2.6, vjust = -0.9)
  }
  if (!is.null(highlight) && nzchar(highlight) && highlight %in% d$gene) {
    hd <- d[match(highlight, d$gene), , drop = FALSE]
    p <- p + geom_point(data = hd, shape = 21, size = 4.5, stroke = 1.1, fill = NA, color = "#111") +
      geom_text(data = hd, aes(label = gene), color = "#111", fontface = "bold",
                size = 3.2, vjust = -1.2)
  }
  p
}
# where a gene sits, in words — shown next to the table so a single gene can be read
# off without hunting for its row
gm_gene_note <- function(gene, panel = "avg") {
  if (is.null(GM) || is.null(gene) || !nzchar(gene)) return(NULL)
  r <- GM[GM$gene == gene, , drop = FALSE]
  if (!nrow(r)) return(div(style = "font-size:13px;color:#777",
    paste0(gene, " is not on the gene map (it needs an association on both axes).")))
  # position in every panel, so "is this gene maturation-linked at P7 but not P0?" can be
  # read off directly instead of by flipping between panels
  one <- function(pn) {
    cols <- gm_cols(pn); if (!all(cols %in% names(r))) return(NULL)
    x <- r[[cols[1]]][1]; y <- r[[cols[2]]][1]
    if (is.na(x) || is.na(y)) return(sprintf("<b>%s</b>: not measured", pn))
    ct <- gm_centre(pn)
    q <- paste0(ifelse(x >= ct[["mat"]], "mature", "immature"), "+",
                ifelse(y >= ct[["met"]], "oxidative", "glycolytic"))
    dd <- sqrt((x - ct[["mat"]])^2 + (y - ct[["met"]])^2)
    sprintf("<b>%s</b>: %s &middot; mat %.3f, met %.3f &middot; dist %.3f",
            if (pn == "avg") "averaged" else pn, q, x, y, dd)
  }
  parts <- Filter(Negate(is.null), lapply(GM_PANELS, one))
  tps <- setdiff(GM_PANELS, "avg")
  qd <- function(pn) {
    cols <- gm_cols(pn); if (!all(cols %in% names(r))) return(NULL)
    x <- r[[cols[1]]][1]; y <- r[[cols[2]]][1]; if (is.na(x) || is.na(y)) return(NULL)
    ct <- gm_centre(pn)
    list(q = paste0(ifelse(x >= ct[["mat"]], "mature", "immature"), "+",
                    ifelse(y >= ct[["met"]], "oxidative", "glycolytic")),
         x = x, y = y, d = sqrt((x - ct[["mat"]])^2 + (y - ct[["met"]])^2))
  }
  qs <- Filter(Negate(is.null), setNames(lapply(tps, qd), tps))
  delta <- NULL; moved <- FALSE
  if (length(qs) == 2) {
    a <- qs[[1]]; b <- qs[[2]]
    delta <- sprintf("%s → %s: maturation %+.3f, metabolic %+.3f",
                     names(qs)[1], names(qs)[2], b$x - a$x, b$y - a$y)
    # A bare quadrant flip is NOT evidence of movement: the panel centres themselves
    # differ by 0.021 on the maturation axis while the median gene's own AUC moves
    # 0.027, so 66% of genes "flip" on centre-line jitter alone. Only call it a move
    # when the gene is clear of the centre in BOTH panels.
    moved <- a$q != b$q && min(a$d, b$d) >= GM_MOVE_MIN
  }
  div(style = "font-size:13px;margin-bottom:6px",
    HTML(paste0("<b>", gene, "</b><br>", paste(parts, collapse = "<br>"))),
    if (!is.null(delta)) div(style = "color:#555", HTML(delta)),
    local({
      cl <- c(r$mat_class[1], r$met_class[1])
      cl <- cl[!is.na(cl) & cl != "ns"]
      if (length(cl)) div(style = "color:#1b5e20",
        HTML(paste0("Confidently labelled on the averaged axis: ", paste(cl, collapse = ", ")))) else NULL
    }),
    if (moved) div(style = "color:#c62828;font-weight:600",
      HTML(sprintf("&#9888; sits on a different side at each age: %s",
                   paste(sprintf("%s %s", names(qs), vapply(qs, `[[`, "", "q")), collapse = " → ")))),
    if (!is.na(r$in_score_set[1])) div(style = "color:#c62828",
      HTML(paste0("In the ", r$in_score_set[1],
                  " scoring set — its position on that axis is partly circular."))))
}

# ---------------------------------------------------------------- UI ----------
# Wrap every plot output in a loading spinner (shown while the output computes /
# on tab switch), by shadowing the two output constructors used across the UI.
.spin <- function(x) shinycssloaders::withSpinner(x, type = 6, color = "#2c3e50", size = 0.6)
plotOutput   <- function(...) .spin(shiny::plotOutput(...))
plotlyOutput <- function(...) .spin(plotly::plotlyOutput(...))

ui <- page_navbar(
  title = "E2F7/8 heart scRNA-seq", theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = if (STUDIO_ON) studio_js(),

  # Navbar grouping, and one coupling worth knowing before rearranging it again:
  # tools/check_docs_coverage.py treats "^  nav_panel(" -- EXACTLY two spaces -- as a
  # top-level tab, which is what stops it demanding a chapter for all ~60 nested panels.
  # nav_menu() children are therefore deliberately NOT indented here. Re-indenting them
  # would make the docs gate stop seeing them, and every chapter claiming one would then
  # fail with "claims tab X, which is not a top-level tab in app.R".
  nav_menu("Whole heart",
  nav_panel("UMAP explorer", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("color_by", "Colour cells by",
                  choices = setNames(c("gene", CAT_COLS, CONT_COLS),
                                     c(labof("gene"), labof(CAT_COLS), labof(CONT_COLS))),
                  selected = "celltype"),
      conditionalPanel("input.color_by == 'gene'",
        selectInput("geneset", "Gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("gene", "Gene", choices = NULL, options = list(maxOptions = 50L))),
      selectInput("split", "Split panels by", c("(none)" = "none",
                  setNames(intersect(c("genotype","timepoint"), CAT_COLS),
                           labof(intersect(c("genotype","timepoint"), CAT_COLS))))),
      sliderInput("ptsize", "Point size", 1, 9, 4.5, 0.5),
      hr(), helpText("Hover a cell to highlight its whole cluster; double-click a legend entry to isolate one; ",
                     "single-click to toggle. Colour by a gene, then split by Genotype to compare KO vs WT."),
      accordion(open = FALSE, accordion_panel("Figure options",
        figure_controls("umap", export = "umap", base = FALSE, labels = TRUE)))),
    card(full_screen = TRUE, card_header(textOutput("umap_title")),
         dl_fig_ui("umapgg", "Download figure (static)"),
         plotlyOutput("umap", height = "640px")))),

  nav_panel("Gene detail", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("geneset2", "Gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
      selectizeInput("g2", "Gene", choices = NULL, options = list(maxOptions = 50L)),
      selectInput("grp", "Group by", setNames(CAT_COLS, labof(CAT_COLS)),
                  selected = if (has("celltype")) "celltype" else CAT_COLS[1]),
      selectInput("sp2", "Split by", c("(none)" = "none",
                  setNames(intersect(c("genotype","timepoint"), CAT_COLS),
                           labof(intersect(c("genotype","timepoint"), CAT_COLS)))), selected = "genotype"),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("Violin options",   figure_controls("vln", palette = TRUE,  rename = TRUE)),
        accordion_panel("Dot plot options", figure_controls("dot", palette = FALSE, rename = TRUE)))),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Expression distribution (violin)"), plotOutput("vln", height = "460px")),
      card(card_header("% expressing & mean expression"), plotOutput("dot", height = "460px"))),
    div(uiOutput("g2_cells"), style = "font-size:12px;color:#666;margin-top:4px"))),

  nav_panel("Composition", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("comp_fill", "Show fractions of", setNames(CAT_COLS, labof(CAT_COLS)),
                  selected = if (has("celltype")) "celltype" else CAT_COLS[1]),
      selectInput("comp_x", "Across groups", setNames(
                  intersect(c("orig.ident","genotype","timepoint"), names(meta)),
                  labof(intersect(c("orig.ident","genotype","timepoint"), names(meta)))),
                  selected = if (has("orig.ident")) "orig.ident" else "genotype"),
      accordion(open = FALSE, accordion_panel("Figure options",
        figure_controls("comp", palette = TRUE, rename = TRUE)))),
    card(full_screen = TRUE, card_header("Cell-type / state proportions"), plotOutput("comp", height = "560px")))),

  nav_panel("DE by cell type", layout_sidebar(
    sidebar = sidebar(width = 300,
      radioButtons("ct_tp", "Timepoint", c("P0","P7"), inline = TRUE),
      selectInput("ct_sel", "Cell type", choices = NULL),
      textInput("ct_search", "Filter genes (substring)", ""),
      checkboxInput("ct_hideconf", "Hide sex/construct genes", FALSE),
      volc_lfc_ui("ct_vlfc"),
      hr(), helpText("KO-vs-WT differential expression within each cell type.",
                     br(), strong("p-axis ranks candidates only — not valid at n = 1."))),
    navset_card_tab(
      nav_panel("Volcano + table",
        helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it and show its info below."),
        layout_columns(col_widths = c(6, 6),
          div(dl_fig_ui("ctvolc", "Download figure (static)"),
              plotlyOutput("ct_volcano", height = "470px")),
          div(uiOutput("ct_pick_ui"), dl_data_ui("ct_table"), DTOutput("ct_table", height = "440px"))),
        uiOutput("ct_geneinfo")),
      nav_panel("Heatmap (top genes × cell types)", dl_fig_ui("ctheat"),
        plotlyOutput("ct_heat", height = "620px"))))),

  nav_panel("Subset & DEGs", layout_sidebar(
    sidebar = sidebar(width = 320,
      tags$b("1. Filter cells"),
      lapply(CAT_COLS, function(c) selectizeInput(paste0("degf_", c), labof(c),
        choices = sort(unique(as.character(DMETA[[c]]))), multiple = TRUE,
        options = list(placeholder = "all"))),
      hr(), tags$b("2. Compare"),
      selectInput("deg_by", "Split groups by",
        c("Genotype (KO vs WT)" = "genotype", "Cell-cycle phase" = "Phase",
          "Cycling vs non-cycling" = "cycling", "Timepoint (P0 vs P7)" = "timepoint")),
      selectizeInput("deg_a", "Group A", choices = NULL, multiple = TRUE),
      selectizeInput("deg_b", "Group B", choices = NULL, multiple = TRUE),
      checkboxInput("deg_hideconf", "Hide sex/construct genes", FALSE),
      volc_lfc_ui("deg_vlfc"),
      actionButton("deg_run", "Compute DEGs", class = "btn-primary"),
      hr(), helpText("Descriptive Wilcoxon (presto) on log-norm expression of the ",
                     "filtered live cells. Hypothesis-generating only (n = 1); ",
                     "for rigorous KO-vs-WT use the precomputed DE tabs.")),
    div(textOutput("deg_n"), style = "font-size:13px;margin-bottom:4px"),
    helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it and show its info below."),
    layout_columns(col_widths = c(6, 6),
      div(dl_fig_ui("degvolc", "Download figure (static)"),
          plotlyOutput("deg_volcano", height = "470px")),
      div(uiOutput("deg_pick_ui"), dl_data_ui("deg_table"), DTOutput("deg_table", height = "440px"))),
    uiOutput("deg_geneinfo"))),

  nav_panel("Pathways & enrichment", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("enr_tp", "Timepoint", c("P0","P7"), selected = "P7"),
      selectInput("enr_ct", "Cell type", choices = NULL),
      hr(), helpText("Pre-computed pathway/GO/TF enrichment of the KO-vs-WT signal ",
                     "(fgsea Hallmark/KEGG/E2F, GO biological process, decoupleR TF activity). ",
                     strong("Descriptive only — n = 1.")),
      # export = "none": each plot already carries dl_fig_ui() above it, which holds the
      # download buttons and the Figure Studio handoff. These accordions add the options
      # those exports read, so what is projected, downloaded and sent to Figure Studio is
      # the same figure.
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("GSEA figure options",
          figure_controls("enrgsea", export = "none", palette = "updown", rename = FALSE,
                          default_base = 12, axis = TRUE, labelchars = TRUE)),
        accordion_panel("GO figure options",
          figure_controls("enrgo", export = "none", palette = "continuous", rename = FALSE,
                          default_base = 11, axis = TRUE, labelchars = TRUE, colourby = TRUE)),
        accordion_panel("TF figure options",
          figure_controls("enrtf", export = "none", palette = "continuous", rename = FALSE,
                          default_base = 11, axis = TRUE, labelchars = TRUE)))),
    navset_card_tab(
      wrapper = function(...) card_body(..., fillable = FALSE),
      nav_panel("GSEA pathways",
        dl_fig_ui("enrgsea"), plotlyOutput("enr_gsea_plot", height = "440px"),
        dl_data_ui("enr_gsea_tab"), DTOutput("enr_gsea_tab", height = "360px")),
      nav_panel("GO biological process",
        dl_fig_ui("enrgo"), plotlyOutput("enr_go_plot", height = "440px"),
        dl_data_ui("enr_go_tab"), DTOutput("enr_go_tab", height = "360px")),
      nav_panel("TF / regulon activity",
        helpText("E2F-family regulon activity across cell types (KO − WT), then the top TFs for the selected cell type."),
        dl_fig_ui("enre2f"), plotlyOutput("enr_e2f_heat", height = "380px"),
        dl_fig_ui("enrtf"), plotlyOutput("enr_tf_top", height = "460px"))))),

  nav_panel("Cell–cell signalling", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("cc_tp", "Timepoint", c("P0","P7"), selected = "P7"),
      selectInput("cc_pathway", "Pathway",
                  choices = c("All" = "All", setNames(commun_pathways(), commun_pathways())),
                  selected = if ("VEGF" %in% commun_pathways()) "VEGF" else "All"),
      radioButtons("cc_metric", "Show",
                   c("KO − WT (delta)" = "delta", "WT" = "WT", "KO" = "KO"), selected = "delta"),
      hr(), helpText("Curated ligand → receptor scoring (mean ligand in the sender × mean receptor ",
                     "in the receiver, on log-norm expression), focused on the E2F7/8 → Vegfa angiogenic axis. ",
                     strong("Not permutation-tested; descriptive, n = 1."), br(),
                     "A full CellChat/LIANA run would need the source Seurat objects.")),
    navset_card_tab(
      nav_panel("Sender × receiver heatmap",
        helpText("Interaction score aggregated (mean) over the ligand→receptor pairs in the chosen pathway."),
        dl_fig_ui("ccheat"), plotlyOutput("cc_heat", height = "540px")),
      nav_panel("Interaction table", dl_data_ui("cc_tab"), DTOutput("cc_tab", height = "520px")))))),

  nav_menu("Cardiomyocytes",
  nav_panel("Cardiomyocyte deep-dive", layout_sidebar(
    sidebar = sidebar(width = 320,
      conditionalPanel("input.cm_tabs == 'de'",
        selectInput("cm_sub", "Subcluster (for DE)", choices = NULL),
        # The four timepoint-specific contrasts come from app$fourgroup. subDE can only
        # ever offer the pooled one: its timepoint dimension was collapsed at build time
        # and cannot be recovered here, which is why splitting the map could never drive
        # this volcano.
        selectInput("cm_contrast", "Comparison",
                    choices = c("KO vs WT (P0 + P7 pooled)" = "pooled", fg_contrast_choices()),
                    selected = "pooled"),
        conditionalPanel("input.cm_contrast != 'pooled'",
          radioButtons("cm_stratum", "Cells used",
                       c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                       selected = "G1"),
          radioButtons("cm_grid", "DE matrix", choices = fg_grid_choices(), selected = "de")),
        checkboxInput("cm_hideconf", "Hide sex/construct genes (DE)", FALSE),
        volc_lfc_ui("cm_vlfc")),
      conditionalPanel("input.cm_tabs == 'subenr'",
        selectInput("cm_enr_contrast", "Comparison",
                    choices = c("KO vs WT (P0 + P7 pooled)" = "pooled", fg_contrast_choices()),
                    selected = "pooled"),
        conditionalPanel("input.cm_enr_contrast != 'pooled'",
          radioButtons("cm_enr_stratum", "Cells used",
                       c("All cells" = "all", "G1 only (phase-matched)" = "G1"), selected = "all"),
          radioButtons("cm_enr_ont", "GO ontology",
                       c("Biological process" = "BP", "Molecular function" = "MF",
                         "Cellular component" = "CC"), selected = "BP", inline = TRUE))),
      selectInput("cm_mapcolor", "Map: colour by",
                  c("Subcluster" = "subcluster", "Cell-cycle phase" = "Phase", "Cycling" = "cycling",
                    "Genotype" = "genotype", "Timepoint" = "timepoint", "Gene" = "gene",
                    # only offered when the flag has been built, so an older bundle does not
                    # present a colour-by that resolves to a missing column
                    if (!is.null(app$immune_contam)) c("Immune contamination" = "immune_contam"))),
      conditionalPanel("input.cm_mapcolor == 'gene'",
        selectInput("cm_geneset", "Gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("cm_gene", "Gene", choices = NULL, options = list(maxOptions = 50L))),
      selectInput("cm_map_split", "Split map by (res 0.2)",
                  c("(none)" = "none", "Genotype (WT|KO)" = "genotype",
                    "Timepoint (P0|P7)" = "timepoint", "Genotype × Timepoint" = "both")),
      conditionalPanel("input.cm_tabs == 'bars'",
        radioButtons("cm_bar_mode", "Y axis", c("Proportion" = "prop", "Count" = "count"), inline = TRUE)),
      conditionalPanel("input.cm_tabs == 'subenr'",
        accordion(open = FALSE,
          accordion_panel("Enrichment figure options",
            # One shared prefix for all four views: they are the same kind of plot and a
            # reader wants them consistent. register_fig()'s opts_prefix lets each keep its
            # own download id while reading these controls.
            figure_controls("cmsubenr", export = "none", palette = "continuous",
                            rename = FALSE, default_base = 11, axis = TRUE, labelchars = TRUE,
                            colourby = TRUE)))),
      conditionalPanel("input.cm_tabs == 'variant'",
        selectInput("clu_var", "Clustering variant", choices = NULL),
        radioButtons("clu_mat", "Matrix for live markers",
                     c("Broad (more genes, fewer cells)" = "deg", "Curated panel (all CM cells)" = "curated"),
                     selected = "deg"),
        selectInput("clu_cl", "Subcluster (for DE / enrichment)", choices = NULL),
        helpText(style = "font-size:12px",
                 "Every variant runs the identical SCTransform, PCA and Harmony — only the",
                 "number of PCs carried into the neighbour graph and UMAP differs.")),
      conditionalPanel("input.cm_tabs == 'objtest'",
        selectInput("objtest_res", "Resolution", choices = NULL),
        helpText(style = "font-size:12px",
                 "Same cells, same pipeline, same seeds — only how the object was built differs.")),
      hr(), helpText("True re-clustering of cardiomyocytes. Explore subcluster identity,",
                     "differential expression per subcluster, and cell-cycle state.",
                     br(), br(),
                     strong("Split map by"), " facets the map only. To split the ",
                     strong("DE"), " by genotype and timepoint use the ", strong("Comparison"),
                     " dropdown — it carries KO-vs-WT at P0 and at P7 separately, and each ",
                     "genotype's own P0-vs-P7, per subcluster.",
                     br(), br(),
                     "Split the map by timepoint to see whether two cycling subclusters ",
                     "separate by P0/P7 or by S vs G2/M phase."),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("Cell-cycle figure options",
          figure_controls("cmphase", palette = TRUE, rename = TRUE)),
        accordion_panel("Composition figure options",
          figure_controls("cmbar", palette = TRUE, rename = FALSE)))),
    # fillable = FALSE on every panel of this tabset, and it is load-bearing.
    #
    # page_navbar() defaults fillable = TRUE, so each nav_panel body is a flex fill
    # container and its children are fill items with min-height: 0. In panels with one
    # child that is what you want -- the plot uses the viewport. These panels have four or
    # five children (help text, a status alert, a control row, then the tabset holding the
    # plot), and flex divides the available height among them. On a short or narrow window
    # the plot's share collapses to a sliver, and because card_body also sets
    # overflow: auto you get a ~60 px box with its own scrollbar rather than a figure.
    #
    # An explicit height= on the plot does NOT win that argument, which is why raising it
    # from 340 to 420 px changed nothing. Opting the panel out of filling makes the
    # heights authoritative again and lets the page scroll instead of compressing.
    navset_card_tab(id = "cm_tabs",
      wrapper = function(...) card_body(..., fillable = FALSE),
      nav_panel("Subcluster map",
        helpText("Hover any cell to highlight all cells of its subcluster; move off to restore the full map.",
                 br(), "When split (res 0.2) is on, panels are coloured by subcluster and “colour by” is ignored."),
        dl_fig_ui("cmmap", "Download figure (static)"),
        plotlyOutput("cm_map", height = "600px")),
      nav_panel("Identity (marker heatmap)", dl_fig_ui("cmmarker"),
        plotlyOutput("cm_markerheat", height = "660px")),
      nav_panel("DE (per subcluster)", value = "de",
        helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it and show its info below."),
        uiOutput("cm_de_note"),
        fluidRow(
          column(6, dl_fig_ui("cmvolc", "Download figure (static)"),
                    plotlyOutput("cm_volcano", height = "440px")),
          column(6, uiOutput("cm_pick_ui"), dl_data_ui("cm_detab"), DTOutput("cm_detab"))),
        uiOutput("cm_geneinfo"),
        div(class = "mt-3",
            h5("DE heatmap — top genes × subclusters"),
            dl_fig_ui("cmlfcheat"), plotlyOutput("cm_lfcheat", height = "560px"),
            uiOutput("cm_lfcheat_note"))),
      nav_panel("Cell cycle", uiOutput("cm_icn_note"), plotOutput("cm_phase", height = "560px")),
      nav_panel("Composition (stacked bars)", value = "bars",
        helpText("Composition of each res-0.2 subcluster, broken down four ways — genotype (WT/KO), ",
                 "timepoint (P0/P7), cell-cycle phase, and cycling status. Each is its own plot; scroll to see all four."),
        dl_fig_ui("cmbargeno", "Download this panel"), plotOutput("cm_bar_geno",  height = "300px"),
        dl_fig_ui("cmbartp", "Download this panel"), plotOutput("cm_bar_tp",    height = "300px"),
        dl_fig_ui("cmbarphase", "Download this panel"), plotOutput("cm_bar_phase", height = "300px"),
        dl_fig_ui("cmbarcyc", "Download this panel"), plotOutput("cm_bar_cyc",   height = "300px")),
      nav_panel("Per-cluster summary", value = "summary",
        helpText("One row per res-0.2 subcluster: identity, composition, top marker & KO-vs-WT genes, and top pathways — scan all clusters without clicking through tabs. Sizes/percentages are from the displayed (sampled) cells."),
        # n_KO_samp / n_WT_samp read as genotype counts and are not. Without this note a
        # reader sees CM12 at n_KO_samp = 0 beside a composition bar showing it 38 % KO and
        # reasonably concludes one of them is wrong; both are right.
        helpText(style = "font-size:12px",
                 HTML(paste0("<b>n_KO_samp / n_WT_samp are pseudobulk SAMPLES, not cells.</b> ",
                 "The KO-vs-WT test aggregates cells into one sample per library &times; lane ",
                 "(eight possible) and keeps only samples with &ge; 20 cells, then needs two per ",
                 "genotype. So <code>n_KO_samp = 0</code> means no KO sample reached 20 cells &mdash; ",
                 "not that the subcluster has no KO cells. CM12 is the clear case: 36 KO and 63 WT ",
                 "cells, but spread ~9 per KO sample, so none clear the floor and the row is ",
                 "<code>skipped_too_few_or_unbalanced</code>."))),
        DTOutput("cm_summary")),
      nav_panel("Top markers", value = "topmarkers",
        helpText("One row per res-0.2 subcluster: top identity markers (by z-scored mean expression) plus ",
                 "each cluster's top identity GO term, top KO-vs-WT GSEA pathway, and top KO-up / KO-down genes — ",
                 "a quick read on what each subcluster is doing biologically."),
        dl_data_ui("cm_topmarkers"), DTOutput("cm_topmarkers")),
      nav_panel("Variant explorer", value = "variant",
        uiOutput("clu_banner"),
        navset_pill(
          nav_panel("Composition & phase",
            div(style = "margin-top:10px", uiOutput("clu_sum_note")),
            layout_columns(col_widths = c(6, 6),
              plotOutput("clu_comp",  height = "420px"),
              plotOutput("clu_phase", height = "420px"))),
          nav_panel("Markers (computed live)",
            div(style = "margin-top:10px", uiOutput("clu_mk_note")),
            dl_data_ui("clu_mk"), DTOutput("clu_mk")),
          nav_panel("KO vs WT (precomputed)",
            div(style = "margin-top:10px", uiOutput("clu_de_note")),
            dl_data_ui("clu_de"), DTOutput("clu_de")),
          nav_panel("Enrichment (precomputed)",
            div(style = "margin-top:10px", uiOutput("clu_enr_note")), DTOutput("clu_enr")),
          nav_panel("Cell cycle", DTOutput("clu_cyc")))),
      nav_panel("Object-mode test", value = "objtest",
        uiOutput("objtest_verdict"),
        dl_fig_ui("objtestmap", "Download figure (static)"),
        plotOutput("objtest_map", height = "420px"),
        div(style = "margin-top:10px", uiOutput("objtest_note")),
        dl_data_ui("objtest_tab"),
        DTOutput("objtest_tab")),
      nav_panel("Subcluster enrichment", value = "subenr",
        helpText("Per res-0.2 subcluster: identity markers and the differential signal, enriched. ",
                 "The ", strong("Comparison"), " dropdown switches between the pooled ",
                 "KO-vs-WT enrichment and the timepoint-specific four-group contrasts; the two ",
                 "directions are enriched ", strong("separately"), ". ",
                 strong("Descriptive only — n = 1.")),
        uiOutput("cm_enr_note"),
        div(class = "d-flex align-items-center gap-3 mb-2",
            radioButtons("cm_enr_mode", NULL, c("Single cluster" = "one", "All clusters" = "all"),
                         selected = "one", inline = TRUE),
            conditionalPanel("input.cm_enr_mode == 'one'",
              selectInput("cm_enr_sub", NULL, choices = NULL, width = "260px"))),
        conditionalPanel("input.cm_enr_mode == 'one'",
          navset_card_tab(
            wrapper = function(...) card_body(..., fillable = FALSE),
            nav_panel("Identity GO",
              dl_fig_ui("cmsubidgo"), plotlyOutput("cm_sub_idgo_plot", height = "440px"),
              dl_data_ui("cm_sub_idgo_tab"), DTOutput("cm_sub_idgo_tab", height = "320px")),
            nav_panel("GO — up",
              dl_fig_ui("cmsubkogo"), plotlyOutput("cm_sub_kogo_plot", height = "440px"),
              dl_data_ui("cm_sub_kogo_tab"), DTOutput("cm_sub_kogo_tab", height = "320px")),
            nav_panel("GO — down",
              dl_fig_ui("cmsubkodn"), plotlyOutput("cm_sub_kodn_plot", height = "440px"),
              dl_data_ui("cm_sub_kodn_tab"), DTOutput("cm_sub_kodn_tab", height = "320px")),
            nav_panel("GSEA",
              dl_fig_ui("cmsubgsea"), plotlyOutput("cm_sub_gsea_plot", height = "440px"),
              dl_data_ui("cm_sub_gsea_tab"), DTOutput("cm_sub_gsea_tab", height = "320px")))),
        conditionalPanel("input.cm_enr_mode == 'all'",
          div(class = "d-flex align-items-center gap-3 mb-2",
              radioButtons("cm_enr_perrow", "Plots per row", c("1" = "1", "2" = "2"),
                           selected = "1", inline = TRUE)),
          navset_card_tab(
            wrapper = function(...) card_body(..., fillable = FALSE),
            nav_panel("Identity GO", uiOutput("cm_grid_idgo")),
            nav_panel("GO — up",     uiOutput("cm_grid_kogo")),
            nav_panel("GO — down",   uiOutput("cm_grid_kodn")),
            nav_panel("GSEA",        uiOutput("cm_grid_gsea")))))))),

  nav_panel("E2F focus", layout_sidebar(
    sidebar = sidebar(width = 320,
      selectInput("e2f_ct", "Cell type",
                  choices = c("All cells" = "All",
                              setNames(sort(unique(as.character(meta$celltype))),
                                       gsub("_", " ", sort(unique(as.character(meta$celltype)))))),
                  selected = if ("Cardiomyocyte" %in% meta$celltype) "Cardiomyocyte" else "All"),
      selectInput("e2f_set", "Downstream gene set", choices = names(GENE_SETS), selected = "E2F targets"),
      hr(), helpText("E2f7/E2f8 are atypical ", strong("repressors"), " of the E2F cell-cycle program, ",
                     "so loss should de-repress (raise) their downstream targets — see the fold-change tab. ",
                     strong("Descriptive only — n = 1, sex-confounded.")),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("E2f7/8 plot options", figure_controls("e2f", palette = TRUE, rename = TRUE)),
        accordion_panel("Fold-change options", figure_controls("e2ffc", palette = FALSE, rename = FALSE)))),
    navset_card_tab(
      nav_panel("E2f7 / E2f8", plotOutput("e2f_expr", height = "560px")),
      nav_panel("Downstream targets — KO vs WT", plotOutput("e2f_fc", height = "560px"))))),

  nav_panel("Four-group (WT/KO × P0/P7)", layout_sidebar(
    sidebar = sidebar(width = 320,
      conditionalPanel("input.fg_tabs == 'de'",
        selectInput("fg_cluster", "Subcluster", choices = NULL),
        selectInput("fg_contrast", "Comparison", choices = NULL),
        radioButtons("fg_grid", "DE matrix", choices = fg_grid_choices(), selected = "de"),
        radioButtons("fg_stratum", "Cells used",
                     c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                     selected = "G1"),
        checkboxInput("fg_hideconf", "Hide sex/construct genes", FALSE),
        volc_lfc_ui("fg_vlfc"),
        hr(),
        # The "answer the email in one click" button: the whole contrast as a
        # workbook, rather than downloading each subcluster's table by hand.
        checkboxInput("fg_book_all", "Include every subcluster", TRUE),
        div(downloadButton("fg_book", "Download contrast workbook (XLSX)",
                           class = "btn-sm btn-outline-primary"), style = "margin-bottom:4px"),
        helpText(tags$small("One sheet per subcluster for the selected comparison and cell set, ",
                            "plus group sizes and a README carrying the caveats."))),
      conditionalPanel("input.fg_tabs == 'counts'",
        radioButtons("fg_count_mode", "Y axis",
                     c("% of subcluster" = "prop", "Cell count" = "count"), inline = TRUE)),
      conditionalPanel("input.fg_tabs == 'enr'",
        selectInput("fg_enr_cluster", "Subcluster", choices = NULL),
        selectInput("fg_enr_contrast", "Comparison", choices = NULL),
        radioButtons("fg_enr_stratum", "Cells used",
                     c("All cells" = "all", "G1 only (phase-matched)" = "G1"), selected = "all"),
        radioButtons("fg_enr_ont", "GO ontology",
                     c("Biological process" = "BP", "Molecular function" = "MF",
                       "Cellular component" = "CC"), selected = "BP", inline = TRUE),
        numericInput("fg_enr_topn", "Terms to plot", 20, 5, 60, 5),
        hr(),
        div(downloadButton("fg_enr_book", "Download enrichment workbook (XLSX)",
                           class = "btn-sm btn-outline-primary"), style = "margin-bottom:4px"),
        helpText(tags$small("Every subcluster for this comparison: GO up, GO down, GSEA and the ",
                            "coverage audit, one sheet each."))),
      conditionalPanel("input.fg_tabs == 'g1'",
        selectizeInput("fg_g1_clusters", "Subclusters", choices = NULL, multiple = TRUE),
        selectInput("fg_score", "Maturation / state score", choices = NULL),
        radioButtons("fg_score_stratum", "Cells used",
                     c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                     selected = "all"),
        dl_data_ui("fg_scores", "Per-cell score summary")),
      hr(),
      helpText(strong("Sort caveat. "), FG_SORT_NOTE),
      conditionalPanel("input.fg_tabs == 'de'",
        helpText(strong("DE matrix. "),
                 "Two grids over the same contrasts. ", strong("Gene coverage"), " tests ~24k genes but ",
                 "on a downsampled 8k cells, so the thinnest arms drop out — CM2's KO-P0 falls to 9 cells ",
                 "and CM4/CM9 lose their G1 strata. ", strong("Cell coverage"), " keeps all 30k cells so ",
                 "every contrast runs, but only over the 2,181-gene curated panel. Neither wins outright; ",
                 "if a contrast is missing here it will say whether the other matrix has it.")),
      helpText(strong("Descriptive only — n = 1 animal per group."),
               " Wilcoxon is run cell-level, so p-values are pseudoreplicated;",
               " tables are ranked by effect size, not by p."),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("Composition figure options",
          figure_controls("fgcount", palette = FALSE, rename = FALSE)),
        accordion_panel("Phase figure options",
          figure_controls("fgphase", palette = FALSE, rename = FALSE)),
        accordion_panel("Score figure options",
          figure_controls("fgscore", palette = FALSE, rename = FALSE)))),
    navset_card_tab(id = "fg_tabs",
      nav_panel("Group sizes", value = "counts",
        helpText("Cell counts and percentages for the four groups in every res-0.2 CM subcluster. ",
                 "The ", strong("underpowered"), " column names any arm too small to support DE ",
                 "— those contrasts are skipped rather than silently reported."),
        plotOutput("fg_counts_plot", height = "380px"),
        dl_data_ui("fg_counts_tab"), DTOutput("fg_counts_tab")),
      nav_panel("Four-group DE", value = "de",
        helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it."),
        uiOutput("fg_de_note"),
        fluidRow(
          column(6, dl_fig_ui("fgvolc", "Download figure (static)"),
                    plotlyOutput("fg_volcano", height = "440px")),
          column(6, uiOutput("fg_pick_ui"), dl_data_ui("fg_detab"), DTOutput("fg_detab"))),
        uiOutput("fg_geneinfo")),
      nav_panel("G1 & maturation", value = "g1",
        helpText("Top: cell-cycle phase composition per group (G1 % printed on each bar). ",
                 "Bottom: per-cell maturation / state scores. ",
                 "The question is whether P7 KO cardiomyocytes sit at a less mature score than P7 WT ",
                 "— compare within the G1 stratum to hold cycling composition fixed."),
        plotOutput("fg_phase_plot", height = "420px"),
        plotOutput("fg_score_plot", height = "440px"),
        div(class = "mt-3",
          h5("Is the P7 KO less mature, and more cycling?"),
          helpText("The KO−WT gap on both readouts, at each age. The question is whether the ",
                   "P7 gap is larger than the P0 gap — i.e. whether the effect is specific to P7 ",
                   "rather than present throughout. Read it in the G1 stratum so cycling composition ",
                   "is held fixed and the gap cannot be the FACS enrichment."),
          dl_fig_ui("fgsumm"), plotOutput("fg_summary_plot", height = "420px"),
          dl_data_ui("fg_summary_tab"), DTOutput("fg_summary_tab"))),
      nav_panel("Enrichment", value = "enr",
        helpText("GO and pathway enrichment for the contrast selected in the sidebar, with the ",
                 "two directions enriched ", strong("separately"), ". ",
                 "This is not the same as the ", strong("Subcluster enrichment"), " tab under the ",
                 "cardiomyocyte deep-dive, which pools P0 and P7 together — these are the ",
                 "timepoint-specific contrasts.",
                 br(), strong("Descriptive only — n = 1.")),
        uiOutput("fg_enr_note"),
        navset_card_tab(
          nav_panel("GO — up",
            dl_fig_ui("fgenrup"), plotlyOutput("fg_enr_up_plot", height = "440px"),
            dl_data_ui("fg_enr_up_tab"), DTOutput("fg_enr_up_tab", height = "320px")),
          nav_panel("GO — down",
            dl_fig_ui("fgenrdn"), plotlyOutput("fg_enr_dn_plot", height = "440px"),
            dl_data_ui("fg_enr_dn_tab"), DTOutput("fg_enr_dn_tab", height = "320px")),
          nav_panel("GSEA",
            dl_fig_ui("fgenrgsea"), plotlyOutput("fg_enr_gsea_plot", height = "440px"),
            dl_data_ui("fg_enr_gsea_tab"), DTOutput("fg_enr_gsea_tab", height = "320px")),
          nav_panel("Coverage audit",
            helpText("One row per cluster × direction × ontology: how many genes went in, how big ",
                     "the universe was, how many terms came out, and which selection rule fired. ",
                     "Read this before concluding a direction is uninformative — an empty result ",
                     "can mean the list was too small to test."),
            dl_data_ui("fg_enr_audit_tab"), DTOutput("fg_enr_audit_tab", height = "420px"))))))),

  nav_panel("Maturation ∩ P7 KO", layout_sidebar(
    sidebar = sidebar(width = 320,
      conditionalPanel("input.mi_tabs != 'candidates'",
        selectInput("mi_cluster", "Subcluster", choices = NULL),
        checkboxGroupInput("mi_quad", "Show quadrants",
          c("Immature genes UP in P7 KO" = "immature_up_in_KO",
            "Mature genes DOWN in P7 KO" = "mature_down_in_KO",
            "Immature DOWN in KO" = "immature_down_in_KO",
            "Mature UP in KO" = "mature_up_in_KO"),
          selected = c("immature_up_in_KO","mature_down_in_KO")),
        checkboxInput("mi_hideconf", "Hide sex/construct genes", TRUE)),
      conditionalPanel("input.mi_tabs == 'candidates'",
        selectInput("mi_source", "Candidate source", choices = FG_CAND_SOURCES, selected = "shortlist"),
        conditionalPanel("input.mi_source != 'shortlist' && input.mi_source != 'intersect'",
          numericInput("mi_topn", "Top N per cluster", 20, 5, 100, 5)),
        selectInput("mi_geneset", "Restrict to gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("mi_genes", "Genes", choices = NULL, multiple = TRUE,
                       options = list(maxOptions = 50L)),
        actionLink("mi_reset_genes", "reset to the shortlist"),
        selectizeInput("mi_cand_clusters", "Subclusters", choices = NULL, multiple = TRUE),
        radioButtons("mi_cand_stratum", "Cells used",
                     c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                     selected = "all"),
        radioButtons("mi_spec_grid", "DE matrix (specificity table)",
                     choices = fg_grid_choices(), selected = "de"),
        dl_data_ui("mi_cand", "Expression grid behind the plot")),
      hr(),
      helpText(strong("Maturation axis. "),
               "Genes are ranked by comparing the most- vs least-mature cardiomyocytes, ",
               "within each timepoint and then averaged, so the axis is maturation and not P0-vs-P7. ",
               "It uses the ", strong("cycle-free"), " maturation score (Mki67 / Top2a / Ccnd1 removed) ",
               "— otherwise “less mature ⇒ more cycling” would be partly circular."),
      helpText(strong("Cycling link. "),
               "The table also carries each gene's cell-cycle association — how strongly it ",
               "marks cycling (S/G2M) over non-cycling cardiomyocytes — and ",
               strong("cyc_resid"), ", the part of that not already explained by its maturation ",
               "position. The residual is the honest version: mature and cycling are ",
               "anti-correlated states, so raw cycling association just restates the maturation ",
               "axis. On this data the hypothesis-quadrant genes sit below the line, not above it."),
      helpText(strong("Descriptive only — n = 1 animal per group.")),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("Quadrant figure options",
          figure_controls("miquad", palette = FALSE, rename = FALSE)),
        accordion_panel("Candidate figure options",
          figure_controls("micand", palette = FALSE, rename = FALSE)))),
    navset_card_tab(id = "mi_tabs",
      nav_panel("Quadrant map", value = "quadrant",
        helpText("x: how strongly a gene marks mature (right) vs immature (left) cardiomyocytes. ",
                 "y: its P7 KO-vs-WT log2 fold change. The two shaded quadrants are the ones ",
                 "consistent with P7 KO cells being held in a less mature state."),
        plotOutput("mi_quadrant", height = "600px")),
      nav_panel("Intersection table", value = "table",
        uiOutput("mi_cyc_note"),
        dl_data_ui("mi_table"), DTOutput("mi_table")),
      nav_panel("Candidate genes", value = "candidates",
        helpText("Any gene, across CM subclusters × the four groups. ",
                 "Size = % of cells expressing, colour = mean expression. ",
                 "The table below reads the KO effect off the precomputed contrasts: ",
                 strong("P7_specificity"), " > 0 means the KO effect is larger at P7 than at P0, and ",
                 strong("priority_concentration"), " > 0 means it is larger inside CM2/CM4/CM5 than outside. ",
                 strong("state_maturation"), " and ", strong("state_cycling"),
                 " say whether the gene marks a less mature or a more cycling state, from the same ",
                 "axes the Gene map uses — so all three of the email's questions sit on one row."),
        plotOutput("mi_candidates", height = "480px"),
        dl_data_ui("mi_spec_tab"), DTOutput("mi_spec_tab"))))),

  nav_panel("Gene-set Venn", layout_sidebar(
    sidebar = sidebar(width = 330,
      selectInput("vn_a", "Set A", choices = vn_choices(), selected = "de:WT_P0_vs_P7:both"),
      selectInput("vn_b", "Set B", choices = vn_choices(), selected = "de:KO_P0_vs_P7:both"),
      selectInput("vn_c", "Set C (optional)",
                  choices = c(list("(none)" = "none"), vn_choices()),
                  selected = "cur:__canonical__"),
      hr(),
      selectInput("vn_cluster", "Subcluster", choices = NULL),
      radioButtons("vn_stratum", "Cells used",
                   c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"), selected = "G1"),
      radioButtons("vn_grid", "DE matrix", choices = fg_grid_choices(), selected = "de"),
      # Same effect-size choice as the WT-programs tab, for the same reason: presto's
      # log2FC scales with expression level, so a symmetric cut cannot see sparse genes
      # (cell cycle, most transcription factors) at all. AUC is rank-based.
      radioButtons("vn_measure", "Effect size",
                   c("log2 fold change" = "lfc", "AUC (rank-based)" = "auc"), selected = "lfc"),
      conditionalPanel("input.vn_measure != 'auc'",
        sliderInput("vn_lfc", "|log2FC| ≥", 0, 2, 0.25, 0.05)),
      conditionalPanel("input.vn_measure == 'auc'",
        sliderInput("vn_auc", "AUC ≥ (0.50 = no effect-size filter)",
                    0.50, 0.85, 0.60, 0.01)),
      numericInput("vn_padj", "padj <", 0.05, 0.001, 1, 0.01),
      hr(),
      helpText("A Venn hides the three things that decide whether an overlap matters — the ",
               "threshold that built each set, the direction of change, and the overlap you would ",
               "get by chance. Every pairwise overlap here is reported against its chance ",
               "expectation on the ", strong("Overlap statistics"), " tab; read the picture with it, ",
               "not instead of it."),
      uiOutput("vn_caveat"),
      helpText(strong("Descriptive only — n = 1 animal per group.")),
      accordion(open = FALSE, accordion_panel("Figure options",
        figure_controls("vn", palette = FALSE, rename = FALSE)))),
    navset_card_tab(id = "vn_tabs",
      nav_panel("Venn", value = "venn",
        uiOutput("vn_note"),
        plotOutput("vn_plot", height = "560px")),
      nav_panel("Overlap statistics", value = "stats",
        helpText("Observed overlap against what independence would predict. ",
                 "Fold near 1 with a non-significant p means the sets are unrelated — ",
                 "which is a real result, not a failed analysis."),
        dl_data_ui("vn_stats"), DTOutput("vn_stats")),
      nav_panel("Region genes", value = "regions",
        helpText("Every disjoint region of the diagram, largest first. ",
                 "Region A is “in set A only”, AB is “in A and B but not C”, and so on."),
        dl_data_ui("vn_regions"), DTOutput("vn_regions"))))),

  nav_panel("WT programs ∩ KO clusters", layout_sidebar(
    sidebar = sidebar(width = 340,
      tags$b("1. WT P0 → P7"),
      selectInput("xc_wt_cluster", "Measured in", choices = NULL),
      selectInput("xc_mat_set", "Maturation category",
                  choices = xc_cat_choices(), selected = "CM maturation"),
      selectInput("xc_cyc_set", "Cell-cycle category",
                  choices = xc_cat_choices(), selected = XC_CANON),
      hr(), tags$b("2. P7 KO vs WT"),
      selectizeInput("xc_mat_clusters", "Maturation clusters", choices = NULL, multiple = TRUE),
      selectizeInput("xc_cyc_clusters", "Cycling clusters", choices = NULL, multiple = TRUE),
      numericInput("xc_minc", "DE in ≥ N of them", 1, 1, 8, 1),
      hr(), tags$b("Shared settings"),
      radioButtons("xc_stratum", "Cells used",
                   c("All cells (raw)" = "all", "G1 only (phase-matched)" = "G1"),
                   selected = "all"),
      radioButtons("xc_grid", "DE matrix", choices = fg_grid_choices(), selected = "de"),
      radioButtons("xc_measure", "Effect size",
                   c("AUC (rank-based)" = "auc", "log2 fold change" = "lfc"), selected = "auc"),
      # Sliders, not boxes: on this data the answer moves a lot across the plausible
      # range (canonical cell-cycle genes up at P7 go 0 -> 3 -> 14 -> 32 as the AUC cut
      # drops 0.65 -> 0.60 -> 0.55 -> 0.50), so the cut is something to sweep, not set once.
      conditionalPanel("input.xc_measure == 'auc'",
        sliderInput("xc_auc", "AUC ≥ (0.50 = no effect-size filter)",
                    0.50, 0.85, 0.60, 0.01)),
      conditionalPanel("input.xc_measure != 'auc'",
        sliderInput("xc_lfc", "|log2FC| ≥", 0, 2, 0.25, 0.05)),
      numericInput("xc_padj", "padj <", 0.05, 0.001, 1, 0.01),
      checkboxInput("xc_hidemt", "Hide mitochondrial (mt-) genes", TRUE),
      hr(),
      div(downloadButton("xc_book", "Download all four comparisons (XLSX)",
                         class = "btn-sm btn-outline-primary"), style = "margin-bottom:4px"),
      helpText(tags$small("One sheet per comparison plus the two step tables, ",
                          "with a README carrying the caveats.")),
      uiOutput("xc_caveat"),
      helpText(strong("Descriptive only — n = 1 animal per group.")),
      accordion(open = FALSE, accordion_panel("Figure options",
        figure_controls("xc", palette = FALSE, rename = FALSE)))),
    navset_card_tab(id = "xc_tabs",
      nav_panel("Step 1 — WT P0→P7", value = "wt",
        helpText("Which genes in each curated category change between WT P0 and WT P7. ",
                 strong("log2FC > 0 means up at P7."), " These categories are small on ",
                 "purpose — a few dozen canonical genes — so single-digit lists are the ",
                 "honest answer, not a failed query."),
        uiOutput("xc_wt_note"),
        dl_fig_ui("xcwt"), plotOutput("xc_wt_plot", height = "300px"),
        dl_data_ui("xc_wt_tab"), DTOutput("xc_wt_tab")),
      nav_panel("Step 2 — P7 KO vs WT", value = "ko",
        helpText("Every named subcluster, both directions. ",
                 strong("log2FC > 0 means up in KO."), " The pivot below the table is the ",
                 "same data as a gene × cluster grid, so “which genes in each” is one glance."),
        uiOutput("xc_ko_note"),
        dl_data_ui("xc_ko_tab"), DTOutput("xc_ko_tab"),
        div(class = "mt-4", h5("Gene × cluster (log2FC, KO/WT)"),
            dl_data_ui("xc_ko_pivot"), DTOutput("xc_ko_pivot"))),
      nav_panel("The four comparisons", value = "venn",
        helpText("Each WT category list crossed against the ", strong("opposite"),
                 " category's cluster group, as asked. Circle areas are not to scale — ",
                 "the numbers carry the counts. Read them with the ",
                 strong("Overlap statistics"), " tab: an overlap only means something ",
                 "against what independence would have predicted."),
        uiOutput("xc_venn_note"),
        fluidRow(
          column(6, dl_fig_ui("xcvenn1"), plotOutput("xc_venn1", height = "380px")),
          column(6, dl_fig_ui("xcvenn2"), plotOutput("xc_venn2", height = "380px"))),
        fluidRow(
          column(6, dl_fig_ui("xcvenn3"), plotOutput("xc_venn3", height = "380px")),
          column(6, dl_fig_ui("xcvenn4"), plotOutput("xc_venn4", height = "380px"))),
        uiOutput("xc_venn_method")),
      nav_panel("Overlap statistics", value = "stats",
        helpText("Observed overlap against the hypergeometric expectation over the shared ",
                 "testable gene space. Fold near 1 with a non-significant p means the two ",
                 "sets are unrelated — a real result, not a failed analysis. ",
                 "No multiple-testing correction across the four rows."),
        dl_data_ui("xc_stats"), DTOutput("xc_stats")),
      nav_panel("Shared genes", value = "genes",
        helpText("One row per gene per comparison, with the WT and KO evidence side by side. ",
                 "Region is “shared”, “WT only” or “KO only”. The ", strong("lfc_CM…"),
                 " columns give the per-cluster KO log2FC so a union never hides which ",
                 "cluster carried a gene."),
        selectInput("xc_gene_cmp", "Comparison", choices = NULL, width = "520px"),
        dl_data_ui("xc_genes"), DTOutput("xc_genes"))))),

  nav_panel("Cell-cycle exit & ploidy", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("cyc_ct", "Cell type", choices = CELLTYPE_CHOICES, selected = CM_DEFAULT_CT),
      hr(), helpText("E2f7/8 govern cardiomyocyte cell-cycle exit and polyploidization. ",
                     "Proliferation vs cytokinesis separates true division from ",
                     "karyokinesis-without-cytokinesis (binucleation / endoreduplication). ",
                     strong("Descriptive only — n = 1, sex-confounded.")),
      accordion(open = FALSE, accordion_panel("Figure options",
        figure_controls("cyc", palette = TRUE, rename = FALSE)))),
    navset_card_tab(
      nav_panel("Score distributions",
        helpText("Proliferation, cytokinesis, cell-cycle-exit and polyploidization-proxy scores ",
                 "across genotype (WT vs KO), faceted by timepoint."),
        plotOutput("cyc_violins", height = "600px"),
        uiOutput("cyc_score_def")),
      nav_panel("Cycling vs cytokinesis", dl_fig_ui("cycsc"),
        plotOutput("cyc_scatter", height = "600px"))))),

  nav_panel("Maturation & metabolism", layout_sidebar(
    sidebar = sidebar(width = 320,
      conditionalPanel("input.matt != 'genemap'",
        selectInput("mat_ct", "Cell type", choices = CELLTYPE_CHOICES, selected = CM_DEFAULT_CT),
        radioButtons("mat_stratum", "Cells used",
                     c("All cells" = "all", "G1 only (phase-matched)" = "G1"), selected = "all")),
      conditionalPanel("input.matt == 'violin'",
        selectInput("mat_score", "Score",
                    c("CM maturation (mature−immature)" = "sig_maturation",
                      "CM maturation, cycle-free" = "sig_maturation_nocc",
                      "Metabolic maturation (FAO−glyc)" = "sig_metabolic",
                      "Mature-CM program" = "sig_mat_mature", "Immature-CM program" = "sig_mat_immature",
                      "Glycolysis" = "sig_glycolysis", "Fatty-acid oxidation" = "sig_faox"))),
      conditionalPanel("input.matt == 'cells'",
        checkboxInput("mat_showcells", "Show individual cells under the contours", TRUE)),
      conditionalPanel("input.matt == 'genemap'",
        radioButtons("gm_panel", "Timepoint",
                     setNames(GM_PANELS, ifelse(GM_PANELS == "avg", "Averaged (P0 & P7)", GM_PANELS)),
                     selected = "avg"),
        sliderInput("gm_dist", "Minimum distance from centre", 0, 0.25, 0, 0.005),
        checkboxGroupInput("gm_quad", "Quadrants", setNames(GM_QUADS, GM_QUADS), selected = GM_QUADS),
        checkboxInput("gm_conf", "Confidently labelled only", FALSE),
        checkboxInput("gm_hidesets", "Hide genes from the scoring sets", TRUE),
        selectInput("gm_geneset", "Restrict to gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("gm_gene", "Find a gene", choices = NULL, options = list(maxOptions = 50L)),
        numericInput("gm_labeln", "Label top N by distance", 20, 0, 200, 5)),
      hr(),
      conditionalPanel("input.matt == 'genemap'",
        helpText("Each point is a ", strong("gene"), ". x = how strongly it marks mature vs ",
                 "immature cardiomyocytes; y = oxidative vs glycolytic. Both are AUCs from a ",
                 "tertile split run within each timepoint and averaged, so neither axis is a ",
                 "restatement of P0-vs-P7.", br(), br(),
                 strong("Distance"), " is how far a gene sits from the centre — the ranking of ",
                 "how strongly it defines the joint program.", br(), br(),
                 strong("Timepoint"), ": the averaged panel is what keeps the axis a maturation ",
                 "axis rather than a P0-vs-P7 axis. The P0 and P7 panels show the same genes scored ",
                 "within one age only — useful for asking whether a gene is maturation-linked at one ",
                 "age and not the other, but they are not equally powered (P0 splits ~1,144 vs 1,144 ",
                 "cells, P7 only ~775 vs 775), so treat a difference between them cautiously.", br(), br(),
                 "Worth knowing before comparing panels: the two axes are ", strong("tightly coupled at P0 ",
                 "and decoupled at P7"), " (correlation +0.68 vs +0.13 over all genes; +0.89 vs +0.16 ",
                 "restricted to genes measured well at that age). So the P7 panel legitimately looks ",
                 "more scattered. Also, a gene near the centre flips side on jitter alone — the two ",
                 "panels' centres differ by about as much as a typical gene moves — so only movement ",
                 "well clear of the centre is flagged.", br(), br(),
                 "Axes split at each one's ", strong("median"), ", not 0.5, because AUC carries a ",
                 "small global offset; splitting at 0.5 would put most genes in one corner.", br(), br(),
                 strong("Scoring-set genes are hidden by default"), " — they sit at the extremes of ",
                 "their own axis by construction, so leaving them in would partly just recover the inputs.", br(), br(),
                 strong("Ringed points are confidently labelled"), " — clear of the centre by 0.10 and ",
                 "significant, judged within the panel shown. That is the same bar the ",
                 em("Maturation ∩ P7 KO"), " tab uses before it will call a gene maturation-linked, and ",
                 "very few genes clear it. Everything else still gets a side so it can be ranked by ",
                 "distance, but treat an unringed point as a position, not a claim.")),
      conditionalPanel("input.matt != 'genemap'",
        helpText("P0→P7 maturation and the glycolysis→fatty-acid-oxidation metabolic switch. ",
                 "The within-genotype P0→P7 axis is the most robust signal in this design. ",
                 strong("KO-vs-WT descriptive only — n = 1."), br(), br(),
                 "P7 was cycling-enriched 4.5–5.2× relative to P0, so any score that tracks cycling ",
                 "is shifted by the sort. Use the G1 stratum to hold that composition fixed.")),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("Violin figure options", figure_controls("mat", palette = TRUE, rename = FALSE)),
        accordion_panel("Cell-scatter figure options", figure_controls("matsc", palette = FALSE, rename = FALSE)))),
    navset_card_tab(id = "matt",
      nav_panel("By four groups", value = "violin",
        plotOutput("mat_violin", height = "560px"),
        uiOutput("mat_score_def"),
        uiOutput("mat_violin_method")),
      nav_panel("Maturation vs metabolism (cells)", value = "cells",
        plotOutput("mat_scatter", height = "600px"),
        uiOutput("mat_scatter_method")),
      nav_panel("Gene map", value = "genemap",
        helpText("Hover a point for the gene and its coordinates; click a point — or a table row — ",
                 "to ring it and show its details."),
        uiOutput("gm_note"),
        dl_fig_ui("gmsc", "Download figure (static)"),
        plotlyOutput("gm_scatter", height = "560px"),
        uiOutput("gm_pick_ui"),
        dl_data_ui("gm_table"), DTOutput("gm_table"),
        uiOutput("gm_geneinfo"),
        uiOutput("gm_method")))))),

  nav_menu("Methods & provenance",
  nav_panel("Precomputed results", layout_sidebar(
    sidebar = sidebar(width = 340,
      helpText(strong("Computed upstream, linked here."), br(),
               "Results the analysis pipeline already produced. Reading them costs the app",
               "nothing, and several are stronger than what the app computes on the fly."),
      selectInput("lk_group", "Category", choices = NULL),
      selectInput("lk_table", "Result", choices = NULL),
      hr(),
      dl_data_ui("lk_tab"),
      uiOutput("lk_source")),
    card(card_header(textOutput("lk_label")),
         uiOutput("lk_note"),
         DTOutput("lk_tab"))
  )),

  nav_panel("PC dimensions", layout_sidebar(
    sidebar = sidebar(width = 340,
      helpText(strong("How many PCs should the UMAP use?"), br(),
               "Both production embeddings use dims 1:30 and neither records why.",
               "SCTransform, PCA and Harmony are identical across the three panels —",
               "only the number of components carried into the neighbour graph and",
               "UMAP changes."),
      selectInput("pcd_obj", "Object", choices = NULL),
      selectInput("pcd_col", "Colour by", choices = NULL),
      hr(),
      dl_fig_ui("pcdmap", "Download figure (static)"),
      dl_data_ui("pcd_tab"),
      uiOutput("pcd_var")),
    card(card_header(textOutput("pcd_label")),
         uiOutput("pcd_verdict"),
         plotOutput("pcd_map", height = "430px"),
         div(style = "margin-top:10px", uiOutput("pcd_note")),
         DTOutput("pcd_tab"))
  )),

  nav_panel("Gene sets & sources", layout_sidebar(
    sidebar = sidebar(width = 340,
      helpText(strong("Where each gene list came from."), br(),
               "Generated from the code itself, not retyped, so it cannot drift out of",
               "step with the sets actually used."),
      selectInput("gsp_type", "Source", choices = NULL),
      textInput("gsp_find", "Find a gene", placeholder = "e.g. Myh6"),
      hr(), dl_data_ui("gsp_tab"), uiOutput("gsp_counts")),
    card(card_header("Gene sets and their provenance"),
         uiOutput("gsp_headline"),
         uiOutput("gsp_drift"),
         navset_card_tab(id = "gsp_tabs",
           nav_panel("Registry", value = "reg",
             DTOutput("gsp_tab"),
             div(style = "margin-top:10px", uiOutput("gsp_caveats"))),
           nav_panel("Benchmark vs published sets", value = "bench",
             uiOutput("gsp_bench_note"),
             DTOutput("gsp_bench")),
           nav_panel("References", value = "refs",
             uiOutput("gsp_refs"))))))),

  nav_spacer(),
  nav_menu("Help",
  nav_panel("QC & normalization", div(style = "max-width:1000px;padding:8px 4px",
    uiOutput("qcfigs"),
    h5("Doublet rate by lane (numbers)"),
    dl_data_ui("doublet_tab"), div(style = "overflow:auto", tableOutput("doublet_tab")),
    h5("Module-score definitions and coverage", class = "mt-4"),
    helpText("What went into each sig_* score: which curated gene sets, which matrix it was ",
             "scored on, how many of the set's genes were actually present, and how many cells ",
             "were scored. Written by build_signature_scores.R."),
    dl_data_ui("score_meta_tab"), DTOutput("score_meta_tab", height = "320px"),
    h5("The curated gene sets themselves", class = "mt-4"),
    helpText("Every list behind those scores, in full. These are hand-curated canonical ",
             "markers — not MSigDB, not GO — and there is no recorded source for the specific ",
             "choices; see the README section “Module scores” before citing them. Read out of ",
             "build_signature_scores.R at startup rather than copied, so this table cannot ",
             "drift from what was actually scored."),
    dl_data_ui("score_sets_tab"), DTOutput("score_sets_tab", height = "320px"))),

  nav_panel("Annotation check", div(style = "max-width:1000px;padding:8px 4px",
    helpText("Each cell is scored against published-style developmental mouse-heart marker panels; ",
             "the argmax lineage (\"predicted\") is cross-tabulated against the existing cell-type label. ",
             "The diagonal should dominate — off-diagonal mass flags populations worth double-checking. ",
             strong("Concordance check, not probabilistic label transfer.")),
    dl_fig_ui("annheat"), plotlyOutput("ann_heat", height = "520px"),
    h5("Confusion (row-normalised)", style = "margin-top:12px"),
    dl_data_ui("ann_tab"), DTOutput("ann_tab", height = "360px"))),

  nav_panel("About / caveats", div(style = "max-width:820px;padding:8px 4px", htmlOutput("about"))))
)


  # ---- "how this plot was made" blocks -------------------------------------------
# Free functions, not renderUI bodies: each takes the selection as an argument so a test
# can call it with a panel or a score directly. testServer snapshots output$ values, so a
# note written inline in renderUI is effectively untestable -- it would look correct
# forever even after it stopped following its dropdown.
mat_violin_method_note <- function(sc) {
    r  <- if (!is.null(SCOREMETA)) SCOREMETA[SCOREMETA$score == sc, , drop = FALSE] else NULL
    method_note(
      tags$p(HTML(sprintf(paste0("Each violin is <b>one point per cell</b>, grouped genotype ",
        "&times; timepoint. The y value is the per-cell module score <code>%s</code>%s."), sc,
        if (!is.null(r) && nrow(r)) sprintf(paste0(", built from the curated set(s) <b>%s</b>",
          " (%d of %d genes found, on the %s matrix; %s cells scored)"),
          r$sets[1], r$n_genes_used[1], r$n_genes_set[1], r$matrix[1],
          format(r$n_cells_scored[1], big.mark = ",")) else ""))),
      tags$ul(
        tags$li(HTML(paste0("<b>The score is a hand-curated gene list, not a database term.</b> ",
          "It is an <code>AddModuleScore</code> equivalent: mean log-normalised expression of ",
          "the set minus the mean of 100 control genes drawn from the same expression bin ",
          "(24 quantile bins, <code>seed = 1</code>). Raw values &mdash; no z-scoring, no ",
          "rescaling. See the README section &ldquo;Module scores&rdquo; for every list."))),
        tags$li(HTML(paste0("<b>Composites are differences, not ratios.</b> ",
          "<code>sig_maturation</code> = mature &minus; immature; <code>sig_metabolic</code> = ",
          "FAO &minus; glycolysis, so positive means more oxidative."))),
        tags$li(HTML(paste0("<b>On the plot:</b> violins use <code>scale = &quot;width&quot;</code>, ",
          "so every group is drawn the same width and shape is comparable but area is <i>not</i> ",
          "proportional to n. The box is the IQR, the white diamond is the <b>mean</b>, and the ",
          "<code>n=</code> label under each group is its cell count."))),
        tags$li(HTML(paste0("<b>Caveat.</b> P7 was FACS cycling-enriched 4.5&ndash;5.2&times; ",
          "relative to P0, so any score that tracks the cell cycle is shifted by the sort. The ",
          "G1 stratum in the sidebar holds that composition fixed.")))),
      code = c("score_violin() &mdash; shiny_app/app.R",
               "shiny_app/build_signature_scores.R"))
}

mat_scatter_method_note <- function() {
    method_note(
      tags$p(HTML(paste0("Each point is <b>one cell</b>, placed by two module scores: ",
        "x = <code>sig_maturation</code> (mature &minus; immature program), ",
        "y = <code>sig_metabolic</code> (FAO &minus; glycolysis, so up = more oxidative)."))),
      tags$ul(
        tags$li(HTML(paste0("<b>Three layers, and only two of them use all the cells.</b> The ",
          "faint dots are a random thin to 6,000 cells, for texture only &mdash; 30k points is ",
          "an unreadable blob. The contours (<code>stat_density_2d</code>, 5 bins) and the ",
          "centroids are computed on <b>every</b> cell in the group. Read the contours, not the ",
          "dot density."))),
        tags$li(HTML(paste0("<b>Large ringed circles are centroids</b> &mdash; the plain ",
          "arithmetic mean of x and y over that group's cells. The grey arrow is the WT&rarr;KO ",
          "displacement within each timepoint, which is the comparison the tab exists for."))),
        tags$li(HTML(paste0("<b>The subtitle</b> gives that arrow's length (Euclidean distance ",
          "between the two centroids) and, in brackets, its maturation component alone. When the ",
          "two are nearly equal the shift is almost entirely along the maturation axis, with ",
          "little metabolic movement."))),
        tags$li(HTML(paste0("<b>The quadrant percentages</b> in the caption are the share of each ",
          "group's cells with <i>both</i> scores &gt; 0. Note this splits at <b>0</b>, whereas the ",
          "<i>Gene map</i> tab splits at the median AUC &mdash; the same words describe different ",
          "geometry on the two tabs, so those percentages are not comparable across them."))),
        tags$li(HTML(paste0("<b>Caveat.</b> The panels share axes, and P7 was FACS ",
          "cycling-enriched 4.5&ndash;5.2&times; relative to P0. Use the G1 stratum to hold ",
          "cycling composition fixed.")))),
      code = c("mat_scatter() &mdash; shiny_app/app.R",
               "shiny_app/build_signature_scores.R"))
}

gm_method_note <- function(panel) {
    ctr <- gm_centre(panel)
    nc  <- if (!is.null(GM_NCELL) && !is.null(GM_NCELL$mat))
             paste(sprintf("%s: %s", names(GM_NCELL$mat),
                           format(GM_NCELL$mat, big.mark = ",")), collapse = ", ") else NULL
    method_note(
      tags$p(HTML(sprintf(paste0("Each point is a <b>gene</b>, not a cell &mdash; this inverts ",
        "the cell scatter. Coordinates are precomputed into ",
        "<code>app$fourgroup$geneaxes</code> (%s genes) and only filtered here."),
        format(if (is.null(GM)) 0L else nrow(GM), big.mark = ",")))),
      tags$ul(
        tags$li(HTML(paste0("<b>How a coordinate is made.</b> Within <i>each timepoint ",
          "separately</i>, CM cells are split into tertiles of the per-cell score and ",
          "<code>presto::wilcoxauc</code> compares the top third against the bottom third; the ",
          "gene's AUC is its coordinate. The two timepoints are then averaged. Splitting within ",
          "timepoint and averaging afterwards is what stops the axis becoming a P0-vs-P7 axis, ",
          "which the FACS sort confounds.", if (is.null(nc)) "" else
            sprintf(" Tertile split used %s cells.", nc)))),
        tags$li(HTML(paste0("<b>x</b> = maturation AUC, from <code>sig_maturation_nocc</code> ",
          "(the cycle-free variant, so the axis is not partly a cell-cycle score). ",
          "<b>y</b> = metabolic AUC, from <code>sig_metabolic</code>."))),
        tags$li(HTML(sprintf(paste0("<b>The crosshair sits at each axis's median, not at 0.5</b> ",
          "&mdash; for this panel, maturation <b>%.3f</b> and metabolic <b>%.3f</b>. ",
          "<code>wilcoxauc</code>'s AUC carries a small global offset because the two tertile ",
          "groups differ in detection rate, and most genes sit within ~0.02 of the median, so a ",
          "hard 0.5 split put 65%% of them in one corner. Quadrant and distance are recomputed ",
          "against whichever panel you select."), ctr[["mat"]], ctr[["met"]]))),
        tags$li(HTML(sprintf(paste0("<b>Ringed points are &ldquo;confidently labelled&rdquo;</b>: ",
          "at least %.2f clear of this panel's centre <i>and</i> padj &lt; 0.05 on that axis. Very ",
          "few genes clear it, which is the point. Everything else still gets a side so it can be ",
          "ranked, but treat an unringed point as a position, not a claim."), GM_CONF_MARGIN))),
        tags$li(HTML(paste0("<b>Labels</b> name the top N by <code>distance</code> from the centre. ",
          "That ranking <i>is</i> the &ldquo;how strongly does this gene define the joint ",
          "program&rdquo; question."))),
        tags$li(HTML(paste0("<b>Scoring-set genes are hidden by default.</b> The 39 genes that ",
          "define the two scores sit at their own axis's extreme by construction, so leaving them ",
          "in would partly just recover the inputs."))),
        tags$li(HTML(paste0("<b>Two things this map does not filter, and the DE tabs do.</b> ",
          "Sex/construct confounders are still here &mdash; <code>Xist</code> is currently one of ",
          "the confidently-labelled genes, which is the sex difference between the two animals ",
          "being read as maturation. And 36 mitochondrially-encoded genes are on the map, 8 of ",
          "them confidently labelled, reflecting the library read-fraction difference between the ",
          "samples rather than metabolism. Treat both as artifacts of the design.")))),
      code = c("gm_df() / gm_plot_ly() &mdash; shiny_app/app.R",
               "axis_association() &mdash; shiny_app/build_fourgroup.R"))
}

xc_venn_method_note <- function(p) {
    method_note(
      tags$p(HTML(paste0("Each diagram crosses <b>two gene lists</b>, built from the same ",
        "four-group DE grid the other cardiomyocyte tabs read (<code>app$fourgroup$de</code>)."))),
      tags$ul(
        tags$li(HTML(sprintf(paste0("<b>Set A (blue)</b> is WT P0&rarr;P7 in <b>%s</b>, filtered ",
          "to one curated gene category. <b>Set B (red)</b> is P7 KO-vs-WT, unioned across the ",
          "named subclusters &mdash; every DE gene there, <i>not</i> filtered by category. ",
          "Filtering both sides by category would intersect the maturation set with the ",
          "cell-cycle set, and those two share only Mki67 and Top2a, so the first two diagrams ",
          "would read empty by construction rather than by biology."), p$wt_cluster))),
        tags$li(HTML(sprintf(paste0("<b>Membership rule.</b> padj &lt; %s and %s, in the direction ",
          "named in the title. A gene joins set B if it passes in at least %d of the named ",
          "subclusters."), p$padj, xc_eff_label(p), p$minc))),
        tags$li(HTML(paste0("<b>Why AUC is the default effect size.</b> <code>presto</code>'s ",
          "log2FC is a difference of mean log-normalised expression, so it scales with expression ",
          "level: from WT P0 to P7 <i>Mcm3</i> quadruples its detection rate (6.6% &rarr; 27.2%) ",
          "for a log2FC of 0.16, while <i>Myh7</i> moves 2.5. A symmetric |log2FC| cut classifies ",
          "maturation genes and can never classify a sparse cell-cycle one. AUC is rank-based and ",
          "means the same thing on both."))),
        tags$li(HTML(paste0("<b>Circle areas are not to scale</b> and the two sets are lopsided ",
          "by design &mdash; one curated category (a few genes) against a whole cluster group ",
          "(hundreds). The numbers carry the counts, and the <i>Overlap statistics</i> tab carries ",
          "the null: an overlap only means something against what independence would predict."))),
        tags$li(HTML(paste0("<b>Geometry.</b> Two equal-radius circles drawn with ",
          "<code>geom_polygon</code> over sampled circle points &mdash; no Venn package. Region ",
          "counts come from set membership, not from the drawing.")))),
      code = c("xc_comparison() / xc_wt_set() / xc_ko_set() &mdash; shiny_app/app.R",
               "vn_plot() / vn_stats() &mdash; shiny_app/app.R",
               "shiny_app/build_fourgroup.R"))
}

# -------------------------------------------------------------- SERVER --------
server <- function(input, output, session) {
  # Figure Studio: when the browser blocks the window.open from studio_js(),
  # offer the same URL as a plain link the user clicks themselves.
  if (STUDIO_ON) observeEvent(input$studio_popup_blocked, {
    showModal(modalDialog(title = "Figure Studio", easyClose = TRUE,
      p("Your browser blocked the pop-up. Open the figure here instead:"),
      tags$a(href = input$studio_popup_blocked, target = "_blank",
             class = "btn btn-primary", "Open in Figure Studio"),
      footer = modalButton("Close")))
  })

  for (id in c("gene","g2","cm_gene")) {
    sel <- if ("Gabbr2" %in% ALL_GENES) "Gabbr2" else ALL_GENES[1]
    updateSelectizeInput(session, id, choices = ALL_GENES, selected = sel, server = TRUE)
  }
  # gene-set quick-pick: each "Gene set" selector repopulates its paired gene input,
  # keeping the current gene if it's still in the chosen set.
  .geneset_pairs <- list(geneset = "gene", geneset2 = "g2", cm_geneset = "cm_gene")
  for (sid in names(.geneset_pairs)) local({
    set_id <- sid; gene_id <- .geneset_pairs[[sid]]
    observeEvent(input[[set_id]], {
      gl  <- genes_for_set(input[[set_id]])
      cur <- input[[gene_id]]
      sel <- if (!is.null(cur) && nzchar(cur) && cur %in% gl) cur else gl[1]
      updateSelectizeInput(session, gene_id, choices = gl, selected = sel, server = TRUE)
    }, ignoreInit = TRUE)
  })

  # ---- publication figure controls: per-plot label-rename stores ----
  # wires an editable rename table for `prefix`; `levels_fn` yields the current
  # category levels, `reset_on` (a reactive) clears the map when the variable changes.
  make_rename <- function(prefix, levels_fn, reset_on = NULL) {
    rv <- reactiveVal(list())
    output[[paste0(prefix, "_renametab")]] <- renderDT(rename_table(levels_fn(), rv()))
    observeEvent(input[[paste0(prefix, "_renametab_cell_edit")]], {
      e <- input[[paste0(prefix, "_renametab_cell_edit")]]; lv <- levels_fn(); i <- e$row + 1
      if (length(lv) && i >= 1 && i <= length(lv)) { m <- rv(); m[[ lv[i] ]] <- e$value; rv(m) }
    })
    if (!is.null(reset_on)) observeEvent(reset_on(), rv(list()), ignoreInit = TRUE)
    rv
  }
  umap_rn    <- make_rename("umap", function() if (input$color_by %in% CAT_COLS) levels(factor(meta[[input$color_by]])) else character(0), reactive(input$color_by))
  vln_rn     <- make_rename("vln",  function() levels(factor(meta[[input$grp]])), reactive(input$grp))
  dot_rn     <- make_rename("dot",  function() levels(factor(meta[[input$grp]])), reactive(input$grp))
  comp_rn    <- make_rename("comp", function() levels(factor(meta[[input$comp_fill]])), reactive(input$comp_fill))
  cmphase_rn <- make_rename("cmphase", function() c("G1","S","G2M"))
  e2f_rn      <- make_rename("e2f",      function() levels(factor(meta$genotype)))

  # cell-type choices depend on timepoint
  observeEvent(input$ct_tp, {
    if (!length(input$ct_tp) || !nzchar(input$ct_tp)) return()
    cts <- sub(paste0("^", input$ct_tp, "_"), "", grep(paste0("^", input$ct_tp, "_"), names(ctDE), value = TRUE))
    updateSelectInput(session, "ct_sel", choices = setNames(cts, gsub("_", " ", cts)),
                      selected = if ("Cardiomyocyte" %in% cts) "Cardiomyocyte" else cts[1])
  }, ignoreNULL = FALSE)
  # populate the DE subcluster picker once (only res 0.2 is exposed)
  local({
    subs <- names(subDE[["res0.2"]])
    subs <- subs[order(as.integer(sub("CM", "", subs)))]
    ch <- setNames(subs, vapply(subs, function(s) sub_label("0.2", s), ""))
    # "All cardiomyocytes" exists only in the four-group grids -- subDE has no pooled-CM
    # entry. Offer it unconditionally and let cm_d() say so if it is picked while the
    # pooled comparison is selected, rather than making the option appear and disappear.
    if ("AllCM" %in% FG_CLUSTERS) ch <- c("All cardiomyocytes" = "AllCM", ch)
    updateSelectInput(session, "cm_sub", choices = ch, selected = subs[1])
  })

  # ---- UMAP explorer ----
  output$umap_title <- renderText({
    if (input$color_by != "gene") return(paste0("Coloured by ", labof(input$color_by)))
    note <- if (!in_panel(input$gene)) "  (broad matrix — shown on the ~8k-cell subset)" else ""
    paste0("Expression of ", input$gene, note)
  })
  output$umap <- renderPlotly({
    cb <- input$color_by; cont <- cb %in% c("gene", CONT_COLS); ps <- input$ptsize
    leg <- input$umap_legend %||% TRUE; pal <- input$umap_palette %||% "Default"; mp <- umap_rn()
    deflab <- if (cb == "gene") input$gene else labof(cb)
    ttl <- if (!is.null(input$umap_title) && nzchar(input$umap_title)) input$umap_title
           else if (cont) deflab else NULL
    p <- if (input$split == "none") {
      if (cont) umap_cont(meta, if (cb == "gene") expr_vec(input$gene, meta$cell) else meta[[cb]],
                          ttl, psize = ps, legend = leg)
      else umap_cat(meta, cb, ttl = ttl, psize = ps, labels = input$umap_labels %||% TRUE,
                    labsize = input$umap_labelsize %||% 12, legend = leg, pal_choice = pal, map = mp)
    } else umap_split(meta, cb, input$split, gene = if (cb == "gene") input$gene else NULL,
                      continuous = cont, psize = max(2, ps - 1), legend = leg, pal_choice = pal, map = mp)
    if (isTRUE(input$umap_export_on))
      config(p, displaylogo = FALSE,
             toImageButtonOptions = list(format = "png", filename = paste0("umap_", Sys.Date()),
               width = input$umap_w %||% 1200, height = input$umap_h %||% 900, scale = input$umap_scale %||% 2))
    else config(p, displaylogo = FALSE)   # default camera download behaviour
  })
  # static ggplot twin of the plotly UMAP (same dispatch), for vector export + studio
  umap_gg_p <- reactive({
    cb <- input$color_by; cont <- cb %in% c("gene", CONT_COLS); ps <- input$ptsize
    leg <- input$umap_legend %||% TRUE; pal <- input$umap_palette %||% "Default"; mp <- umap_rn()
    deflab <- if (cb == "gene") input$gene else labof(cb)
    ttl <- if (!is.null(input$umap_title) && nzchar(input$umap_title)) input$umap_title
           else if (cont) deflab else NULL
    if (input$split == "none")
      umap_gg(meta, cb, gene = if (cb == "gene") input$gene else NULL, continuous = cont,
              ttl = ttl, psize = ps, labels = input$umap_labels %||% TRUE,
              labsize = input$umap_labelsize %||% 12, legend = leg, pal_choice = pal, map = mp)
    else
      umap_gg(meta, cb, splitvar = input$split, gene = if (cb == "gene") input$gene else NULL,
              continuous = cont, ttl = ttl, psize = max(2, ps - 1), legend = leg,
              pal_choice = pal, map = mp)
  })
  register_fig(output, "umapgg", umap_gg_p, input)

  # ---- Gene detail ----
  vln_plot <- reactive({
    bs <- input$vln_basesize %||% 13; ch <- input$vln_palette %||% "Default"
    df <- meta; df$expr <- expr_vec(input$g2, df$cell)
    df$grp <- factor(relab(df[[input$grp]], vln_rn()), levels = relab(levels(factor(df[[input$grp]])), vln_rn()))
    base <- theme_minimal(base_size = bs) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    if (input$sp2 != "none") {
      df$splitv <- factor(df[[input$sp2]])
      p <- ggplot(df, aes(grp, expr, fill = splitv)) +
        geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2, position = position_dodge(.9)) +
        base + labs(x = labof(input$grp), y = paste0(input$g2, " (log-norm)"), fill = labof(input$sp2))
      if (ch != "Default") p <- p + scale_fill_manual(values = disc_pal(levels(df$splitv), ch))
    } else {
      p <- ggplot(df, aes(grp, expr, fill = grp)) + geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2) +
        base + guides(fill = "none") + labs(x = labof(input$grp), y = paste0(input$g2, " (log-norm)"))
      if (ch != "Default") p <- p + scale_fill_manual(values = disc_pal(levels(df$grp), ch))
    }
    p
  })
  output$vln <- renderPlot(apply_fig_opts(vln_plot(), "vln", input))
  register_fig(output, "vln", vln_plot, input)

  dot_plot <- reactive({
    bs <- input$dot_basesize %||% 13
    df <- meta; df$expr <- expr_vec(input$g2, df$cell); sp <- input$sp2
    df$grp <- factor(relab(df[[input$grp]], dot_rn()), levels = relab(levels(factor(df[[input$grp]])), dot_rn()))
    key <- if (sp != "none") interaction(df$grp, df[[sp]], sep = " · ", drop = TRUE) else df$grp
    agg <- do.call(rbind, lapply(split(df, key), function(s) data.frame(
      grp = s$grp[1], split = if (sp != "none") as.character(s[[sp]][1]) else "all",
      pct = 100 * mean(s$expr > 0), mean = mean(s$expr))))
    if (sp != "none") agg$split <- factor(agg$split, levels = levels(factor(df[[sp]])))   # honour WT-before-KO
    ggplot(agg, aes(grp, split, size = pct, color = mean)) + geom_point() +
      scale_color_viridis_c(option = "magma", direction = -1) + scale_size_area(max_size = 12) +
      theme_minimal(base_size = bs) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
      labs(x = labof(input$grp), y = if (sp != "none") labof(sp) else "", size = "% expr", color = "mean", title = input$g2)
  })
  output$dot <- renderPlot(apply_fig_opts(dot_plot(), "dot", input))
  register_fig(output, "dot", dot_plot, input)
  output$g2_cells <- renderUI({
    g <- input$g2; req(g)
    n <- sum(!is.na(expr_vec(g, meta$cell)))
    if (in_panel(g))
      HTML(sprintf("Showing <b>%s</b> across <b>%s</b> cells (curated panel, full resolution).",
                   g, format(n, big.mark = ",")))
    else
      HTML(sprintf("Showing <b>%s</b> across <b>%s</b> cells — broad subset (gene not in the curated panel, so fewer cells are stored).",
                   g, format(n, big.mark = ",")))
  })

  # ---- Composition ----
  comp_plot <- reactive({
    bs <- input$comp_basesize %||% 13; ch <- input$comp_palette %||% "Default"
    df <- meta; x <- input$comp_x; f <- input$comp_fill
    tb <- as.data.frame(prop.table(table(df[[x]], df[[f]]), margin = 1)); names(tb) <- c("x","fill","prop")
    tb$fill <- factor(relab(tb$fill, comp_rn()), levels = relab(levels(factor(df[[f]])), comp_rn()))
    p <- ggplot(tb, aes(x, prop, fill = fill)) + geom_col() + theme_minimal(base_size = bs) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1)) + labs(x = labof(x), y = "proportion", fill = labof(f))
    if (ch != "Default") p <- p + scale_fill_manual(values = disc_pal(levels(tb$fill), ch))
    p
  })
  output$comp <- renderPlot(apply_fig_opts(comp_plot(), "comp", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("comp_dl_", f)]] <- dl_ggplot("comp", comp_plot, input, f) })

  # ---- DE by cell type ----
  ct_d    <- reactive({ req(input$ct_tp, input$ct_sel); ctDE[[paste(input$ct_tp, input$ct_sel, sep = "_")]] })
  ct_pick <- reactiveVal(NULL)                                     # gene selected on the volcano or in the table
  ct_tab  <- reactive(de_table(drop_conf(ct_d(), input$ct_hideconf), input$ct_search))  # full table
  ct_dt_proxy <- DT::dataTableProxy("ct_table")
  output$ct_volcano <- renderPlotly(de_volcano_ly(drop_conf(ct_d(), input$ct_hideconf),
                         paste0(input$ct_tp, " ", gsub("_", " ", input$ct_sel)), "ct_volcano",
                         highlight = ct_pick(), lfc_cut = input$ct_vlfc %||% DE_LFC_CUT))
  outputOptions(output, "ct_volcano", suspendWhenHidden = FALSE)  # render at startup so plotly_click source registers before its click observer fires
  # static twin for vector export + studio (same data, title and colour cut)
  ct_volcano_gg <- reactive(de_volcano(drop_conf(ct_d(), input$ct_hideconf),
    paste0(input$ct_tp, " ", gsub("_", " ", input$ct_sel)),
    lfc_cut = input$ct_vlfc %||% DE_LFC_CUT, highlight = ct_pick()))
  register_fig(output, "ctvolc", ct_volcano_gg, input)
  output$ct_table   <- renderDT(de_datatable(ct_tab()))
  output$ct_pick_ui  <- renderUI(pick_banner(ct_pick(), "ct_clear"))
  output$ct_geneinfo <- renderUI(gene_info_card(ct_pick(), "ct_infoclose"))
  observeEvent(event_data("plotly_click", source = "ct_volcano"),
    ct_pick(event_data("plotly_click", source = "ct_volcano")$customdata))
  observeEvent(input$ct_table_rows_selected, {                          # table row -> selected gene
    r <- input$ct_table_rows_selected; g <- if (length(r)) ct_tab()$gene[r] else NULL
    if (!identical(g, ct_pick())) ct_pick(g)
  }, ignoreNULL = FALSE)
  observeEvent(ct_pick(), {                                             # selected gene -> highlight table row
    g <- ct_pick(); rows <- if (!is.null(g)) which(ct_tab()$gene == g) else integer(0)
    if (!identical(as.integer(rows), as.integer(input$ct_table_rows_selected)))
      DT::selectRows(ct_dt_proxy, if (length(rows)) rows else NULL)
  }, ignoreNULL = FALSE)
  observeEvent(input$ct_clear, ct_pick(NULL))
  observeEvent(input$ct_infoclose, ct_pick(NULL))
  observeEvent(list(input$ct_tp, input$ct_sel, input$ct_search), ct_pick(NULL), ignoreInit = TRUE)
  ct_heat_p <- reactive({
    req(input$ct_tp)
    keys <- grep(paste0("^", input$ct_tp, "_"), names(ctDE), value = TRUE)
    dl <- setNames(ctDE[keys], gsub("_", " ", sub(paste0("^", input$ct_tp, "_"), "", keys)))
    lfc_heat(dl, 24, paste0("Top KO-vs-WT genes across cell types — ", input$ct_tp))
  })
  output$ct_heat <- renderPlotly(ggheat(ct_heat_p()))
  register_fig(output, "ctheat", ct_heat_p, input)

  # ---- CM deep-dive ----
  output$cm_map <- renderPlotly({
    req(input$cm_mapcolor); cb <- input$cm_mapcolor; df <- cmm
    sp <- input$cm_map_split %||% "none"
    if (sp != "none") {                          # split panels — always res 0.2, coloured by subcluster
      df$mapval <- factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
      if (sp == "both") {
        gl <- levels(factor(df$genotype)); tl <- levels(factor(df$timepoint))
        df$splitgrp <- factor(paste(df$genotype, df$timepoint, sep = " · "),
                              levels = as.vector(t(outer(gl, tl, paste, sep = " · "))))
        return(umap_split(df, "mapval", "splitgrp", psize = 4, nrows = 2))
      }
      return(umap_split(df, "mapval", sp, psize = 4))
    }
    if (cb == "gene") return(umap_cont(df, expr_vec(input$cm_gene, df$cell), input$cm_gene, psize = 5))
    df$mapval <- if (cb == "subcluster")
      factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2")) else factor(df[[cb]])
    umap_cat(df, "mapval", ttl = paste0(labof(cb), if (cb == "subcluster") " — res 0.2" else ""),
             psize = 5, labels = (cb == "subcluster"))
  })
  # static ggplot twin of the CM map (same dispatch), for vector export + studio
  cm_map_gg_p <- reactive({
    req(input$cm_mapcolor); cb <- input$cm_mapcolor; df <- cmm
    sp <- input$cm_map_split %||% "none"
    if (sp != "none") {
      df$mapval <- factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
      if (sp == "both") {
        gl <- levels(factor(df$genotype)); tl <- levels(factor(df$timepoint))
        df$splitgrp <- factor(paste(df$genotype, df$timepoint, sep = " · "),
                              levels = as.vector(t(outer(gl, tl, paste, sep = " · "))))
        return(umap_gg(df, "mapval", splitvar = "splitgrp", psize = 4, nrows = 2))
      }
      return(umap_gg(df, "mapval", splitvar = sp, psize = 4))
    }
    if (cb == "gene") return(umap_gg(df, NULL, gene = input$cm_gene, continuous = TRUE,
                                     ttl = input$cm_gene, psize = 5))
    df$mapval <- if (cb == "subcluster")
      factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2")) else factor(df[[cb]])
    umap_gg(df, "mapval", ttl = paste0(labof(cb), if (cb == "subcluster") " — res 0.2" else ""),
            psize = 5, labels = (cb == "subcluster"))
  })
  register_fig(output, "cmmap", cm_map_gg_p, input)

  # ---- immune-contamination note (app$immune_contam, from build_immune_flag.R) ----
  output$cm_icn_note <- renderUI({
    if (is.null(ICN)) return(NULL)
    div(class = "alert alert-warning", style = "font-size:12px",
        HTML(sprintf(paste0(
          "<b>%d of these cells are immune cells, not cardiomyocytes.</b> They carry ",
          "pan-leukocyte markers and were labelled Cardiomyocyte because annotate.R had no ",
          "mast-cell or lymphocyte panel to assign them to. They cycle at <b>%.0f%%</b> against ",
          "%.0f%% for the compartment, so a per-subcluster cycling row that is mostly these ",
          "cells is not a cardiomyocyte number. Excluding them moves the compartment figure ",
          "only %.1f%% &rarr; %.1f%%. Colour the map by <i>Immune contamination</i> to see them."),
          ICN$n_flagged_bundle, ICN$cycling_flagged, ICN$cycling_all,
          ICN$cycling_all, ICN$cycling_excl)))
  })

  # ---- Clustering variants (app$clusterings, from build_clusterings.R) ----
  observe({
    req(CLU)
    updateSelectInput(session, "clu_var", choices = clu_choices(), selected = CLU$production)
  })
  clu_v <- reactive({ clu_ok(); req(input$clu_var)
    v <- CLU$variants[[input$clu_var]]
    validate(need(!is.null(v), "That variant is not in this data build.")); v })
  observeEvent(input$clu_var, {
    v <- CLU$variants[[input$clu_var]]; req(v)
    cl <- sort(unique(as.integer(v$labels)))
    updateSelectInput(session, "clu_cl", choices = paste0("CM", cl))
  }, ignoreNULL = TRUE)

  # Reading a non-production variant and reading production must never look the same in a
  # screenshot, so the banner is loud and states which labelling produced the numbers.
  output$clu_banner <- renderUI({
    v <- clu_v()
    if (isTRUE(v$is_production))
      div(class = "alert alert-success", style = "font-size:13px",
          HTML(sprintf("Showing the <b>production</b> labelling (%s, %d subclusters) &mdash; these are the numbers the rest of the app and the book report.", v$label, v$n_clusters)))
    else
      div(class = "alert alert-danger", style = "font-size:13px",
          HTML(sprintf("<b>Not the published labelling.</b> Showing %s (%d subclusters); production is %s. Every number on this page belongs to the selected variant. Subcluster IDs are not comparable across variants &mdash; CM3 here is not CM3 in production.",
                       v$label, v$n_clusters, CLU$variants[[CLU$production]]$label)))
  })

  clu_lab_cells <- reactive({ v <- clu_v(); v$labels })

  # Live markers. presto::wilcoxauc is the same call the "Subset & DEGs" tab uses; cached
  # on (variant, matrix) so flipping back to a variant already seen is instant.
  clu_mk_df <- reactive({
    v <- clu_v(); req(input$clu_mat)
    M <- if (input$clu_mat == "deg") EXPR else expr
    lab <- clu_lab_cells()
    cells <- intersect(names(lab), colnames(M))
    validate(need(length(cells) > 50, "Too few of this variant's cells are in the selected matrix."))
    g <- paste0("CM", unname(lab[cells]))
    keep <- names(table(g))[table(g) >= 10]
    validate(need(length(keep) >= 2, "Fewer than two subclusters have >= 10 cells here."))
    sel <- g %in% keep
    w <- suppressWarnings(presto::wilcoxauc(M[, cells[sel], drop = FALSE], g[sel]))
    w <- w[w$logFC > 0 & w$padj < 0.05, c("group","feature","auc","logFC","pct_in","pct_out","padj")]
    w <- w[order(w$group, -w$auc), ]
    names(w)[1:2] <- c("subcluster", "gene")
    w
  }) |> bindCache(input$clu_var, input$clu_mat)
  output$clu_mk <- renderDT(enr_dt(clu_mk_df(), scroll = "380px"))
  output$clu_mk_note <- renderUI({
    v <- clu_v(); M <- if ((input$clu_mat %||% "deg") == "deg") EXPR else expr
    cells <- intersect(names(clu_lab_cells()), colnames(M))
    helpText(style = "font-size:12px", sprintf(
      "presto::wilcoxauc, one-vs-rest, computed now on %s cells x %s genes. Positive logFC, padj < 0.05. Subclusters with < 10 cells in this matrix are dropped.",
      format(length(cells), big.mark = ","), format(nrow(M), big.mark = ",")))
  })

  # Composition and phase: cheap enough to recompute, so no precompute for these either.
  clu_meta <- reactive({
    lab <- clu_lab_cells()
    m <- meta[match(names(lab), meta$cell), , drop = FALSE]
    m$sub <- factor(paste0("CM", unname(lab)),
                    levels = paste0("CM", sort(unique(as.integer(lab)))))
    m[!is.na(m$cell), , drop = FALSE]
  })
  output$clu_sum_note <- renderUI({
    v <- clu_v()
    helpText(style = "font-size:12px", sprintf(
      "%s — %d subclusters. Left: genotype composition. Right: cell-cycle phase. Both recomputed from the selected labelling.",
      v$label, v$n_clusters))
  })
  output$clu_comp <- renderPlot({
    d <- clu_meta(); validate(need(nrow(d) && "genotype" %in% names(d), "No genotype column."))
    ggplot(d, aes(sub, fill = genotype)) + geom_bar(position = "fill") +
      theme_bw() + labs(x = NULL, y = "fraction", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  output$clu_phase <- renderPlot({
    d <- clu_meta(); validate(need("Phase" %in% names(d), "No Phase column in this build."))
    ggplot(d[!is.na(d$Phase), ], aes(sub, fill = Phase)) + geom_bar(position = "fill") +
      theme_bw() + labs(x = NULL, y = "fraction", fill = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # Precomputed halves: DESeq2 and clusterProfiler are not in the runtime image, so these
  # are looked up rather than computed, and say so when a variant has not been built.
  clu_de_df <- reactive({
    v <- clu_v(); req(input$clu_cl)
    validate(need(!is.null(v$de), "Pseudobulk DE has not been built for this variant — run cm_subcluster_analyze.R --variant=... then build_clusterings.R."))
    d <- v$de[[input$clu_cl]]
    validate(need(!is.null(d) && nrow(d), "No DE table for that subcluster (it may have been too thin to test)."))
    d
  })
  output$clu_de <- renderDT(enr_dt(clu_de_df(), scroll = "380px"))
  output$clu_de_note <- renderUI({
    v <- clu_v()
    helpText(style = "font-size:12px", HTML(paste0(
      "Pseudobulk DESeq2, apeglm-shrunken, precomputed for <b>", v$label, "</b>. ",
      "n = 1 per group and sex-confounded, so these are <b>descriptive log2 fold changes with no valid p-values</b> ",
      "&mdash; the same caveat the production tables carry.")))
  })
  clu_enr_df <- reactive({
    v <- clu_v()
    # Deliberately NO fallback to app$enrich$sub. That belongs to the shipped labelling in
    # app$cm$meta, which is a different clustering from the registry's cm_dims30_res0.2
    # even though both are "dims 30, res 0.2" (ARI 0.68, 13 subclusters against 11).
    # Falling back would silently pair one labelling's enrichment with another's clusters.
    e <- app$enrich$sub_by_variant[[v$variant_id]]
    validate(need(!is.null(e) && !is.null(e$go),
                  "Enrichment has not been built for this variant — run build_subcluster_enrichment.R --variant=... and redeploy."))
    e$go
  })
  output$clu_enr <- renderDT(enr_dt(clu_enr_df(), scroll = "380px"))
  output$clu_enr_note <- renderUI({
    helpText(style = "font-size:12px",
             "KO-vs-WT GO biological process per subcluster, precomputed with clusterProfiler against the DE tables of this same variant.")
  })
  output$clu_cyc <- renderDT({
    v <- clu_v()
    validate(need(!is.null(v$cellcycle), "No cell-cycle table for this variant."))
    enr_dt(v$cellcycle, scroll = "340px")
  })

  # ---- Object-mode test (app$cmtest, from build_cm_objectmode.R) ----
  # Does the CM subset need to be its own object? Three builds of the same cells, the
  # same pipeline and seeds on each, diffed. The verdict line states the answer; the
  # table carries the numbers behind it.
  observe({
    req(CMTEST)
    updateSelectInput(session, "objtest_res", choices = CMTEST$res,
                      selected = if ("0.2" %in% CMTEST$res) "0.2" else CMTEST$res[1])
  })
  output$objtest_verdict <- renderUI({
    cmtest_ok()
    div(class = paste("alert", if (isTRUE(CMTEST$ab_same)) "alert-success" else "alert-warning"),
        style = "font-size:13px", CMTEST$verdict)
  })
  objtest_map_p <- reactive({
    cmtest_ok(); req(input$objtest_res)
    rc <- paste0("SCT_snn_res.", input$objtest_res)
    validate(need(rc %in% names(CMTEST$percell), "That resolution is not in this data build."))
    objtest_gg(CMTEST$percell, rc)
  })
  output$objtest_map <- renderPlot(objtest_map_p())
  register_fig(output, "objtestmap", objtest_map_p, input)
  output$objtest_note <- renderUI({
    cmtest_ok()
    tagList(
      div(style = "font-size:12px;color:#555",
          tags$b("How each arm was built."),
          tags$ul(style = "margin-bottom:4px", lapply(CMTEST$variants, function(v)
            tags$li(tags$b(unname(CMTEST$labels[v])), " — ", unname(CMTEST$blurb[v]))))),
      helpText(style = "font-size:11px", CMTEST$note))
  })
  objtest_df <- reactive({ cmtest_ok(); CMTEST$metrics })
  output$objtest_tab <- renderDT(enr_dt(objtest_df(), scroll = "360px"))

  # ---- PC dimensions (app$pcdims, from build_pcdims.R) ----
  observe({
    req(PCD)
    updateSelectInput(session, "pcd_obj",
      choices = setNames(PCD$objects, vapply(PCD$objects, function(o) PCD[[o]]$label, "")))
  })
  pcd_g <- reactive({ pcd_ok(); req(input$pcd_obj)
    g <- PCD[[input$pcd_obj]]; validate(need(!is.null(g), "That object is not in this build.")); g })
  observeEvent(input$pcd_obj, {
    g <- PCD[[input$pcd_obj]]; req(g)
    updateSelectInput(session, "pcd_col", choices = g$colby, selected = unname(g$colby)[1])
  }, ignoreNULL = TRUE)
  output$pcd_label <- renderText(paste0(pcd_g()$label, " — PC dimensions carried into the UMAP"))
  # The two numbers that answer the question: does the UMAP shape move, and do the
  # clusters people actually read move with it.
  output$pcd_verdict <- renderUI({
    g <- pcd_g()
    fmt <- function(a, b) sprintf("<b>%s</b>: %.0f%% of neighbours kept, shape m&sup2; %.2f, clusters ARI %.2f",
                                  a, 100 * g$knn30[[b]], g$m2[[b]], g$ari[[b]])
    div(class = "alert alert-warning", style = "font-size:13px",
        HTML(paste0("Against the production cut of 30 &mdash; ",
                    fmt("10 vs 30", "10v30"), "<br>", fmt("30 vs 50", "30v50"),
                    "<br><span style='font-size:12px;color:#555'>",
                    "<b>Neighbours kept</b>: share of each cell's 30 nearest neighbours that are still ",
                    "neighbours under the other cut &mdash; 100% would mean the local structure is identical. ",
                    "<b>m&sup2;</b>: Procrustes shape distance, 0 = identical, 1 = unrelated. ",
                    "<b>ARI</b>: 1 = identical cluster assignments, at ",
                    sub("SCT_snn_res.", "resolution ", g$prod_res), ".</span>")))
  })
  output$pcd_var <- renderUI({
    g <- pcd_g()
    helpText(style = "font-size:11px",
      HTML(paste0("Share of the top-50 PC variance<br>",
                  paste(sprintf("&nbsp;dims 1:%s &mdash; %.1f%%", names(g$varpct), g$varpct),
                        collapse = "<br>"),
                  "<br><br>", format(nrow(g$percell) / length(g$dims), big.mark = ","),
                  " cells per panel")))
  })
  pcd_map_p <- reactive({ g <- pcd_g(); req(input$pcd_col); pcdims_gg(g, input$pcd_col) })
  output$pcd_map <- renderPlot(pcd_map_p())
  register_fig(output, "pcdmap", pcd_map_p, input)
  output$pcd_note <- renderUI({ pcd_ok(); helpText(style = "font-size:11px", PCD$note) })
  pcd_df <- reactive({ pcd_g()$metrics })
  output$pcd_tab <- renderDT(enr_dt(pcd_df(), scroll = "300px"))

  # ---- Gene sets & sources (app$genesets, from build_gene_provenance.R) ----
  observe({
    req(GSP)
    updateSelectInput(session, "gsp_type",
      choices = c("All sources" = "__all__", sort(unique(GSP$registry$source_type))))
  })
  output$gsp_headline <- renderUI({
    gsp_ok()
    div(class = "alert alert-warning", style = "font-size:13px", GSP$headline)
  })
  output$gsp_drift <- renderUI({
    gsp_ok(); d <- GSP$drift
    if (is.null(d) || !any(d$drifted))
      return(div(class = "alert alert-success", style = "font-size:12px",
                 "Duplicate-name check: every panel defined in more than one place holds the same genes."))
    bad <- d[d$drifted, , drop = FALSE]
    div(class = "alert alert-danger", style = "font-size:12px",
        HTML(paste0("<b>", nrow(bad), " panel(s) share a name but not their genes</b> \u2014 ",
                    "the same label means different things in different parts of this app.<ul>",
                    paste0("<li><b>", bad$set_a, "</b> (", bad$n_a, " genes) vs <b>", bad$set_b,
                           "</b> (", bad$n_b, " genes); only in the second: <code>",
                           bad$only_in_b, "</code></li>", collapse = ""), "</ul>")))
  })
  gsp_df <- reactive({
    gsp_ok(); d <- GSP$registry
    if (!is.null(input$gsp_type) && input$gsp_type != "__all__")
      d <- d[d$source_type == input$gsp_type, , drop = FALSE]
    q <- trimws(input$gsp_find %||% "")
    if (nzchar(q))    # word-boundary match so "Ckm" does not also hit "Ckmt2"
      d <- d[grepl(paste0("(^|, )", q, "($|,)"), d$genes, ignore.case = TRUE), , drop = FALSE]
    d
  })
  output$gsp_tab <- renderDT(enr_dt(gsp_df(), scroll = "420px"))
  output$gsp_counts <- renderUI({
    gsp_ok()
    helpText(style = "font-size:11px",
      HTML(sprintf("%d sets shown of %d<br>%d external database<br>%d hand-curated",
                   nrow(gsp_df()), nrow(GSP$registry), GSP$n_ext, GSP$n_hand)))
  })
  output$gsp_caveats <- renderUI({
    gsp_ok()
    div(style = "font-size:12px;color:#555",
        tags$b("Things to know before quoting a score built on these:"),
        tags$ul(lapply(GSP$caveats, tags$li)))
  })
  output$gsp_bench_note <- renderUI({
    gsp_ok()
    validate(need(!is.null(GSP$benchmark),
                  "Benchmark not in this data build — run our_analysis/05_analyses/gene_set_benchmark.R, then build_gene_provenance.R."))
    n_un <- sum(GSP$benchmark$reference == "(none found)")
    tagList(
      div(class = "alert alert-secondary", style = "font-size:12px", GSP$bench_note),
      if (n_un > 0) div(class = "alert alert-danger", style = "font-size:12px",
        HTML(paste0("<b>", n_un, " panel(s) have no published counterpart at all.</b> ",
                    "The CM maturation panels are the important case: the maturation score ",
                    "in this app rests on a list with no external anchor. Uosaki et al. 2015 ",
                    "(References tab) is the atlas to reconcile them against."))))
  })
  gsp_bench_df <- reactive({
    gsp_ok(); validate(need(!is.null(GSP$benchmark), "Benchmark not in this data build."))
    GSP$benchmark
  })
  output$gsp_bench <- renderDT(enr_dt(gsp_bench_df(), scroll = "420px"))
  output$gsp_refs <- renderUI({
    gsp_ok()
    validate(need(!is.null(GSP$refs), "References not in this data build."))
    r <- GSP$refs
    div(style = "font-size:13px",
      p(style = "color:#555",
        "Checked against the actual record, not recalled. Where a panel has no entry here, ",
        "it has no source \u2014 that is the finding, not an omission."),
      tags$ul(lapply(seq_len(nrow(r)), function(i) tags$li(
        style = "margin-bottom:9px",
        tags$b(r$topic[i]), tags$br(),
        tags$a(href = r$link[i], target = "_blank", rel = "noopener", r$reference[i]),
        tags$br(), tags$span(style = "color:#555;font-size:12px", r$relevance[i])))))
  })
  cm_markerheat_p <- reactive({
    h <- heat[["res0.2"]]; validate(need(!is.null(h), "No marker heatmap for this resolution."))
    long <- h$long; long$gene <- factor(long$gene, levels = rev(h$genes)); long$cluster <- factor(long$cluster, levels = h$clusters)
    p <- ggplot(long, aes(cluster, gene, fill = z)) + geom_tile() +
      scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
      theme_minimal(base_size = 11) + theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = "subcluster", y = "marker gene", fill = "z-score\nmean expr",
           title = "Subcluster identity markers — res 0.2")
    p
  })
  output$cm_markerheat <- renderPlotly(ggheat(cm_markerheat_p()))
  register_fig(output, "cmmarker", cm_markerheat_p, input)
  # Which comparison the DE panel is showing. NULL = the pooled KO-vs-WT tables (subDE);
  # anything else is one of build_fourgroup.R's timepoint-specific contrasts, read from
  # the same grids the Four-group tab uses so the two tabs cannot disagree.
  cm_ct <- reactive(if (identical(input$cm_contrast %||% "pooled", "pooled")) NULL
                    else fg_ct(input$cm_contrast))
  cm_d    <- reactive({ req(input$cm_sub); ct <- cm_ct()
    if (is.null(ct)) {
      validate(need(!identical(input$cm_sub, "AllCM"),
        paste("Pooled KO-vs-WT was not computed over all cardiomyocytes together.",
              "Pick a subcluster, or one of the timepoint-specific comparisons.")))
      subDE[["res0.2"]][[input$cm_sub]]
    } else fg_de(input$cm_sub, ct$key, input$cm_stratum %||% "G1", input$cm_grid %||% "de") })
  cm_pick <- reactiveVal(NULL)
  cm_tab  <- reactive(de_table(drop_conf(cm_d(), input$cm_hideconf)))   # full table (selection never filters it)
  cm_dt_proxy <- DT::dataTableProxy("cm_detab")
  # Same caveats the Four-group tab prints, from the same helper.
  output$cm_de_note <- renderUI({ req(input$cm_sub); ct <- cm_ct()
    if (is.null(ct))
      return(div(style = "font-size:13px;margin-bottom:6px",
                 HTML(if (identical(input$cm_sub, "AllCM"))
                        paste0("<b>All cardiomyocytes</b> has no pooled KO-vs-WT table",
                               " &middot; pick a timepoint-specific <i>Comparison</i>, or a subcluster.")
                      else
                        paste0("<b>KO vs WT</b> in ", input$cm_sub,
                               " &middot; P0 and P7 pooled. Use <i>Comparison</i> in the sidebar",
                               " for the timepoint-specific contrasts."))))
    d <- try(cm_d(), silent = TRUE)
    de_context_note(input$cm_sub, ct, input$cm_stratum %||% "G1", input$cm_grid %||% "de",
                    if (inherits(d, "try-error")) NULL else d) })
  output$cm_volcano <- renderPlotly({ ct <- cm_ct()
    ttl <- if (is.null(ct)) paste0(input$cm_sub, " — ", sub_label("0.2", input$cm_sub))
           else paste0(ct$label, " — ",
                       if (identical(input$cm_sub, "AllCM")) "all CM" else input$cm_sub,
                       if (identical(input$cm_stratum, "G1")) " (G1)" else "")
    de_volcano_ly(drop_conf(cm_d(), input$cm_hideconf), ttl, "cm_volcano",
                  pos = ct$pos %||% "up in KO", neg = ct$neg %||% "up in WT",
                  xlab = ct$xlab %||% "log2 fold change (KO / WT)",
                  highlight = cm_pick(), lfc_cut = input$cm_vlfc %||% DE_LFC_CUT) })
  outputOptions(output, "cm_volcano", suspendWhenHidden = FALSE)  # render at startup so plotly_click source registers before its click observer fires
  # static twin for vector export + studio (same contrast, title and colour cut)
  cm_volcano_gg <- reactive({ ct <- cm_ct()
    ttl <- if (is.null(ct)) paste0(input$cm_sub, " — ", sub_label("0.2", input$cm_sub))
           else paste0(ct$label, " — ",
                       if (identical(input$cm_sub, "AllCM")) "all CM" else input$cm_sub,
                       if (identical(input$cm_stratum, "G1")) " (G1)" else "")
    de_volcano(drop_conf(cm_d(), input$cm_hideconf), ttl,
               pos = ct$pos %||% "up in KO", neg = ct$neg %||% "up in WT",
               xlab = ct$xlab %||% "log2 fold change (KO / WT)",
               lfc_cut = input$cm_vlfc %||% DE_LFC_CUT, highlight = cm_pick()) })
  register_fig(output, "cmvolc", cm_volcano_gg, input)
  output$cm_detab   <- renderDT(de_datatable(cm_tab(), scroll = NULL))
  output$cm_pick_ui  <- renderUI(pick_banner(cm_pick(), "cm_clear"))
  output$cm_geneinfo <- renderUI(gene_info_card(cm_pick(), "cm_infoclose"))
  observeEvent(event_data("plotly_click", source = "cm_volcano"),
    cm_pick(event_data("plotly_click", source = "cm_volcano")$customdata))
  observeEvent(input$cm_detab_rows_selected, {                          # table row -> selected gene
    r <- input$cm_detab_rows_selected; g <- if (length(r)) cm_tab()$gene[r] else NULL
    if (!identical(g, cm_pick())) cm_pick(g)
  }, ignoreNULL = FALSE)
  observeEvent(cm_pick(), {                                             # selected gene -> highlight table row
    g <- cm_pick(); rows <- if (!is.null(g)) which(cm_tab()$gene == g) else integer(0)
    if (!identical(as.integer(rows), as.integer(input$cm_detab_rows_selected)))
      DT::selectRows(cm_dt_proxy, if (length(rows)) rows else NULL)
  }, ignoreNULL = FALSE)
  observeEvent(input$cm_clear, cm_pick(NULL))
  observeEvent(input$cm_infoclose, cm_pick(NULL))
  observeEvent(list(input$cm_sub, input$cm_contrast, input$cm_stratum, input$cm_grid),
               cm_pick(NULL), ignoreInit = TRUE)
  # The heatmap answers "this comparison across every subcluster at once", so it follows
  # the same dropdown. lfc_heat() drops absent tables silently, which would read as
  # "these are all the subclusters" -- cm_lfcheat_note names the ones that fell out.
  cm_lfc_list <- reactive({ ct <- cm_ct()
    if (is.null(ct)) return(subDE[["res0.2"]])
    key <- paste0(ct$key, "__", input$cm_stratum %||% "G1")
    g   <- fg_grid(input$cm_grid %||% "de")
    cl  <- setdiff(FG_CLUSTERS, "AllCM")
    setNames(lapply(cl, function(x) g[[x]][[key]]), cl) })
  cm_lfcheat_p <- reactive({ ct <- cm_ct()
    lfc_heat(cm_lfc_list(), 22,
      if (is.null(ct)) "KO-vs-WT log2FC across CM subclusters — res 0.2"
      else paste0(ct$label, " — log2FC across CM subclusters (",
                  if (identical(input$cm_stratum, "G1")) "G1 cells" else "all cells", ")"),
      fill_lab = if (is.null(ct)) "log2FC\n(KO/WT)"
                 else paste0("log2FC\n", sub("^[^(]*", "", ct$xlab))) })
  output$cm_lfcheat <- renderPlotly(ggheat(cm_lfcheat_p()))
  output$cm_lfcheat_note <- renderUI({ ct <- cm_ct(); if (is.null(ct)) return(NULL)
    l <- cm_lfc_list(); gone <- names(l)[vapply(l, is.null, logical(1))]
    if (!length(gone)) return(NULL)
    div(style = "font-size:12px;color:#c62828;margin-top:4px",
        HTML(sprintf("Not shown — no DE table for this comparison and stratum: <b>%s</b>. %s",
                     paste(gone, collapse = ", "),
                     if (identical(input$cm_stratum, "G1"))
                       "Try the “All cells” stratum, or the other DE matrix."
                     else "Try the other DE matrix."))) })
  register_fig(output, "cmlfcheat", cm_lfcheat_p, input)
  cm_phase_plot <- reactive({
    bs <- input$cmphase_basesize %||% 13; ch <- input$cmphase_palette %||% "Default"; rn <- cmphase_rn()
    df <- cmm; df$sub <- factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
    validate(need("Phase" %in% names(df), "No cell-cycle phase data."))
    has_tp <- "timepoint" %in% names(df)
    split_cols <- if (has_tp) c("genotype","timepoint") else "genotype"
    # fraction of each Phase within sub × (split groups): margin = every dim but Phase (dim 2)
    tab <- table(df[c("sub","Phase", split_cols)])
    tb <- as.data.frame(prop.table(tab, setdiff(seq_along(dim(tab)), 2L)))
    names(tb)[match("Freq", names(tb))] <- "prop"
    tb$prop[is.nan(tb$prop)] <- 0     # empty subcluster × group cells -> 0, not NaN
    tb$Phase <- factor(relab(tb$Phase, rn), levels = relab(c("G1","S","G2M"), rn))
    fw <- if (length(split_cols) == 2) facet_grid(reformulate(split_cols[2], split_cols[1])) else facet_wrap(reformulate(split_cols))
    fillv <- if (ch == "Default") setNames(c("#bdbdbd","#1565c0","#c62828"), relab(c("G1","S","G2M"), rn))
             else disc_pal(relab(c("G1","S","G2M"), rn), ch)
    ggplot(tb, aes(sub, prop, fill = Phase)) + geom_col() + fw +
      scale_fill_manual(values = fillv, na.value = "grey90") +
      theme_minimal(base_size = bs) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      labs(x = "subcluster", y = "fraction of cells",
           title = "Cell-cycle phase by subcluster — res 0.2")
  })
  output$cm_phase <- renderPlot(apply_fig_opts(cm_phase_plot(), "cmphase", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("cmphase_dl_", f)]] <- dl_ggplot("cmphase", cm_phase_plot, input, f) })

  # ---- CM composition stacked bars — one plot per breakdown (res 0.2) ----
  cm_bar_one <- function(v, nm) {                                 # single-dimension stacked bar (own legend)
    mode <- input$cm_bar_mode %||% "prop"
    bs <- input$cmbar_basesize %||% 13; ch <- input$cmbar_palette %||% "Default"
    df <- cmm; df$sub <- factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
    validate(need(v %in% names(df), paste0("No ", nm, " data.")))
    fv <- if (v == "cycling") factor(ifelse(as.logical(df$cycling), "cycling", "non-cycling"),
                                     levels = c("non-cycling","cycling")) else factor(df[[v]])
    tab <- table(df$sub, fv)
    tb  <- if (mode == "prop") as.data.frame(prop.table(tab, 1)) else as.data.frame(tab)
    names(tb) <- c("sub","cat","y"); tb$y[is.nan(tb$y)] <- 0
    p <- ggplot(tb, aes(sub, y, fill = cat)) + geom_col() +
      theme_minimal(base_size = bs) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      labs(x = "subcluster", y = if (mode == "prop") "proportion" else "cell count", fill = NULL, title = nm)
    if (ch != "Default") p <- p + scale_fill_manual(values = disc_pal(levels(tb$cat), ch))
    p
  }
  # combined faceted version — used only for the figure download button
  cm_bar_plot <- reactive({
    ps <- list(cm_bar_one("genotype", "By genotype"), cm_bar_one("timepoint", "By timepoint"),
               cm_bar_one("Phase", "By cell-cycle phase"), cm_bar_one("cycling", "By cycling"))
    dfs <- Map(function(p, nm) { d <- p$data; d$panel <- nm; d },
               ps, c("By genotype","By timepoint","By cell-cycle phase","By cycling"))
    tb <- do.call(rbind, dfs); tb$panel <- factor(tb$panel, levels = unique(tb$panel))
    tb$cat <- factor(tb$cat, levels = unique(tb$cat)); ch <- input$cmbar_palette %||% "Default"
    p <- ggplot(tb, aes(sub, y, fill = cat)) + geom_col() + facet_wrap(~ panel, ncol = 1, scales = "free_y") +
      theme_minimal(base_size = input$cmbar_basesize %||% 13) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      labs(x = "subcluster", y = if ((input$cm_bar_mode %||% "prop") == "prop") "proportion" else "cell count",
           fill = NULL, title = "Subcluster composition — res 0.2")
    if (ch != "Default") p <- p + scale_fill_manual(values = disc_pal(levels(tb$cat), ch))
    p
  })
  # Each panel is exportable on its own, not only via the combined `cmbar` figure --
  # the combined one exists for download convenience, but someone looking at the
  # phase panel wants the phase panel.
  cm_bar_ps <- list(cmbargeno  = reactive(cm_bar_one("genotype",  "By genotype (WT / KO)")),
                    cmbartp    = reactive(cm_bar_one("timepoint", "By timepoint (P0 / P7)")),
                    cmbarphase = reactive(cm_bar_one("Phase",     "By cell-cycle phase")),
                    cmbarcyc   = reactive(cm_bar_one("cycling",   "By cycling status")))
  output$cm_bar_geno  <- renderPlot(apply_fig_opts(cm_bar_ps$cmbargeno(),  "cmbar", input))
  output$cm_bar_tp    <- renderPlot(apply_fig_opts(cm_bar_ps$cmbartp(),    "cmbar", input))
  output$cm_bar_phase <- renderPlot(apply_fig_opts(cm_bar_ps$cmbarphase(), "cmbar", input))
  output$cm_bar_cyc   <- renderPlot(apply_fig_opts(cm_bar_ps$cmbarcyc(),   "cmbar", input))
  # Registered explicitly rather than in a loop: the coverage check in
  # tools/check_download_coverage.R reads the source, and a loop variable hides
  # the prefixes from it.
  register_fig(output, "cmbargeno",  cm_bar_ps$cmbargeno,  input, opts_prefix = "cmbar")
  register_fig(output, "cmbartp",    cm_bar_ps$cmbartp,    input, opts_prefix = "cmbar")
  register_fig(output, "cmbarphase", cm_bar_ps$cmbarphase, input, opts_prefix = "cmbar")
  register_fig(output, "cmbarcyc",   cm_bar_ps$cmbarcyc,   input, opts_prefix = "cmbar")
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("cmbar_dl_", f)]] <- dl_ggplot("cmbar", cm_bar_plot, input, f) })

  # ---- CM per-cluster summary sheet (res 0.2) ----
  cm_summary_df <- reactive({
    ss <- subSum[["res0.2"]]; st <- subType[["res0.2"]]; hl <- heat[["res0.2"]]$long
    df <- cmm; sub <- factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
    lv  <- levels(sub)
    n   <- as.integer(table(sub)[lv])
    pKO <- round(100 * tapply(df$genotype == "KO", sub, mean)[lv], 1)
    pP0 <- round(100 * tapply(df$timepoint == "P0", sub, mean)[lv], 1)
    mk  <- vapply(lv, function(cl) { s <- hl[hl$cluster == cl, ]; paste(head(s$gene[order(-s$z)], 5), collapse = ", ") }, "")
    topOf <- function(tbl, cl, key, out) {
      if (is.null(tbl)) return(NA_character_)
      d <- tbl[tbl$subcluster == cl & !is.na(tbl[[key]]), , drop = FALSE]
      if (!nrow(d)) return(NA_character_)
      as.character(d[[out]][which.min(d[[key]])])
    }
    idgo <- if (!is.null(ENR$sub)) ENR$sub$identity_go else NULL
    gsea <- if (!is.null(ENR$sub)) ENR$sub$gsea else NULL
    data.frame(
      cluster = lv,
      subtype = st$nearest_CM_subtype[match(lv, st$subcluster)],
      n_cells = n,
      pct_KO = pKO, pct_WT = round(100 - pKO, 1),
      pct_P0 = pP0, pct_P7 = round(100 - pP0, 1),
      top_markers = mk,
      n_DE = ss$n_DE_absLFC_gt1[match(lv, ss$subcluster)],
      top_KO_up = ss$top_KO_up[match(lv, ss$subcluster)],
      top_KO_down = ss$top_KO_down[match(lv, ss$subcluster)],
      top_identity_GO    = vapply(lv, function(c) topOf(idgo, c, "p.adjust", "Description"), ""),
      top_KOvsWT_pathway = vapply(lv, function(c) topOf(gsea, c, "padj", "pathway"), ""),
      check.names = FALSE, stringsAsFactors = FALSE)
  })
  output$cm_summary <- renderDT(DT::datatable(cm_summary_df(), rownames = FALSE,
    options = list(pageLength = 13, scrollX = TRUE, dom = "ft"), class = "compact stripe hover"))

  # ---- CM top markers + biology per subcluster (stylized; res 0.2) ----
  # Built as a plain-text frame first, then decorated with HTML chips for the
  # screen -- the chips are unreadable in a CSV, so the download takes the frame.
  cm_topmarkers_df <- reactive({
    hl <- heat[["res0.2"]]$long; st <- subType[["res0.2"]]; ss <- subSum[["res0.2"]]
    idgo <- if (!is.null(ENR$sub)) ENR$sub$identity_go else NULL
    gsea <- if (!is.null(ENR$sub)) ENR$sub$gsea else NULL
    flat <- function(x) {
      if (length(x) == 1 && is.character(x)) x <- trimws(strsplit(x, ",")[[1]])
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) "" else paste(x, collapse = ", ")
    }
    topOf <- function(tbl, cl, key, out) {
      if (is.null(tbl)) return("")
      d <- tbl[tbl$subcluster == cl & !is.na(tbl[[key]]), , drop = FALSE]
      if (!nrow(d)) return("")
      as.character(d[[out]][which.min(d[[key]])])
    }
    do.call(rbind, lapply(cm_subs("0.2"), function(cl) {
      m <- hl$cluster == cl; gs <- head(hl$gene[m][order(-hl$z[m])], 8)
      inSS <- !is.null(ss) && cl %in% ss$subcluster
      data.frame(
        Subcluster = cl,
        Subtype = if (!is.null(st) && cl %in% st$subcluster) st$nearest_CM_subtype[match(cl, st$subcluster)] else "",
        Cells   = if (inSS) ss$n_cells[match(cl, ss$subcluster)] else NA_integer_,
        `Top markers` = flat(gs),
        `Identity (top GO BP)` = topOf(idgo, cl, "p.adjust", "Description"),
        `KO-vs-WT (top GSEA pathway)` = topOf(gsea, cl, "padj", "pathway"),
        `KO-up genes`   = if (inSS) flat(ss$top_KO_up[match(cl, ss$subcluster)]) else "",
        `KO-down genes` = if (inSS) flat(ss$top_KO_down[match(cl, ss$subcluster)]) else "",
        check.names = FALSE, stringsAsFactors = FALSE)
    }))
  })
  output$cm_topmarkers <- renderDT({
    hl <- heat[["res0.2"]]$long; st <- subType[["res0.2"]]; ss <- subSum[["res0.2"]]
    idgo <- if (!is.null(ENR$sub)) ENR$sub$identity_go else NULL
    gsea <- if (!is.null(ENR$sub)) ENR$sub$gsea else NULL
    chip  <- function(g, bg) paste0("<span style='display:inline-block;background:", bg,
      ";border:1px solid #d0d7e2;border-radius:10px;padding:1px 8px;margin:2px;font-size:12px'>", g, "</span>")
    chips <- function(x, bg = "#eef2f7") {                        # gene vector OR comma-string -> chips
      if (length(x) == 1 && is.character(x)) x <- trimws(strsplit(x, ",")[[1]])
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return("—")
      paste(vapply(x, chip, "", bg = bg), collapse = " ")
    }
    topOf <- function(tbl, cl, key, out) {                        # top-ranked term for a subcluster
      if (is.null(tbl)) return("—")
      d <- tbl[tbl$subcluster == cl & !is.na(tbl[[key]]), , drop = FALSE]
      if (!nrow(d)) return("—")
      as.character(d[[out]][which.min(d[[key]])])
    }
    rows <- lapply(cm_subs("0.2"), function(cl) {
      m <- hl$cluster == cl; gs <- head(hl$gene[m][order(-hl$z[m])], 8)
      inSS <- !is.null(ss) && cl %in% ss$subcluster
      data.frame(
        Subcluster = cl,
        Subtype = if (!is.null(st) && cl %in% st$subcluster) st$nearest_CM_subtype[match(cl, st$subcluster)] else "—",
        Cells   = if (inSS) ss$n_cells[match(cl, ss$subcluster)] else NA_integer_,
        `Top markers` = chips(gs),
        `Identity (top GO BP)` = topOf(idgo, cl, "p.adjust", "Description"),
        `KO-vs-WT (top GSEA pathway)` = topOf(gsea, cl, "padj", "pathway"),
        `KO-up genes`   = if (inSS) chips(ss$top_KO_up[match(cl, ss$subcluster)], "#fdecea") else "—",
        `KO-down genes` = if (inSS) chips(ss$top_KO_down[match(cl, ss$subcluster)], "#e8f0fe") else "—",
        check.names = FALSE, stringsAsFactors = FALSE)
    })
    DT::datatable(do.call(rbind, rows), rownames = FALSE, escape = FALSE, selection = "none",
      options = list(pageLength = 15, dom = "t", ordering = FALSE, scrollX = TRUE),
      class = "compact stripe hover")
  })

  # ---- CM per-subcluster enrichment (ENR$sub pooled, or FGE per contrast; res 0.2) ----
  # Which enrichment source the tab is reading. NULL = the pooled KO-vs-WT tables.
  cm_enr_ct <- reactive(if (identical(input$cm_enr_contrast %||% "pooled", "pooled")) NULL
                        else fg_ct(input$cm_enr_contrast))
  cm_enr_a  <- reactive({ ct <- cm_enr_ct()
    list(ct = ct, st = if (is.null(ct)) "all" else (input$cm_enr_stratum %||% "all"),
         ont = if (is.null(ct)) "BP" else (input$cm_enr_ont %||% "BP"),
         lab = enr_src_labs(ct)) })
  # Availability differs by source, so the cluster list follows the comparison. isolate()
  # on the current value: this observer writes the input it reads.
  observe({
    ct <- cm_enr_ct()
    avail <- if (is.null(ct)) {
      if (!is.null(ENR$sub)) sort(unique(c(ENR$sub$identity_go$subcluster,
                                           ENR$sub$go$subcluster, ENR$sub$gsea$subcluster))) else character(0)
    } else if (!is.null(FGE$go)) unique(FGE$go$cluster[FGE$go$contrast == ct$key]) else character(0)
    subs <- cm_subs("0.2"); subs <- if (length(avail)) subs[subs %in% avail] else subs
    if (!length(subs)) return()
    cur  <- isolate(input$cm_enr_sub)
    updateSelectInput(session, "cm_enr_sub",
                      choices = setNames(subs, vapply(subs, function(s) sub_label("0.2", s), "")),
                      selected = if (!is.null(cur) && cur %in% subs) cur else subs[1])
  })
  output$cm_enr_note <- renderUI({ a <- cm_enr_a(); req(input$cm_enr_sub)
    div(class = "alert alert-light border py-2 px-3 mb-2",
        tags$small(strong(a$lab$label), " · ", input$cm_enr_sub,
                   if (is.null(a$ct)) NULL else paste0(" · ", a$st, " cells · GO ", a$ont),
                   sprintf(" · “%s” vs “%s”", a$lab$up, a$lab$dn),
                   if (is.null(a$ct))
                     span(" · pooled over P0 and P7 — this cannot answer what changes in the KO at P7.")
                   else NULL)) })
  # ---- identity GO: contrast-independent, always the pooled marker enrichment ----
  cm_sub_idgo_p <- reactive({ req(input$cm_enr_sub)
    do.call(go_dotplot_gg, c(list(enr_sub_df("identity_go", input$cm_enr_sub),
                                  paste0("Identity GO BP — ", input$cm_enr_sub)), cm_sub_opts())) })
  output$cm_sub_idgo_plot <- renderPlotly(
    ggplotly(cm_sub_idgo_p(), tooltip = "text") |> layout(margin = list(l = 0, t = 40)))
  cm_sub_idgo_df <- reactive({ req(input$cm_enr_sub)
    d <- enr_sub_df("identity_go", input$cm_enr_sub)
    d[order(d$p.adjust), intersect(c("ID","Description","FoldEnrichment","p.adjust","Count","geneID"), names(d))] })
  output$cm_sub_idgo_tab <- renderDT(enr_dt(cm_sub_idgo_df()))
  # ---- contrast GO, the two directions enriched separately ----
  CM_GO_COLS <- c("ID","Description","direction","direction_label","FoldEnrichment",
                  "p.adjust","qvalue","Count","geneID","n_input","n_universe","input_rule")
  cm_sub_go <- function(dir) { a <- cm_enr_a(); req(input$cm_enr_sub)
    enr_src_df("go", input$cm_enr_sub, a$ct, a$st, a$ont, dir) }
  cm_sub_go_ttl <- function(dir) { a <- cm_enr_a()
    sprintf("GO %s — %s, %s", a$ont, input$cm_enr_sub,
            if (identical(dir, "up")) a$lab$up else a$lab$dn) }
  # An empty GO table is ambiguous: "tested, nothing passed" and "list too small to test"
  # look identical. The four-group audit knows which; say so.
  cm_sub_go_empty <- function(dir) { a <- cm_enr_a()
    lab <- if (identical(dir, "up")) a$lab$up else a$lab$dn
    if (is.null(a$ct)) return(paste0("No GO BP term for the ", lab, " list in this subcluster."))
    paste0("No GO ", a$ont, " term reached padj < 0.05 / q < 0.2 for the ", lab, " list. ",
           fg_enr_why(input$cm_enr_sub, a$ct$key, a$st, a$ont,
                      if (identical(dir, "up")) "A_up" else "B_up")) }
  # One options bundle for all four subcluster-enrichment views (prefix "cmsubenr").
  # In the "All clusters" grid the panels are half or full width rather than a whole
  # card, so the label budget is tighter there; `compact` applies that without needing a
  # second set of controls.
  cm_sub_opts <- function(compact = FALSE) {
    lc <- input$cmsubenr_labelchars %||% 46
    list(base_size  = input$cmsubenr_basesize %||% 11,
         pal_choice = input$cmsubenr_palette,
         axis_scale = input$cmsubenr_axisscale %||% 1,
         colour_by  = input$cmsubenr_colourby %||% "padj",
         label_chars = if (compact) min(lc, if ((input$cm_enr_perrow %||% "1") == "2") 26 else 38) else lc)
  }
  cm_sub_go_plot <- function(dir) { d <- cm_sub_go(dir)
    validate(need(nrow(d), cm_sub_go_empty(dir)))
    do.call(go_dotplot_gg, c(list(d, cm_sub_go_ttl(dir)), cm_sub_opts())) }
  cm_sub_go_ly <- function(dir)                      # same body go_dotplot_df() uses
    ggplotly(cm_sub_go_plot(dir), tooltip = "text") |> layout(margin = list(l = 0, t = 40))
  output$cm_sub_kogo_plot <- renderPlotly(cm_sub_go_ly("up"))
  output$cm_sub_kodn_plot <- renderPlotly(cm_sub_go_ly("down"))
  cm_sub_kogo_df <- reactive({ d <- cm_sub_go("up")
    d[order(d$p.adjust), intersect(CM_GO_COLS, names(d))] })
  cm_sub_kodn_df <- reactive({ d <- cm_sub_go("down")
    d[order(d$p.adjust), intersect(CM_GO_COLS, names(d))] })
  output$cm_sub_kogo_tab <- renderDT(enr_dt(cm_sub_kogo_df()))
  output$cm_sub_kodn_tab <- renderDT(enr_dt(cm_sub_kodn_df()))
  # ---- GSEA: same source switch; up/down wording from the contrast ----
  cm_sub_gsea_dat <- reactive({ a <- cm_enr_a(); req(input$cm_enr_sub)
    enr_src_df("gsea", input$cm_enr_sub, a$ct, a$st) })
  cm_sub_gsea_p <- reactive({ a <- cm_enr_a(); o <- cm_sub_opts()
    o$pal_choice <- NULL; o$colour_by <- NULL      # bar chart: up/down pair, no ramp
    do.call(gsea_barplot_gg,
            c(list(cm_sub_gsea_dat(),
                   paste0("GSEA — ", input$cm_enr_sub, " · ", a$lab$label),
                   20, up_lab = a$lab$up, down_lab = a$lab$dn), o)) })
  output$cm_sub_gsea_plot <- renderPlotly(
    ggplotly(cm_sub_gsea_p(), tooltip = "text") |> layout(margin = list(l = 0, t = 40)))
  cm_sub_gsea_df <- reactive({ d <- cm_sub_gsea_dat()
    d[order(d$padj), intersect(c("pathway","NES","padj","pval","size","leadingEdge"), names(d))] })
  output$cm_sub_gsea_tab <- renderDT(enr_dt(cm_sub_gsea_df()))
  # filename stem for the enrichment downloads -- carries the comparison, not a fixed "pooled"
  cm_enr_dl_base <- function(kind, suffix) { a <- cm_enr_a()
    paste0(kind, "_",
           if (is.null(a$ct)) "KOvsWT_pooled" else paste0(a$ct$key, "_", a$st),
           if (identical(kind, "GO") && !is.null(a$ct)) paste0("_", a$ont) else "",
           "_", input$cm_enr_sub, suffix) }
  # the same reactive the panel renders, so download and Figure Studio get what is shown
  register_fig(output, "cmsubidgo", cm_sub_idgo_p, input)
  register_fig(output, "cmsubkogo", reactive(cm_sub_go_plot("up")), input)
  register_fig(output, "cmsubkodn", reactive(cm_sub_go_plot("down")), input)
  register_fig(output, "cmsubgsea", cm_sub_gsea_p, input)
  for (.k in c("idgo", "kogo", "kodn", "gsea")) local({
    k <- .k
    output[[paste0("cm_grid_", k)]] <- renderUI(
      cm_enr_grid(k, per_row = as.integer(input$cm_enr_perrow %||% "1"),
                  height = if ((input$cm_enr_perrow %||% "1") == "2") "360px" else "420px"))
  })
  # "All clusters" view: one plot per subcluster, 1 or 2 per row (see cm_enr_perrow).
  # These read cm_enr_a() too, so the grid honours the comparison like the single view.
  local({
    for (.cl in cm_subs("0.2")) local({
      cl <- .cl; ttl <- sub_label("0.2", cl)
      go_all <- function(dir) { a <- cm_enr_a()
        d <- enr_src_df("go", cl, a$ct, a$st, a$ont, dir)
        lab <- if (identical(dir, "up")) a$lab$up else a$lab$dn
        validate(need(nrow(d), paste0(cl, ": no GO term for the ", lab, " list.")))
        do.call(go_dotplot_gg, c(list(d, paste0(ttl, " — ", lab), topn = 8), cm_sub_opts(TRUE))) }
      output[[paste0("cm_idgo_all_", cl)]] <- renderPlot(
        do.call(go_dotplot_gg, c(list(enr_sub_df("identity_go", cl), ttl, topn = 8), cm_sub_opts(TRUE))))
      output[[paste0("cm_kogo_all_", cl)]] <- renderPlot(go_all("up"))
      output[[paste0("cm_kodn_all_", cl)]] <- renderPlot(go_all("down"))
      output[[paste0("cm_gsea_all_", cl)]] <- renderPlot({ a <- cm_enr_a()
        o <- cm_sub_opts(TRUE); o$pal_choice <- NULL; o$colour_by <- NULL   # bar chart: up/down pair
        do.call(gsea_barplot_gg, c(list(enr_src_df("gsea", cl, a$ct, a$st), ttl, topn = 10,
                                        up_lab = a$lab$up, down_lab = a$lab$dn), o)) })
    })
  })

  # ---- E2F focus (E2f7 / E2f8 expression by genotype x timepoint) ----
  e2f_plot <- reactive({
    req(input$e2f_ct)
    bs <- input$e2f_basesize %||% 13; ch <- input$e2f_palette %||% "Default"; rn <- e2f_rn()
    df <- meta
    if (input$e2f_ct != "All" && has("celltype")) df <- df[df$celltype == input$e2f_ct, ]
    validate(need(nrow(df) > 0, "No cells for this cell type."))
    eg <- intersect(c("E2f7","E2f8"), c(rownames(expr), if (!is.null(EXPR)) rownames(EXPR)))
    validate(need(length(eg) > 0, "E2f7/E2f8 not present in this data build."))
    long <- do.call(rbind, lapply(eg, function(g) {
      v <- expr_vec(g, df$cell)
      data.frame(gene = g, expr = v, genotype = df$genotype,
                 timepoint = if (has("timepoint", df)) df$timepoint else "all",
                 stringsAsFactors = FALSE)
    }))
    long <- long[is.finite(long$expr), ]
    validate(need(nrow(long) > 0, "No expression values (gene only in the broad subset with no overlap here)."))
    glv <- levels(factor(meta$genotype))
    long$genotype <- factor(relab(long$genotype, rn), levels = relab(glv, rn))
    fillv <- if (ch == "Default") setNames(c("#c62828","#1565c0"), relab(c("KO","WT"), rn))
             else disc_pal(relab(glv, rn), ch)
    ggplot(long, aes(genotype, expr, fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, alpha = .55, linewidth = .2) +
      stat_summary(fun = mean, geom = "point", size = 2.4, color = "black") +
      facet_grid(gene ~ timepoint, scales = "free_y") +
      scale_fill_manual(values = fillv, na.value = "grey70") +
      theme_minimal(base_size = bs) + guides(fill = "none") +
      labs(x = "genotype", y = "log-norm expression",
           title = paste0("E2f7 / E2f8 — ", gsub("_", " ", input$e2f_ct),
                          " (black dot = mean; descriptive, n = 1)"))
  })
  output$e2f_expr <- renderPlot(apply_fig_opts(e2f_plot(), "e2f", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("e2f_dl_", f)]] <- dl_ggplot("e2f", e2f_plot, input, f) })

  # downstream E2F targets: gene set (from the shared GENE_SETS list), panel genes only
  e2f_targets <- reactive({
    s <- input$e2f_set %||% "E2F targets"; if (!s %in% names(GENE_SETS)) s <- "E2F targets"
    intersect(GENE_SETS[[s]], rownames(expr))
  })
  # per-gene KO-vs-WT fold change of the downstream targets (from precomputed DE)
  e2f_fc_plot <- reactive({
    tg <- e2f_targets(); validate(need(length(tg) >= 1, "No target genes."))
    ct <- if (input$e2f_ct == "All" || is.null(input$e2f_ct)) "Cardiomyocyte" else input$e2f_ct
    bs <- input$e2ffc_basesize %||% 13
    rows <- do.call(rbind, lapply(c("P0","P7"), function(tp) {
      d <- ctDE[[paste(tp, ct, sep = "_")]]; if (is.null(d)) return(NULL)
      i <- match(tg, d$gene); i <- i[!is.na(i)]
      if (!length(i)) return(NULL)
      data.frame(gene = d$gene[i], timepoint = tp, lfc = d$log2FoldChange[i], padj = d$padj[i])
    }))
    validate(need(!is.null(rows) && nrow(rows), sprintf("No DE for %s targets.", ct)))
    rows$dir <- ifelse(rows$lfc >= 0, "up in KO", "up in WT")
    rows$gene <- factor(rows$gene, levels = unique(rows$gene[order(ave(rows$lfc, rows$gene, FUN = mean))]))
    ggplot(rows, aes(lfc, gene, color = dir)) +
      geom_vline(xintercept = 0, color = "grey70") +
      geom_segment(aes(x = 0, xend = lfc, yend = gene), linewidth = .5) + geom_point(size = 2) +
      facet_wrap(~timepoint) + scale_color_manual(values = VOLC_PAL, guide = guide_legend(title = NULL)) +
      theme_minimal(base_size = bs) +
      labs(x = "log2 fold change (KO / WT)", y = NULL,
           title = paste0(input$e2f_set, " — KO vs WT in ", gsub("_", " ", ct),
                          if (input$e2f_ct == "All") " (all-cells → CM shown)" else ""))
  })
  output$e2f_fc <- renderPlot(apply_fig_opts(e2f_fc_plot(), "e2ffc", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("e2ffc_dl_", f)]] <- dl_ggplot("e2ffc", e2f_fc_plot, input, f) })

  # ---- Four-group (WT/KO x P0/P7) within CM subclusters ----
  # Populate the selectors from the bundle; an un-rebuilt bundle leaves them empty
  # and every output falls through to the "run build_fourgroup.R" validate().
  observe({
    cc <- fg_cluster_choices()
    if (!length(cc)) return()
    dflt <- if ("CM2" %in% cc) "CM2" else cc[[1]]
    updateSelectInput(session, "fg_cluster", choices = cc, selected = dflt)
    updateSelectInput(session, "mi_cluster", choices = cc, selected = dflt)
    updateSelectInput(session, "fg_contrast", choices = fg_contrast_choices(),
                      selected = "P7_KO_vs_WT")
    pri <- intersect(c("AllCM", FG_PRIORITY, "CM9"), FG_CLUSTERS)
    updateSelectizeInput(session, "fg_g1_clusters", choices = cc, selected = pri)
    updateSelectizeInput(session, "mi_cand_clusters", choices = cc, selected = pri)
    sc <- intersect(c("sig_maturation_nocc","sig_maturation","sig_mat_mature","sig_mat_immature",
                      "sig_metabolic","sig_prolif","sig_cytokinesis","sig_ccexit"), names(cmm))
    if (length(sc)) updateSelectInput(session, "fg_score",
                      choices = setNames(sc, labof(sc)), selected = sc[[1]])
  })
  observe({
    pool <- fg_candidate_pool(input$mi_source %||% "shortlist", input$mi_topn %||% 20,
                              input$mi_spec_grid %||% "de")
    if (!is.null(input$mi_geneset) && input$mi_geneset != "__all__")
      pool <- intersect(pool, genes_for_set(input$mi_geneset))
    updateSelectizeInput(session, "mi_genes", choices = genes_for_set(input$mi_geneset),
                         selected = head(pool, 60), server = TRUE)
  })
  observeEvent(input$mi_reset_genes, {
    updateSelectInput(session, "mi_source", selected = "shortlist")
    updateSelectizeInput(session, "mi_genes", selected = FG_SHORTLIST, server = TRUE)
  })

  # group sizes
  fg_counts_p <- reactive(fg_counts_plot(input$fg_count_mode %||% "prop",
                                         input$fgcount_basesize %||% 13))
  output$fg_counts_plot <- renderPlot(apply_fig_opts(fg_counts_p(), "fgcount", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("fgcount_dl_", f)]] <- dl_ggplot("fgcount", fg_counts_p, input, f) })
  output$fg_counts_tab <- renderDT(DT::datatable(fg_counts_wide(), rownames = FALSE,
    options = list(pageLength = 14, scrollX = TRUE, dom = "ft"),
    class = "compact stripe hover"))

  # four-group DE: volcano <-> table <-> gene card (same wiring as the other DE tabs)
  fg_d    <- reactive(drop_conf(fg_de(input$fg_cluster, input$fg_contrast, input$fg_stratum,
                                      input$fg_grid %||% "de"), input$fg_hideconf))
  fg_tab  <- reactive(de_table(fg_d()))
  fg_pick <- reactiveVal(NULL)
  fg_dt_proxy <- DT::dataTableProxy("fg_detab")
  # An arm can clear the 10-cell floor and still be far too thin to trust; de_context_note()
  # says so here, and on the CM deep-dive DE tab, from a single definition.
  output$fg_de_note <- renderUI({
    ct <- fg_ct(input$fg_contrast); req(ct)
    d <- try(fg_d(), silent = TRUE)
    de_context_note(input$fg_cluster, ct, input$fg_stratum, input$fg_grid,
                    if (inherits(d, "try-error")) NULL else d)
  })
  output$fg_volcano <- renderPlotly({
    ct <- fg_ct(input$fg_contrast); req(ct)
    de_volcano_ly(fg_d(),
      paste0(ct$label, " — ", if (input$fg_cluster == "AllCM") "all CM" else input$fg_cluster,
             if (input$fg_stratum == "G1") " (G1)" else ""),
      "fg_volcano", pos = ct$pos, neg = ct$neg, xlab = ct$xlab, highlight = fg_pick(),
      lfc_cut = input$fg_vlfc %||% DE_LFC_CUT)
  })
  outputOptions(output, "fg_volcano", suspendWhenHidden = FALSE)  # register the click source at startup
  # static twin for vector export + studio (same contrast, title and colour cut)
  fg_volcano_gg <- reactive({ ct <- fg_ct(input$fg_contrast); req(ct)
    de_volcano(fg_d(),
      paste0(ct$label, " — ", if (input$fg_cluster == "AllCM") "all CM" else input$fg_cluster,
             if (input$fg_stratum == "G1") " (G1)" else ""),
      pos = ct$pos, neg = ct$neg, xlab = ct$xlab,
      lfc_cut = input$fg_vlfc %||% DE_LFC_CUT, highlight = fg_pick()) })
  register_fig(output, "fgvolc", fg_volcano_gg, input)
  output$fg_detab    <- renderDT(de_datatable(fg_tab(), scroll = NULL))
  output$fg_pick_ui  <- renderUI(pick_banner(fg_pick(), "fg_clear"))
  output$fg_geneinfo <- renderUI(gene_info_card(fg_pick(), "fg_infoclose"))
  observeEvent(event_data("plotly_click", source = "fg_volcano"),
    fg_pick(event_data("plotly_click", source = "fg_volcano")$customdata))
  observeEvent(input$fg_detab_rows_selected, {
    r <- input$fg_detab_rows_selected; g <- if (length(r)) fg_tab()$gene[r] else NULL
    if (!identical(g, fg_pick())) fg_pick(g)
  }, ignoreNULL = FALSE)
  observeEvent(fg_pick(), {
    g <- fg_pick(); rows <- if (!is.null(g)) which(fg_tab()$gene == g) else integer(0)
    if (!identical(as.integer(rows), as.integer(input$fg_detab_rows_selected)))
      DT::selectRows(fg_dt_proxy, if (length(rows)) rows else NULL)
  }, ignoreNULL = FALSE)
  observeEvent(input$fg_clear, fg_pick(NULL))
  observeEvent(input$fg_infoclose, fg_pick(NULL))
  observeEvent(list(input$fg_cluster, input$fg_contrast, input$fg_stratum, input$fg_grid),
               fg_pick(NULL), ignoreInit = TRUE)

  # G1 proportion + maturation scores
  fg_phase_p <- reactive(fg_phase_plot(input$fg_g1_clusters, input$fgphase_basesize %||% 13))
  output$fg_phase_plot <- renderPlot(apply_fig_opts(fg_phase_p(), "fgphase", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("fgphase_dl_", f)]] <- dl_ggplot("fgphase", fg_phase_p, input, f) })
  fg_score_p <- reactive({ req(input$fg_score)
    fg_score_plot(input$fg_score, input$fg_g1_clusters, input$fg_score_stratum %||% "all",
                  input$fgscore_basesize %||% 13) })
  output$fg_score_plot <- renderPlot(apply_fig_opts(fg_score_p(), "fgscore", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("fgscore_dl_", f)]] <- dl_ggplot("fgscore", fg_score_p, input, f) })
  fg_summary_p <- reactive(fg_summary_plot(input$fg_score %||% "sig_maturation_nocc",
    input$fg_score_stratum %||% "all", input$fg_g1_clusters, input$fgscore_basesize %||% 13))
  output$fg_summary_plot <- renderPlot(apply_fig_opts(fg_summary_p(), "fgscore", input))
  register_fig(output, "fgsumm", fg_summary_p, input, opts_prefix = "fgscore")
  output$fg_summary_tab <- renderDT({
    d <- fg_summary_df(input$fg_score %||% "sig_maturation_nocc", input$fg_score_stratum %||% "all")
    DT::datatable(d, rownames = FALSE,
      options = list(pageLength = 14, scrollX = TRUE, dom = "ft"),
      class = "compact stripe hover") |>
      DT::formatSignif(intersect(c("mat_WT_P7","mat_KO_P7","mat_gap_P7","d_P7",
                                   "mat_gap_P0","d_P0"), names(d)), 3)
  })
  # No table on screen for this one -- it is the per-cell score summary behind the
  # G1/maturation plots -- so it gets the same csv+xlsx pair rather than a bespoke
  # handler. Same for mi_cand below (the expression grid behind the dot plot).
  register_dl(output, "fg_scores", function() {
    validate(need(!is.null(FG$scores), "No score summary in this data build."))
    FG$scores
  }, "fourgroup_scores")

  # Whole-contrast workbook. Built on demand rather than precomputed: it is one
  # click for a user and ~10 s of work, and precomputing 8 contrasts x 2 strata
  # of xlsx into the bundle would be absurd.
  # One sheet per subcluster for a contrast. A named function rather than an
  # inline handler body so tools/test_downloads.R can build the workbook without
  # a browser.
  fg_book_sheets <- function(contrast, stratum, grid = "de", clusters = NULL, hideconf = FALSE) {
    fg_ok(); req(contrast, stratum)
    cls <- clusters %||% FG_CLUSTERS
    out <- list()
    for (cl in cls) {
      d <- fg_grid(grid)[[cl]][[paste0(contrast, "__", stratum)]]
      if (is.null(d) || !nrow(d)) next
      out[[paste0("DE_", cl)]] <- drop_conf(d, hideconf)
    }
    out$Group_sizes <- fg_counts_wide()
    if (!is.null(FG$skipped)) out$Skipped <- fg_skipset(grid)
    out
  }
  output$fg_book <- dl_book(
    basename = function() paste0("fourgroup_", input$fg_contrast %||% "DE", "_",
                                 input$fg_stratum %||% "G1",
                                 if (identical(input$fg_grid, "de2")) "_curated" else ""),
    sheets_fn = function()
      fg_book_sheets(input$fg_contrast, input$fg_stratum, input$fg_grid %||% "de",
                     clusters = if (isTRUE(input$fg_book_all)) NULL else input$fg_cluster,
                     hideconf = isTRUE(input$fg_hideconf)),
    title = function() {
      lab <- if (!is.null(FG_CTAB)) FG_CTAB$label[match(input$fg_contrast, FG_CTAB$key)] else input$fg_contrast
      paste0("Four-group cardiomyocyte DE - ", lab, " (", input$fg_stratum, " cells)")
    },
    notes = c(
      "One sheet per CM subcluster. log2FoldChange > 0 means up in the first arm of the comparison.",
      "Row gate: these tables come from build_fourgroup.R, which keeps a gene when it is expressed in >= 5% of one arm AND (padj < 0.05 OR |log2FC| >= 0.5). They are therefore NOT complete gene lists -- for every tested gene, use the analysis/ pipeline's CSVs.",
      "'G1' sheets are phase-matched between the arms; 'all' sheets use every cell.")
  )

  # ---- Four-group enrichment (precomputed by build_fourgroup_enrichment.R) ----
  # Its own cluster/contrast selectors rather than reusing the DE tab's: you
  # routinely want to read the DE table for one subcluster while comparing
  # enrichment across another, and sharing the inputs would couple them.
  observe({
    cc <- fg_cluster_choices()
    if (!length(cc)) return()                       # un-rebuilt bundle: leave empty
    avail <- if (!is.null(FGE$go)) unique(FGE$go$cluster) else character(0)
    sel   <- if (length(avail)) cc[cc %in% avail] else cc
    updateSelectInput(session, "fg_enr_cluster", choices = sel,
                      selected = if ("CM2" %in% sel) "CM2" else sel[1])
    updateSelectInput(session, "fg_enr_contrast", choices = fg_contrast_choices(),
                      selected = "P7_KO_vs_WT")
  })
  fg_enr_args <- reactive({
    req(input$fg_enr_cluster, input$fg_enr_contrast, input$fg_enr_stratum)
    ct <- fg_ct(input$fg_enr_contrast)
    list(cl = input$fg_enr_cluster, ck = input$fg_enr_contrast, st = input$fg_enr_stratum,
         ont = input$fg_enr_ont %||% "BP", topn = input$fg_enr_topn %||% 20,
         up = ct$pos %||% "up in A", dn = ct$neg %||% "up in B", label = ct$label %||% input$fg_enr_contrast)
  })
  fg_enr_go_df <- function(direction) {
    a <- fg_enr_args()
    d <- fg_enr_df("go", a$cl, a$ck, a$st, a$ont, direction)
    d[order(d$p.adjust), intersect(c("ID","Description","direction_label","FoldEnrichment",
                                     "p.adjust","qvalue","Count","geneID","n_input","n_universe",
                                     "n_mt_input","input_rule"), names(d))]
  }
  fg_enr_up_df    <- reactive(fg_enr_go_df("A_up"))
  fg_enr_dn_df    <- reactive(fg_enr_go_df("B_up"))
  fg_enr_gsea_dat <- reactive({ a <- fg_enr_args()
    d <- fg_enr_df("gsea", a$cl, a$ck, a$st)
    d[order(d$padj), intersect(c("pathway","NES","padj","pval","size","leadingEdge"), names(d))] })
  fg_enr_audit_dat <- reactive({ fg_enr_ok(); a <- fg_enr_args()
    d <- FGE$audit
    d[d$contrast == a$ck & d$stratum == a$st, , drop = FALSE] })

  # An empty GO table is ambiguous -- "tested, nothing passed" and "list too
  # small to test" look identical. Say which.
  fg_enr_empty <- function(direction) {
    a <- fg_enr_args()
    paste0("No GO ", a$ont, " term reached padj < 0.05 / q < 0.2 for the ", 
           if (direction == "A_up") a$up else a$dn, " list. ",
           fg_enr_why(a$cl, a$ck, a$st, a$ont, direction))
  }
  fg_enr_go_plot <- function(direction) {
    a <- fg_enr_args()
    d <- fg_enr_df("go", a$cl, a$ck, a$st, a$ont, direction)
    validate(need(nrow(d), fg_enr_empty(direction)))
    go_dotplot_df(d, sprintf("GO %s — %s, %s (%s cells)", a$ont, a$cl,
                             if (direction == "A_up") a$up else a$dn, a$st), a$topn)
  }
  output$fg_enr_up_plot <- renderPlotly(fg_enr_go_plot("A_up"))
  output$fg_enr_dn_plot <- renderPlotly(fg_enr_go_plot("B_up"))
  output$fg_enr_up_tab  <- renderDT(enr_dt(fg_enr_up_df()))
  output$fg_enr_dn_tab  <- renderDT(enr_dt(fg_enr_dn_df()))
  output$fg_enr_gsea_plot <- renderPlotly({ a <- fg_enr_args()
    # up_lab/down_lab from the contrast table: "up in KO" is wrong for WT: P0 vs P7.
    gsea_barplot_df(fg_enr_df("gsea", a$cl, a$ck, a$st),
                    sprintf("GSEA (Hallmark + KEGG) — %s, %s", a$cl, a$label),
                    a$topn, up_lab = a$up, down_lab = a$dn) })
  output$fg_enr_gsea_tab  <- renderDT(enr_dt(fg_enr_gsea_dat()))
  output$fg_enr_audit_tab <- renderDT(enr_dt(fg_enr_audit_dat(), scroll = "420px"))
  output$fg_enr_note <- renderUI({ a <- fg_enr_args()
    n <- if (is.null(FGE$go)) 0 else sum(FGE$go$cluster == a$cl & FGE$go$contrast == a$ck &
                                         FGE$go$stratum == a$st & FGE$go$ontology == a$ont)
    div(class = "alert alert-light border py-2 px-3 mb-2",
        tags$small(strong(a$label), " · ", a$cl, " · ", a$st, " cells · ",
                   sprintf("%d GO %s terms · ", n, a$ont),
                   sprintf("“%s” vs “%s”", a$up, a$dn)))
  })
  register_fig(output, "fgenrup", reactive({ a <- fg_enr_args()
    d <- fg_enr_df("go", a$cl, a$ck, a$st, a$ont, "A_up")
    validate(need(nrow(d), fg_enr_empty("A_up")))
    go_dotplot_gg(d, sprintf("GO %s — %s, %s", a$ont, a$cl, a$up), a$topn) }), input)
  register_fig(output, "fgenrdn", reactive({ a <- fg_enr_args()
    d <- fg_enr_df("go", a$cl, a$ck, a$st, a$ont, "B_up")
    validate(need(nrow(d), fg_enr_empty("B_up")))
    go_dotplot_gg(d, sprintf("GO %s — %s, %s", a$ont, a$cl, a$dn), a$topn) }), input)
  register_fig(output, "fgenrgsea", reactive({ a <- fg_enr_args()
    gsea_barplot_gg(fg_enr_df("gsea", a$cl, a$ck, a$st),
                    sprintf("GSEA — %s, %s", a$cl, a$label), a$topn,
                    up_lab = a$up, down_lab = a$dn) }), input)

  # Whole-contrast enrichment workbook: every subcluster, both directions.
  fg_enr_book_sheets <- function(contrast, stratum, ont) {
    fg_enr_ok()
    g <- FGE$go[FGE$go$contrast == contrast & FGE$go$stratum == stratum &
                FGE$go$ontology == ont, , drop = FALSE]
    s <- if (!is.null(FGE$gsea)) FGE$gsea[FGE$gsea$contrast == contrast &
                                          FGE$gsea$stratum == stratum, , drop = FALSE] else NULL
    a <- FGE$audit[FGE$audit$contrast == contrast & FGE$audit$stratum == stratum, , drop = FALSE]
    out <- list()
    for (cl in unique(g$cluster)) {
      x <- g[g$cluster == cl, , drop = FALSE]
      out[[paste0("GO_", cl)]] <- x[order(x$direction, x$p.adjust), ]
    }
    if (!is.null(s) && nrow(s)) out$GSEA <- s[order(s$cluster, s$padj), ]
    out$Coverage_audit <- a
    out
  }
  output$fg_enr_book <- dl_book(
    basename = function() paste0("fourgroup_enrichment_", input$fg_enr_contrast %||% "DE",
                                 "_", input$fg_enr_stratum %||% "all",
                                 "_", input$fg_enr_ont %||% "BP"),
    sheets_fn = function() fg_enr_book_sheets(input$fg_enr_contrast, input$fg_enr_stratum,
                                              input$fg_enr_ont %||% "BP"),
    title = function() { a <- fg_enr_args()
      paste0("GO/GSEA enrichment - ", a$label, " (", a$st, " cells, ", a$ont, ")") },
    notes = c(
      "One sheet per CM subcluster. The two directions are enriched SEPARATELY and are distinguished by the 'direction_label' column.",
      "Universe = genes detected in >= 5% of one arm in that cluster, recomputed from the expression matrix rather than taken from the (row-gated) DE table.",
      "'n_mt_input' counts mitochondrially-encoded genes in the input list. Where it is non-zero on an up-in-KO list, treat oxidative-phosphorylation and electron-transport terms with care - see the README caveats.",
      "The Coverage_audit sheet says how many genes went into each test and which selection rule fired, so an empty result can be told apart from an untested one.")
  )

  # ---- WT programs x P7 KO clusters (the four crossings from the email) ----
  observe({
    cc <- fg_cluster_choices()
    if (!length(cc)) return()
    updateSelectInput(session, "xc_wt_cluster", choices = cc, selected = "AllCM")
    subs <- setdiff(FG_CLUSTERS, "AllCM")
    lab  <- setNames(subs, vapply(subs, function(s) sub_label("0.2", s), ""))
    updateSelectizeInput(session, "xc_mat_clusters", choices = lab,
                         selected = intersect(XC_MAT_CLUSTERS, subs))
    updateSelectizeInput(session, "xc_cyc_clusters", choices = lab,
                         selected = intersect(XC_CYC_CLUSTERS, subs))
  })
  xc_p <- reactive({
    fg_ok()
    p <- list(
      wt_cluster   = input$xc_wt_cluster %||% "AllCM",
      mat_clusters = input$xc_mat_clusters %||% XC_MAT_CLUSTERS,
      cyc_clusters = input$xc_cyc_clusters %||% XC_CYC_CLUSTERS,
      mat_set      = input$xc_mat_set %||% "CM maturation",
      cyc_set      = input$xc_cyc_set %||% XC_CANON,
      stratum      = input$xc_stratum %||% "all",
      grid         = input$xc_grid %||% "de",
      padj         = input$xc_padj %||% 0.05,
      measure      = input$xc_measure %||% "auc",
      eff_auc      = input$xc_auc %||% 0.60,
      eff_lfc      = input$xc_lfc %||% 0.25,
      minc         = input$xc_minc %||% 1,
      hide_mt      = isTRUE(input$xc_hidemt))
    p$eff <- if (identical(p$measure, "auc")) p$eff_auc else p$eff_lfc
    validate(need(length(p$mat_clusters) && length(p$cyc_clusters),
                  "Pick at least one subcluster in each group."))
    p
  })
  xc_all <- reactive({ p <- xc_p(); lapply(XC_COMPARISONS, xc_comparison, p = p) })

  # Every way this tab can be misread, said out loud. Each of these has a specific
  # failure mode behind it, not a general disclaimer.
  output$xc_caveat <- renderUI({
    p <- try(xc_p(), silent = TRUE); if (inherits(p, "try-error")) return(NULL)
    red <- function(...) helpText(style = "color:#c62828", ...)
    gone <- c(setdiff(p$mat_clusters, xc_ko_present(p$mat_clusters, p$stratum, p$grid)),
              setdiff(p$cyc_clusters, xc_ko_present(p$cyc_clusters, p$stratum, p$grid)))
    thin <- unique(c(xc_ko_thin(p$mat_clusters, p$stratum, p$grid),
                     xc_ko_thin(p$cyc_clusters, p$stratum, p$grid)))
    both <- intersect(p$mat_clusters, p$cyc_clusters)
    shared_cat <- intersect(xc_category(p$mat_set), xc_category(p$cyc_set))
    tagList(
      if (length(gone)) red(strong("Dropped: "), paste(unique(gone), collapse = ", "),
        " has no P7 KO-vs-WT table in the ", p$stratum, " stratum, so it is not in the union. ",
        if ("CM4" %in% gone) "CM4 has no G1 stratum for any contrast." else ""),
      if (identical(p$stratum, "all")) red(strong("Sort-confounded. "),
        "P7 was FACS cycling-enriched 4.5–5.2× relative to P0, so the WT P0→P7 side reads ",
        "out the sort as well as development — most of all for the cell-cycle category. ",
        "The G1 stratum holds that fixed, at the cost of thinning the cell-cycle signal ",
        "by construction and dropping CM4."),
      if (length(thin)) red(strong("Thin arms: "), paste(thin, collapse = ", "),
        " have fewer than ", FG$built$thin_cells %||% 50,
        " cells on one side of the P7 contrast. They clear the DE floor but should not ",
        "carry a conclusion on their own."),
      if (length(both)) helpText(strong("Overlapping groups. "),
        paste(both, collapse = ", "), " is in both cluster groups, so the two KO unions ",
        "are not independent — a gene can reach the shared region of more than one ",
        "comparison through it alone."),
      if (length(shared_cat)) helpText(strong("Overlapping categories. "),
        paste(shared_cat, collapse = ", "),
        " belong to both curated categories, so they can appear on both WT sides. ",
        "The in_maturation_set / in_cellcycle_set columns make that visible."),
      if (p$hide_mt) helpText(strong("mt- genes hidden. "),
        "They are up in the KO in every subcluster and down in none, which is a library ",
        "read-fraction difference rather than biology. Untick to include them."))
  })

  # ---- Step 1: WT P0->P7 by category and direction ----
  xc_wt_tab_df <- reactive(xc_wt_df(xc_p()))
  output$xc_wt_tab <- renderDT(DT::datatable(xc_wt_tab_df(), rownames = FALSE,
    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"),
    class = "compact stripe hover") |>
    DT::formatSignif(c("log2FoldChange","padj","pct_P7","pct_P0"), 3))
  xc_wt_p <- reactive({
    d <- xc_wt_tab_df()
    n <- as.data.frame(table(category = d$category, direction = d$direction),
                       stringsAsFactors = FALSE)
    ggplot(n, aes(category, Freq, fill = direction)) +
      geom_col(position = position_dodge(width = .8), width = .7) +
      geom_text(aes(label = Freq), position = position_dodge(width = .8),
                vjust = -0.3, size = 4) +
      scale_fill_manual(values = c("up at P7" = "#c62828", "up at P0" = "#1565c0")) +
      theme_minimal(base_size = input$xc_basesize %||% 13) +
      labs(x = NULL, y = "genes", fill = NULL,
           title = sprintf("WT P0→P7 in %s — genes per category and direction", xc_p()$wt_cluster))
  })
  output$xc_wt_plot <- renderPlot(apply_fig_opts(xc_wt_p(), "xc", input))
  register_fig(output, "xcwt", xc_wt_p, input, opts_prefix = "xc")
  output$xc_wt_note <- renderUI({
    p <- xc_p()
    d <- fg_grid(p$grid)[[p$wt_cluster]][[paste0(XC_WT_KEY, "__", p$stratum)]]
    ct <- fg_ct(XC_WT_KEY); req(ct)
    a <- xc_measure_audit(p)
    # an empty list under one measure and a non-empty one under the other is the
    # difference between "nothing changed" and "the cut could not see it"
    disagree <- !is.null(a) && any(pmin(a$n_AUC, a$n_log2FC) == 0 & pmax(a$n_AUC, a$n_log2FC) > 0)
    tagList(
      de_context_note(p$wt_cluster, ct, p$stratum, p$grid, d),
      div(style = "font-size:13px",
          sprintf("Categories: %s (%d genes) · %s (%d genes), %d in both. Cuts: padj < %s, %s.",
                  xc_cat_label(p$mat_set), length(xc_category(p$mat_set)),
                  xc_cat_label(p$cyc_set), length(xc_category(p$cyc_set)),
                  length(intersect(xc_category(p$mat_set), xc_category(p$cyc_set))),
                  p$padj, xc_eff_label(p))),
      if (!is.null(a)) div(style = "font-size:12px;margin-top:2px",
        HTML(paste0("Genes passing by measure &mdash; ",
          paste(sprintf("<b>%s, %s</b>: AUC %d / log2FC %d", a$category, a$direction,
                        a$n_AUC, a$n_log2FC), collapse = " &middot; ")))),
      if (disagree) div(style = "color:#c62828;font-size:12px;margin-top:2px",
        HTML(paste("&#9888; The two measures disagree on whether a list is empty, so the",
                   "choice in the sidebar decides the answer. presto's log2FC is a difference",
                   "of mean log-normalised expression and therefore scales with expression",
                   "level: from WT P0 to P7 <i>Mcm3</i> quadruples its detection rate",
                   "(6.6% &rarr; 27.2%) for a log2FC of 0.16, while <i>Myh7</i> moves 2.5.",
                   "A symmetric |log2FC| cut can classify maturation genes and can never",
                   "classify a cell-cycle one. AUC is rank-based and means the same thing on",
                   "both sides, which is why build_fourgroup.R classifies its maturation axis",
                   "on AUC too."))))
  })

  # ---- Step 2: P7 KO vs WT across the two cluster groups ----
  xc_ko_tab_df    <- reactive(xc_ko_df(xc_p()))
  xc_ko_pivot_dat <- reactive(xc_ko_pivot_df(xc_p()))
  output$xc_ko_tab <- renderDT(DT::datatable(xc_ko_tab_df(), rownames = FALSE,
    filter = "top", options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"),
    class = "compact stripe hover") |>
    DT::formatSignif(c("log2FoldChange","padj","pct_KO","pct_WT"), 3))
  output$xc_ko_pivot <- renderDT({ d <- xc_ko_pivot_dat()
    DT::datatable(d, rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"),
      class = "compact stripe hover") |>
      DT::formatSignif(setdiff(names(d), c("gene","n_clusters")), 3) })
  output$xc_ko_note <- renderUI({ p <- xc_p()
    div(style = "font-size:13px;margin-bottom:6px",
        HTML(sprintf(paste0("<b>P7: KO vs WT</b> &middot; maturation clusters <b>%s</b>",
                            " &middot; cycling clusters <b>%s</b> &middot; %s cells%s.",
                            " log2FC &gt; 0 is up in KO."),
                     paste(p$mat_clusters, collapse = ", "),
                     paste(p$cyc_clusters, collapse = ", "), p$stratum,
                     if (identical(p$grid, "de2")) " &middot; curated panel" else "")))
  })

  # ---- the four crossings ----
  # vn_plot()/vn_stats()/vn_regions() are reused unchanged: xc_comparison() returns
  # exactly the list(sets, universe) shape they already take.
  xc_venn_p <- function(i) { force(i)
    reactive({ vs <- xc_all()[[i]]
      vn_plot(vs, input$xc_basesize %||% 13, xc_label(vs$cmp, xc_p())) }) }
  xc_vp1 <- xc_venn_p(1); xc_vp2 <- xc_venn_p(2)
  xc_vp3 <- xc_venn_p(3); xc_vp4 <- xc_venn_p(4)
  output$xc_venn1 <- renderPlot(apply_fig_opts(xc_vp1(), "xc", input))
  output$xc_venn2 <- renderPlot(apply_fig_opts(xc_vp2(), "xc", input))
  output$xc_venn3 <- renderPlot(apply_fig_opts(xc_vp3(), "xc", input))
  output$xc_venn4 <- renderPlot(apply_fig_opts(xc_vp4(), "xc", input))
  # spelled out rather than looped: check_download_coverage.R matches figure prefixes
  # by literal string, so a paste0() here would read as "no download registered"
  register_fig(output, "xcvenn1", xc_vp1, input, opts_prefix = "xc")
  register_fig(output, "xcvenn2", xc_vp2, input, opts_prefix = "xc")
  register_fig(output, "xcvenn3", xc_vp3, input, opts_prefix = "xc")
  register_fig(output, "xcvenn4", xc_vp4, input, opts_prefix = "xc")
  output$xc_venn_note <- renderUI({
    p <- xc_p(); vs <- xc_all()
    empt <- vapply(vs, function(v) any(vapply(v$sets, function(s) !length(s$genes), TRUE)), TRUE)
    sz   <- vapply(vs, function(v) length(v$sets[[2]]$genes) /
                     max(length(v$sets[[1]]$genes), 1), 0)
    # The KO side is not symmetric on this data: with AUC 0.60 CM1 has ~1,700 genes
    # down in KO and ~50 up. Seven independently clustered populations do not agree
    # that hard by biology -- it is the same one-directional signature the 2026-08-21
    # deliverable traced to a library read-fraction difference. Worth knowing before
    # reading a big B circle as a big result.
    ko <- do.call(rbind, lapply(unique(c(p$mat_clusters, p$cyc_clusters)), function(cl) {
      u <- xc_ko_long(cl, "up",   p$stratum, p$grid, p$padj, p$eff, p$measure, p$hide_mt)
      d <- xc_ko_long(cl, "down", p$stratum, p$grid, p$padj, p$eff, p$measure, p$hide_mt)
      data.frame(cluster = cl, up = if (is.null(u)) 0L else nrow(u),
                 down = if (is.null(d)) 0L else nrow(d)) }))
    lop <- !is.null(ko) && nrow(ko) &&
      median(pmax(ko$up, ko$down) / pmax(pmin(ko$up, ko$down), 1)) >= 5
    tagList(
      div(style = "font-size:13px;margin-bottom:4px",
          sprintf("Universe: %s genes testable in common. Set sizes are lopsided by design — the WT side is one curated category, the KO side is every DE gene in a cluster group — so read the fold and p in %s, not the picture.",
                  format(length(vs[[1]]$universe), big.mark = ","), "Overlap statistics")),
      if (any(empt)) div(style = "color:#c62828;font-size:13px",
        sprintf("%d of 4 comparisons have an empty circle at these cuts. The counts under each title say which side; Step 1 reports what the other effect-size measure would have given.",
                sum(empt))),
      if (lop) div(style = "color:#c62828;font-size:12px",
        HTML(paste0("&#9888; The KO contrast is strongly one-directional here (",
          paste(sprintf("%s %d&uarr;/%d&darr;", ko$cluster, ko$up, ko$down), collapse = ", "),
          "). Seven independently clustered populations do not agree that hard by biology; ",
          "the 2026-08-21 deliverable traced the same one-directional signature to a library ",
          "read-fraction difference between the two samples. A large KO circle is partly that."))))
  })
  xc_stats_df <- reactive({
    p <- xc_p()
    do.call(rbind, lapply(xc_all(), function(vs)
      cbind(comparison = xc_label(vs$cmp, p), vn_stats(vs), stringsAsFactors = FALSE)))
  })
  output$xc_stats <- renderDT(DT::datatable(xc_stats_df(), rownames = FALSE,
    options = list(pageLength = 10, scrollX = TRUE, dom = "t"),
    class = "compact stripe hover") |>
    DT::formatSignif(c("expected","fold","p_hypergeom"), 3))

  # ---- the gene tables behind each Venn ----
  xc_gene_all <- reactive({ p <- xc_p()
    d <- xc_rbind(lapply(xc_all(), xc_gene_df, p = p))
    validate(need(!is.null(d) && nrow(d),
                  "No genes in any of the four comparisons at these cuts."))
    d })
  observe({
    p <- try(xc_p(), silent = TRUE); if (inherits(p, "try-error")) return()
    ch <- vapply(XC_COMPARISONS, xc_label, "", p = p)
    ch <- c("All four" = "__all__", setNames(ch, ch))
    cur <- isolate(input$xc_gene_cmp)
    updateSelectInput(session, "xc_gene_cmp", choices = ch,
                      selected = if (!is.null(cur) && cur %in% ch) cur else "__all__")
  })
  xc_genes_df <- reactive({ d <- xc_gene_all()
    sel <- input$xc_gene_cmp %||% "__all__"
    if (!identical(sel, "__all__")) d <- d[d$comparison == sel, , drop = FALSE]
    validate(need(nrow(d), "No genes in this comparison at these cuts."))
    d })
  output$xc_genes <- renderDT({ d <- xc_genes_df()
    num <- intersect(c("wt_log2FC","wt_padj","ko_log2FC_top","ko_padj_top",
                       grep("^lfc_", names(d), value = TRUE)), names(d))
    DT::datatable(d, rownames = FALSE, filter = "top",
      options = list(pageLength = 20, scrollX = TRUE, dom = "ftip"),
      class = "compact stripe hover") |> DT::formatSignif(num, 3) })

  output$xc_book <- dl_book(
    basename = function() paste0("wt_programs_x_ko_clusters_",
                                 (xc_p())$stratum,
                                 if (identical((xc_p())$grid, "de2")) "_curated" else ""),
    sheets_fn = function() {
      p <- xc_p(); sh <- list()
      sh[["Step1_WT_P0_vs_P7"]] <- xc_wt_tab_df()
      sh[["Step2_P7_KO_vs_WT"]] <- xc_ko_tab_df()
      sh[["Step2_gene_x_cluster"]] <- xc_ko_pivot_dat()
      sh[["Overlap_stats"]] <- xc_stats_df()
      for (k in seq_along(XC_COMPARISONS)) {
        g <- xc_gene_df(xc_all()[[k]], p)
        # sheet names are capped at 31 chars, so key not label
        sh[[paste0("C", k, "_", XC_COMPARISONS[[k]]$key)]] <-
          if (is.null(g)) data.frame(note = "No genes at these cuts.") else g
      }
      sh
    },
    title = function() { p <- xc_p()
      paste0("WT P0-vs-P7 programs crossed with P7 KO-vs-WT clusters (", p$stratum, " cells)") },
    notes = c(
      "Step 1 is WT_P0_vs_P7: log2FoldChange > 0 means UP AT P7. Step 2 is P7_KO_vs_WT: log2FoldChange > 0 means UP IN KO.",
      "The WT side is filtered to a curated gene category. The KO side is NOT filtered by category - it is every DE gene in the named subclusters, unioned. Filtering both sides by category would intersect the maturation set with the cell-cycle set, which share only Mki67 and Top2a.",
      "Source tables come from build_fourgroup.R and are row-gated: a gene is kept when expressed in >= 5% of one arm AND (padj < 0.05 OR |log2FC| >= 0.5). These are NOT complete gene lists; for every tested gene use the analysis/ pipeline CSVs.",
      "P7 was FACS cycling-enriched 4.5-5.2x relative to P0, so the 'all cells' WT P0-vs-P7 contrast reads out the sort as well as development. The G1 sheets hold cycling composition fixed but thin the cell-cycle category by construction and drop CM4, which has no G1 stratum.",
      "Mitochondrially-encoded (mt-) genes are up in the KO in every subcluster and down in none - a library read-fraction difference, not biology. They are excluded unless the sidebar box was unticked.",
      "n = 1 animal per genotype x timepoint, sexes differ between genotypes, and the KO is not transcript-confirmed. Cell-level Wilcoxon p-values are pseudoreplicated. Descriptive and hypothesis-generating only.")
  )

  # ---- "how this plot was made" blocks (bodies are free functions, see above) ----
  output$mat_score_def <- renderUI({ sc <- input$mat_score %||% "sig_maturation"
    method_note(score_def_ui(sc),
                title = paste0("What ", sc, " is, and which genes are in it")) })
  # the ploidy tab shows four scores at once, so it documents all four
  output$cyc_score_def <- renderUI(method_note(
    tagList(lapply(intersect(c("sig_prolif","sig_cytokinesis","sig_ccexit","sig_ploidy"),
                             if (is.null(SCOREMETA)) character(0) else SCOREMETA$score),
      function(sc) tagList(tags$h6(labof(sc), style = "margin-top:10px;font-weight:700"),
                           score_def_ui(sc)))),
    title = "What these scores are, and which genes are in them"))
  output$mat_violin_method  <- renderUI(mat_violin_method_note(input$mat_score %||% "sig_maturation"))
  output$mat_scatter_method <- renderUI(mat_scatter_method_note())
  output$gm_method          <- renderUI(gm_method_note(input$gm_panel %||% "avg"))
  output$xc_venn_method     <- renderUI({ p <- try(xc_p(), silent = TRUE)
    if (inherits(p, "try-error")) NULL else xc_venn_method_note(p) })

  # ---- Gene-set Venn ----
  observe({
    cc <- fg_cluster_choices()
    if (length(cc)) updateSelectInput(session, "vn_cluster", choices = cc, selected = "AllCM")
  })
  vn_v <- reactive({ m <- input$vn_measure %||% "lfc"
    vn_sets(c(input$vn_a, input$vn_b, input$vn_c), input$vn_cluster %||% "AllCM",
            input$vn_stratum %||% "G1", input$vn_grid %||% "de", input$vn_padj %||% 0.05,
            if (identical(m, "auc")) input$vn_auc %||% 0.60 else input$vn_lfc %||% 0.25, m) })
  vn_p <- reactive(vn_plot(vn_v(), input$vn_basesize %||% 13))
  output$vn_plot <- renderPlot(apply_fig_opts(vn_p(), "vn", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("vn_dl_", f)]] <- dl_ggplot("vn", vn_p, input, f) })
  output$vn_note <- renderUI({
    v <- try(vn_v(), silent = TRUE)
    if (inherits(v, "try-error")) return(NULL)
    empt <- vapply(v$sets, function(x) length(x$genes) == 0L, TRUE)
    tagList(
      div(style = "font-size:13px;margin-bottom:4px",
          sprintf("Universe: %s genes testable in common. %s",
                  format(length(v$universe), big.mark = ","),
                  paste(vapply(v$sets, function(x) sprintf("%s: %d", x$label, length(x$genes)), ""),
                        collapse = " · "))),
      if (any(empt)) div(style = "color:#c62828;font-size:13px",
        sprintf("%s empty for this cluster/stratum — try the other DE matrix or loosen the cuts.",
                paste(vapply(v$sets[empt], `[[`, "", "label"), collapse = ", "))))
  })
  # the two caveats that must not be left to the reader to reconstruct
  output$vn_caveat <- renderUI({
    ids <- c(input$vn_a, input$vn_b, input$vn_c)
    tagList(
      if (any(grepl("^ax:cyc:", ids))) helpText(style = "color:#c62828",
        strong("Data-driven cycling set: "), "531 genes, of which only 46 are canonical. ",
        "The rest include Ran, Nap1l1, Calm1 and Ppia — cycling cells are globally more ",
        "transcriptionally active, so the axis partly measures output, not cell cycle. ",
        "The curated canonical set is the safer circle to label “cycling genes”."),
      if (sum(grepl("^de:(WT|KO)_P0_vs_P7", ids)) >= 2) helpText(style = "color:#c62828",
        strong("WT-only vs KO-only: "), "the KO contributes more P7 cells than the WT ",
        "(3,674 vs 2,496 in G1 pooled over cardiomyocytes), so a larger KO-only region is ",
        "partly a power difference rather than biology."))
  })
  output$vn_stats <- renderDT(DT::datatable(vn_stats(vn_v()), rownames = FALSE,
    options = list(pageLength = 10, scrollX = TRUE, dom = "t"),
    class = "compact stripe hover") |>
    DT::formatSignif(c("expected","fold","p_hypergeom"), 3))
  output$vn_regions <- renderDT(DT::datatable(vn_region_df(vn_v()), rownames = FALSE,
    options = list(pageLength = 10, scrollX = TRUE, dom = "ftip",
                   columnDefs = list(list(targets = 3, render = DT::JS(
                     "function(d,t,r,m){return t==='display'&&d&&d.length>140?",
                     "'<span title=\"'+d+'\">'+d.substr(0,140)+'…</span>':d;}")))),
    class = "compact stripe hover", escape = FALSE))

  # ---- Maturation intersection + candidate genes ----
  mi_quad_p <- reactive(fg_quadrant_plot(input$mi_cluster, input$mi_hideconf,
                                         bs = input$miquad_basesize %||% 13))
  output$mi_quadrant <- renderPlot(apply_fig_opts(mi_quad_p(), "miquad", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("miquad_dl_", f)]] <- dl_ggplot("miquad", mi_quad_p, input, f) })
  mi_tab_df <- reactive(fg_intersect_df(input$mi_cluster, input$mi_quad, input$mi_hideconf))
  output$mi_cyc_note <- renderUI({
    d <- try(fg_intersect_df(input$mi_cluster, NULL, input$mi_hideconf), silent = TRUE)
    if (inherits(d, "try-error") || !"cyc_resid" %in% names(d)) return(NULL)
    h <- d[d$quadrant %in% c("immature_up_in_KO","mature_down_in_KO"), , drop = FALSE]
    if (!nrow(h)) return(NULL)
    nc <- sum(!is.na(h$cyc_class) & h$cyc_class == "cycling-associated")
    rv <- h$cyc_resid[!is.na(h$cyc_resid)]
    div(style = "font-size:13px;margin:6px 0;padding:8px 10px;background:#f6f8fa;border-left:4px solid #90a4ae",
      HTML(sprintf(paste0(
        "<b>Do these genes link maturation to cycling?</b> Of the %d genes in the two ",
        "hypothesis quadrants here, <b>%d</b> are cycling-associated. Median residual ",
        "cycling association is <b>%+.3f</b> — the part of a gene's cycling link not already ",
        "explained by where it sits on the maturation axis. Zero would mean “exactly as ",
        "expected”; %s."),
        nrow(h), nc, if (length(rv)) stats::median(rv) else NA_real_,
        if (length(rv) && stats::median(rv) < 0)
          "negative means these genes are <i>less</i> cycling-linked than their maturation position predicts"
        else "positive would mean more")))
  })
  output$mi_table <- renderDT(DT::datatable(mi_tab_df(), rownames = FALSE,
    options = list(pageLength = 25, scrollX = TRUE, scrollY = "420px",
                   scrollCollapse = TRUE, dom = "ftip"),
    class = "compact stripe hover") |>
    DT::formatSignif(intersect(c("mat_log2FC","mat_auc","cyc_auc","cyc_resid",
                                 "p7ko_log2FC","p7ko_padj"), names(mi_tab_df())), 3))

  mi_cand_p <- reactive(fg_candidate_plot(input$mi_genes, input$mi_cand_clusters,
                                          input$mi_cand_stratum %||% "all",
                                          input$micand_basesize %||% 13))
  output$mi_candidates <- renderPlot(apply_fig_opts(mi_cand_p(), "micand", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("micand_dl_", f)]] <- dl_ggplot("micand", mi_cand_p, input, f) })
  mi_spec_df <- reactive(fg_specificity_df(input$mi_genes, input$mi_cand_stratum %||% "all",
                                           input$mi_spec_grid %||% "de"))
  output$mi_spec_tab <- renderDT(DT::datatable(mi_spec_df(), rownames = FALSE,
    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"),
    class = "compact stripe hover") |>
    DT::formatSignif(intersect(c("P7_KO_vs_WT","P0_KO_vs_WT","P7_specificity",
                                 "CM2_4_5_mean_absLFC","other_clusters_mean_absLFC",
                                 "priority_concentration","mat_auc","cyc_auc"),
                               names(mi_spec_df())), 3))
  register_dl(output, "mi_cand", function()
    fg_candidate_df(input$mi_genes, input$mi_cand_clusters, input$mi_cand_stratum %||% "all"),
    "candidate_expression")

  # ---- Subset & DEGs (interactive descriptive DE) ----
  observeEvent(input$deg_by, {
    req(input$deg_by)   # ignoreNULL=FALSE + no UI round-trip (testServer) => NULL reaches here
    lv <- sort(unique(as.character(DMETA[[input$deg_by]])))
    defA <- defB <- NULL
    if (input$deg_by == "genotype")      { defA <- "KO";  defB <- "WT" }
    else if (input$deg_by == "timepoint"){ defA <- "P0";  defB <- "P7" }
    else if (input$deg_by == "cycling")  { defA <- "TRUE"; defB <- "FALSE" }
    else if (input$deg_by == "Phase")    { defA <- intersect("G2M", lv); defB <- intersect("G1", lv) }
    defA <- intersect(defA, lv); defB <- intersect(defB, lv)
    updateSelectizeInput(session, "deg_a", choices = lv, selected = defA)
    updateSelectizeInput(session, "deg_b", choices = lv, selected = defB)
  }, ignoreNULL = FALSE)

  deg_lab <- function(x) if (length(x)) paste(x, collapse = "/") else "?"
  deg_res <- eventReactive(input$deg_run, {
    req(input$deg_by, input$deg_a, input$deg_b)
    validate(need(!length(intersect(input$deg_a, input$deg_b)),
                  "Groups A and B overlap — pick distinct levels."))
    filters <- setNames(lapply(CAT_COLS, function(c) input[[paste0("degf_", c)]]), CAT_COLS)
    deg_compute(deg_mask(filters), input$deg_by, input$deg_a, input$deg_b)
  }, ignoreNULL = FALSE)
  deg_pick <- reactiveVal(NULL)
  observeEvent(input$deg_run, deg_pick(NULL))                 # reset selection on a new run
  deg_tab <- reactive(de_table(drop_conf(deg_res(), input$deg_hideconf)))   # full table
  deg_dt_proxy <- DT::dataTableProxy("deg_table")
  output$deg_n <- renderText({
    d <- deg_res()
    sprintf("Group A (%s): %d cells   |   Group B (%s): %d cells   |   %d genes tested",
            deg_lab(isolate(input$deg_a)), d$n_A[1], deg_lab(isolate(input$deg_b)), d$n_B[1], nrow(d))
  })
  output$deg_volcano <- renderPlotly({
    d <- drop_conf(deg_res(), input$deg_hideconf)
    de_volcano_ly(d, paste0(deg_lab(isolate(input$deg_a)), "  vs  ", deg_lab(isolate(input$deg_b))),
                  "deg_volcano", pos = paste0("up in ", deg_lab(isolate(input$deg_a))),
                  neg = paste0("up in ", deg_lab(isolate(input$deg_b))),
                  xlab = "logFC  (A / B)", highlight = deg_pick(),
                  lfc_cut = input$deg_vlfc %||% DE_LFC_CUT)
  })
  outputOptions(output, "deg_volcano", suspendWhenHidden = FALSE)  # render at startup so plotly_click source registers before its click observer fires
  # static twin for vector export + studio (same groups, labels and colour cut)
  deg_volcano_gg <- reactive({
    d <- drop_conf(deg_res(), input$deg_hideconf)
    de_volcano(d, paste0(deg_lab(isolate(input$deg_a)), "  vs  ", deg_lab(isolate(input$deg_b))),
               pos = paste0("up in ", deg_lab(isolate(input$deg_a))),
               neg = paste0("up in ", deg_lab(isolate(input$deg_b))),
               xlab = "logFC  (A / B)", lfc_cut = input$deg_vlfc %||% DE_LFC_CUT,
               highlight = deg_pick()) })
  register_fig(output, "degvolc", deg_volcano_gg, input)
  output$deg_table   <- renderDT(de_datatable(deg_tab()))
  output$deg_pick_ui <- renderUI(pick_banner(deg_pick(), "deg_clear"))
  output$deg_geneinfo <- renderUI(gene_info_card(deg_pick(), "deg_infoclose"))
  observeEvent(event_data("plotly_click", source = "deg_volcano"),
    deg_pick(event_data("plotly_click", source = "deg_volcano")$customdata))
  observeEvent(input$deg_table_rows_selected, {                          # table row -> selected gene
    r <- input$deg_table_rows_selected; g <- if (length(r)) deg_tab()$gene[r] else NULL
    if (!identical(g, deg_pick())) deg_pick(g)
  }, ignoreNULL = FALSE)
  observeEvent(deg_pick(), {                                             # selected gene -> highlight table row
    g <- deg_pick(); rows <- if (!is.null(g)) which(deg_tab()$gene == g) else integer(0)
    if (!identical(as.integer(rows), as.integer(input$deg_table_rows_selected)))
      DT::selectRows(deg_dt_proxy, if (length(rows)) rows else NULL)
  }, ignoreNULL = FALSE)
  observeEvent(input$deg_clear, deg_pick(NULL))
  observeEvent(input$deg_infoclose, deg_pick(NULL))

  # ---- Pathways & enrichment (precomputed) ----
  if (!is.null(ENR)) {
    cts <- enr_celltypes()
    updateSelectInput(session, "enr_ct", choices = cts,
                      selected = if ("Cardiomyocyte" %in% cts) "Cardiomyocyte" else cts[1])
  }
  output$enr_gsea_plot <- renderPlotly({ req(input$enr_ct); enr_gsea_plot(input$enr_ct, input$enr_tp,
                            base_size = input$enrgsea_basesize %||% 12,
                            pal_choice = input$enrgsea_palette,
                            label_chars = input$enrgsea_labelchars %||% 46,
                            axis_scale = input$enrgsea_axisscale %||% 1) })
  output$enr_gsea_tab  <- renderDT({ req(input$enr_ct); enr_gsea_table(input$enr_ct, input$enr_tp) })
  output$enr_go_plot   <- renderPlotly({ req(input$enr_ct); enr_go_plot(input$enr_ct, input$enr_tp,
                            base_size = input$enrgo_basesize %||% 11,
                            pal_choice = input$enrgo_palette,
                            label_chars = input$enrgo_labelchars %||% 46,
                            axis_scale = input$enrgo_axisscale %||% 1,
                            colour_by = input$enrgo_colourby %||% "padj") })
  output$enr_go_tab    <- renderDT({ req(input$enr_ct); enr_go_table(input$enr_ct, input$enr_tp) })
  output$enr_e2f_heat  <- renderPlotly(enr_e2f_heat())
  output$enr_tf_top    <- renderPlotly({ req(input$enr_ct); enr_tf_top(input$enr_ct,
                            base_size = input$enrtf_basesize %||% 11,
                            pal_choice = input$enrtf_palette,
                            axis_scale = input$enrtf_axisscale %||% 1) })
  register_fig(output, "enre2f", reactive(enr_e2f_heat_gg()), input)
  register_fig(output, "enrtf",  reactive({ req(input$enr_ct)
    enr_tf_top_gg(input$enr_ct,
                  base_size = input$enrtf_basesize %||% 11,
                  pal_choice = input$enrtf_palette,
                  axis_scale = input$enrtf_axisscale %||% 1) }), input)
  register_fig(output, "enrgsea", reactive({ req(input$enr_ct)
    gsea_barplot_gg(enr_gsea(input$enr_ct, input$enr_tp),
                    paste0("GSEA — ", input$enr_ct, " ", input$enr_tp),
                    base_size = input$enrgsea_basesize %||% 12,
                    pal_choice = input$enrgsea_palette,
                    label_chars = input$enrgsea_labelchars %||% 46,
                    axis_scale = input$enrgsea_axisscale %||% 1) }), input)
  register_fig(output, "enrgo", reactive({ req(input$enr_ct)
    go_dotplot_gg(enr_go(input$enr_ct, input$enr_tp),
                  paste0("GO BP enriched in KO-up genes — ", input$enr_ct, " ", input$enr_tp),
                  base_size = input$enrgo_basesize %||% 11,
                  pal_choice = input$enrgo_palette,
                  label_chars = input$enrgo_labelchars %||% 46,
                  axis_scale = input$enrgo_axisscale %||% 1,
                  colour_by = input$enrgo_colourby %||% "padj") }), input)

  # ---- QC & normalization (embedded figures) ----
  figcard <- function(uri, title, desc) {
    if (is.null(uri) || is.na(uri)) return(NULL)
    div(style = "margin-bottom:22px",
        tags$h5(title, style = "margin-bottom:4px"),
        tags$img(src = uri, style = "max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px"),
        tags$p(desc, style = "color:#444;font-size:13px;margin-top:4px"))
  }
  output$qcfigs <- renderUI({
    validate(need(!is.null(figs), "Figures not available in this data build."))
    tagList(
      tags$p(tags$b("How the raw counts become analysis-ready data — the QC & normalization steps, with plots."),
             style = "font-size:14px"),
      figcard(figs$filtering, "1. QC filtering — cells removed by reason (per sample/lane)",
        "Droplets removed by each quality filter: likely doublets (red), high-mitochondrial cells >20% (orange), and the few trimmed by the upper gene-count cap (blue). Most cells pass."),
      figcard(figs$qc_violins, "1. Per-cell QC distributions after filtering",
        "Genes per cell, sequencing depth (UMIs, log scale), and % mitochondrial reads for each sample. Dashed lines mark the thresholds (genes >= 1500; mito <= 20%)."),
      figcard(figs$doublet, "2. Doublet detection — scDblFinder vs Scrublet",
        "Fraction of cells flagged as two-cell droplets by each method, per lane. scDblFinder recovers ~5-8% (near the ~8% expected) where the original Scrublet calls (~1-2%) under-called; doublets are removed before analysis."),
      figcard(figs$hvg, "3-4. Normalization (SCTransform) & feature selection",
        "After variance-stabilizing normalization, each gene's variability vs mean expression; the ~2,000 highly-variable genes (red) carry the biological signal used for clustering."),
      figcard(figs$harmony, "5. Batch integration (Harmony) — before vs after",
        "UMAP coloured by library before and after Harmony integration. Before, cells split by sample (a technical batch effect); after, libraries intermix while biological structure is preserved, so clusters reflect cell type.")
    )
  })
  score_meta_df <- reactive({
    validate(need(!is.null(SCOREMETA),
      "Module-score definitions aren't in this data build — run build_signature_scores.R and redeploy."))
    SCOREMETA
  })
  output$score_meta_tab <- renderDT(enr_dt(score_meta_df()))
  output$score_sets_tab <- renderDT(enr_dt(score_sets_df()))
  output$doublet_tab <- renderTable({ validate(need(!is.null(tabs$doublet), "No doublet table.")); tabs$doublet },
                                    striped = TRUE, hover = TRUE)

  output$about <- renderUI(HTML(paste0(
    "<h3>E2F7/8 knockout mouse-heart single-cell RNA-seq</h3>",
    "<p>Single-cell RNA-seq of E2F7/8 knockout (KO) vs wild-type (WT) developing mouse heart at P0 and P7. ",
    "Explore ~", format(nrow(meta), big.mark = ","), " cells: colour the map by any gene or metadata, compare KO vs WT, ",
    "inspect differential expression by cell type and by cardiomyocyte subcluster, and view subcluster identity / cell-cycle.</p>",
    "<h4>How the data were processed (normalization &amp; preprocessing)</h4>",
    "<ol style='font-size:13.5px'>",
    "<li><b>QC filtering</b> &mdash; per lane, keep cells with &ge;1,500 genes and &le;20% mitochondrial reads (mouse <code>mt-</code>), with an upper gene-count cap to drop likely multiplets. The mito cap is generous on purpose: heart cells are genuinely mitochondria-rich.</li>",
    "<li><b>Doublet removal</b> &mdash; two-cell droplets detected (Scrublet per lane, cross-checked with scDblFinder) and removed.</li>",
    "<li><b>Normalization</b> &mdash; SCTransform (glmGamPoi) variance-stabilizing transform so cells of different sequencing depth are comparable; mitochondrial percent deliberately <i>not</i> regressed out.</li>",
    "<li><b>Feature selection</b> &mdash; ~2,000 most-variable genes drive dimensionality reduction.</li>",
    "<li><b>Integration &amp; embedding</b> &mdash; the four libraries integrated with Harmony, then clustered and laid out as the UMAP shown here. Differential expression is computed on the raw / pseudobulk counts, not the normalized/integrated values.</li>",
    "</ol>",
    "<p style='font-size:12px;color:#777'>Full method detail and the comparison to the original analysis are in the accompanying report and <code>00_DOCS/NORMALIZATION.md</code>.</p>",
    # The methods book (docs/) is the long form of every tab in this app: parameters,
    # rationale, and what each result cannot support. Kept as a repo-relative pointer
    # rather than a URL because the book is a build artifact, not a deployed site.
    "<p style='font-size:12px;color:#777'>Chapter-by-chapter methods for every tab of this app, ",
    "including the parameters and the reasoning behind each choice, are in the repository's ",
    "methods book: <code>docs/</code> (build with <code>docs/render.sh</code>).</p>",
    "<div style='background:#fff3e0;border-left:4px solid #e65100;padding:10px 14px;border-radius:4px'>",
    "<b>Critical caveats — read before interpreting:</b><ul>",
    "<li><b>n = 1 animal per condition.</b> Two lanes per sample are the same library sequenced twice, not biological replicates.</li>",
    "<li><b>Sex confound:</b> KO and WT animals are different sexes (Y-genes top the KO-up list; flagged as 'sex/construct').</li>",
    "<li><b>KO not confirmed:</b> E2f7/E2f8 are not reduced at the transcript level (likely a conditional allele a 3' assay cannot see).</li>",
    "<li><b>Cycling fractions are a sorting artefact.</b> P7 was FACS cycling-enriched 4.5&ndash;5.2&times; while P0 was ",
    "essentially unenriched, so the raw cycling fraction <i>rises</i> from P0 to P7 in this data and <i>falls</i> in reality. ",
    "Any P0-vs-P7 contrast on all cells partly reads out the sort rather than development &mdash; prefer the ",
    "phase-matched (G1-only) contrasts on the <b>Four-group</b> tab. Only within-timepoint comparisons ",
    "(P7 KO vs WT) are free of this.</li>",
    "<li>All KO-vs-WT differences are <b>descriptive / hypothesis-generating only.</b> ",
    "Valid inference needs a replicated, sex-matched cohort (&ge;3 animals/condition).</li>",
    "<li><b>Volcano p-axis is for ranking, not significance.</b> Because the n = 1 design treats technical ",
    "replicates as biological ones (pseudoreplication), the test under-estimates variance and returns ",
    "extreme p-values &mdash; many so small they underflow to 0 in double precision. The plot floors these at ",
    "1e-300, so points pinned at <code>-log10 p = 300</code> simply mean &ldquo;smaller than the computer can ",
    "represent,&rdquo; not a real significance level. Use the vertical axis only to rank candidate genes.</li></ul></div>",
    "<p style='color:#777;font-size:12px'>Live views colour a stratified sample of ", format(nrow(meta), big.mark = ","),
    " cells over a ", length(genes), "-gene curated panel at full resolution. Any other gene seen in a volcano ",
    "can still be plotted from a broader ~", format(if (!is.null(EXPR)) ncol(EXPR) else 0L, big.mark = ","),
    "-cell subset (the UMAP title flags when this fallback is in use). Differential-expression tables and heatmaps ",
    "are computed from the FULL data. Built ", app$built, ".</p>")))

  # ---- Cell-cycle exit & ploidy (module scores from build_signature_scores.R) ----
  cyc_violins_plot <- reactive({
    bs <- input$cyc_basesize %||% 13; ch <- input$cyc_palette %||% "Default"; ct <- input$cyc_ct
    cols <- intersect(c("sig_prolif","sig_cytokinesis","sig_ccexit","sig_ploidy"), names(meta))
    validate(need(length(cols), "Run build_signature_scores.R for these scores, then redeploy."))
    df <- meta; if (ct != "All" && has("celltype")) df <- df[as.character(df$celltype) == ct, ]
    has_tp <- "timepoint" %in% names(df)
    long <- do.call(rbind, lapply(cols, function(c) {
      d <- df[!is.na(df[[c]]), , drop = FALSE]; if (!nrow(d)) return(NULL)
      data.frame(genotype = d$genotype, timepoint = if (has_tp) as.character(d$timepoint) else "all",
                 score = labof(c), value = d[[c]], stringsAsFactors = FALSE)
    }))
    validate(need(!is.null(long) && nrow(long), "No scored cells for this selection."))
    long$score <- factor(long$score, levels = labof(cols))
    p <- ggplot(long, aes(genotype, value, fill = timepoint)) +
      geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2, position = position_dodge(.9)) +
      facet_wrap(~ score, scales = "free_y") + theme_minimal(base_size = bs) +
      labs(x = NULL, y = "module score", fill = "timepoint", title = paste0("Cycle-exit / ploidy scores — ", ct))
    if (ch != "Default") p <- p + scale_fill_manual(values = disc_pal(levels(factor(long$timepoint)), ch))
    p
  })
  output$cyc_violins <- renderPlot(apply_fig_opts(cyc_violins_plot(), "cyc", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("cyc_dl_", f)]] <- dl_ggplot("cyc", cyc_violins_plot, input, f) })
  cyc_scatter_p <- reactive(ploidy_scatter(input$cyc_ct, input$cyc_basesize %||% 13))
  output$cyc_scatter <- renderPlot(apply_fig_opts(cyc_scatter_p(), "cyc", input))
  register_fig(output, "cycsc", cyc_scatter_p, input, opts_prefix = "cyc")

  # ---- Maturation & metabolism ----
  mat_violin_plot <- reactive(score_violin(input$mat_score, input$mat_ct,
    input$mat_basesize %||% 13, input$mat_palette %||% "Default", input$mat_stratum %||% "all"))
  output$mat_violin <- renderPlot(apply_fig_opts(mat_violin_plot(), "mat", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("mat_dl_", f)]] <- dl_ggplot("mat", mat_violin_plot, input, f) })
  # the cell scatter previously had no export path at all — it is the one plot in the
  # app that could not be downloaded; now wired like every other ggplot panel
  mat_scatter_plot <- reactive(mat_scatter(input$mat_ct, input$matsc_basesize %||% 13,
    input$mat_stratum %||% "all", input$mat_showcells %||% TRUE))
  output$mat_scatter <- renderPlot(apply_fig_opts(mat_scatter_plot(), "matsc", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("matsc_dl_", f)]] <- dl_ggplot("matsc", mat_scatter_plot, input, f) })

  # ---- Gene map (maturation axis x metabolic axis) ----
  observe(updateSelectizeInput(session, "gm_gene",
    choices = c("", if (!is.null(GM)) sort(GM$gene) else character(0)), server = TRUE))
  gm_d    <- reactive(gm_df(input$gm_panel %||% "avg", input$gm_quad, input$gm_dist %||% 0,
                            input$gm_hidesets %||% TRUE, input$gm_geneset %||% "__all__",
                            input$gm_conf %||% FALSE))
  gm_tab  <- reactive(gm_table(gm_d()))
  gm_pick <- reactiveVal(NULL)
  gm_dt_proxy <- DT::dataTableProxy("gm_table")
  output$gm_scatter <- renderPlotly(gm_plot_ly(gm_d(), input$gm_labeln %||% 20, gm_pick(),
                                                input$gm_panel %||% "avg"))
  outputOptions(output, "gm_scatter", suspendWhenHidden = FALSE)  # register the click source at startup
  # static twin for vector export + studio (same panel, filters and labels)
  gm_scatter_gg <- reactive(gm_plot_gg(gm_d(), input$gm_labeln %||% 20, gm_pick(),
                                       input$gm_panel %||% "avg"))
  register_fig(output, "gmsc", gm_scatter_gg, input)
  output$gm_note <- renderUI({
    d <- try(gm_d(), silent = TRUE)
    n <- if (inherits(d, "try-error")) 0 else nrow(d)
    tagList(
      div(style = "font-size:13px;margin-bottom:4px",
          sprintf("%d genes shown of %d on the map — %s. %s confidently labelled.",
                  n, if (is.null(GM)) 0L else nrow(GM),
                  if (identical(input$gm_panel %||% "avg", "avg")) "P0 and P7 averaged"
                  else paste0(input$gm_panel, " cells only"),
                  if (inherits(d, "try-error")) "?" else sum(d$confidence != "neither")),
          if (isTRUE(input$gm_hidesets)) " Scoring-set genes hidden." else
            span(style = "color:#c62828", " Scoring-set genes shown — those sit at their own axis's extreme by construction.")),
      gm_gene_note(gm_pick(), input$gm_panel %||% "avg"))
  })
  output$gm_table    <- renderDT(de_datatable(gm_tab(), scroll = "340px"))
  output$gm_pick_ui  <- renderUI(pick_banner(gm_pick(), "gm_clear"))
  output$gm_geneinfo <- renderUI(gene_info_card(gm_pick(), "gm_infoclose"))
  observeEvent(event_data("plotly_click", source = "gm_scatter"),
    gm_pick(event_data("plotly_click", source = "gm_scatter")$customdata))
  observeEvent(input$gm_gene, if (nzchar(input$gm_gene %||% "")) gm_pick(input$gm_gene))
  observeEvent(input$gm_table_rows_selected, {
    r <- input$gm_table_rows_selected; g <- if (length(r)) gm_tab()$gene[r] else NULL
    if (!identical(g, gm_pick())) gm_pick(g)
  }, ignoreNULL = FALSE)
  observeEvent(gm_pick(), {
    g <- gm_pick(); rows <- if (!is.null(g)) which(gm_tab()$gene == g) else integer(0)
    if (!identical(as.integer(rows), as.integer(input$gm_table_rows_selected)))
      DT::selectRows(gm_dt_proxy, if (length(rows)) rows else NULL)
  }, ignoreNULL = FALSE)
  observeEvent(input$gm_clear, gm_pick(NULL))
  observeEvent(input$gm_infoclose, gm_pick(NULL))

  # ---- Cell-cell signalling (curated L-R; build_communication.R) ----
  output$cc_heat <- renderPlotly(commun_heat(input$cc_pathway, input$cc_tp, input$cc_metric))
  register_fig(output, "ccheat",
    reactive(commun_heat_gg(input$cc_pathway, input$cc_tp, input$cc_metric)), input)
  output$cc_tab  <- renderDT(commun_table(input$cc_pathway, input$cc_tp))

  # ---- Annotation check (reference-marker concordance; build_refmap.R) ----
  output$ann_heat <- renderPlotly(refmap_heat())
  register_fig(output, "annheat", reactive(refmap_heat_gg()), input)
  output$ann_tab  <- renderDT(refmap_table())

  # ---- table downloads, in one place ---------------------------------------
  # Every DT in the app appears here, so "can I get this as a file?" has one
  # answer and adding a table means adding a row. The coverage test in
  # tools/check_download_coverage.R fails the build if a table is missing.
  #
  # `df` returns the FULL frame with the user's filters applied but no display
  # column-subsetting -- someone taking data away wants the columns, not the
  # eight that fitted on screen. `base` may be a function when the filename
  # should carry the current selection.
  # ---- Precomputed results (app$linked, from build_linked_results.R) ----
  observe({
    validate(need(!is.null(LKM), LK_MSG))
    updateSelectInput(session, "lk_group", choices = unique(LKM$group))
  })
  observeEvent(input$lk_group, {
    req(LKM, input$lk_group)
    m <- LKM[LKM$group == input$lk_group, , drop = FALSE]
    updateSelectInput(session, "lk_table", choices = setNames(m$key, m$label))
  }, ignoreNULL = TRUE)
  lk_row <- reactive({
    validate(need(!is.null(LKM), LK_MSG)); req(input$lk_table)
    r <- LKM[LKM$key == input$lk_table, , drop = FALSE]
    validate(need(nrow(r), "Pick a result.")); as.list(r[1, ])
  })
  lk_df <- reactive({
    validate(need(!is.null(LK), LK_MSG)); req(input$lk_table)
    d <- LK[[input$lk_table]]
    validate(need(!is.null(d) && nrow(d), "That result is not in this data build."))
    d
  })
  output$lk_label  <- renderText(lk_row()$label)
  output$lk_note   <- renderUI(div(style = "font-size:13px;color:#555;margin-bottom:8px",
                                   lk_row()$note))
  output$lk_source <- renderUI({ r <- try(lk_row(), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    helpText(style = "font-size:11px",
             sprintf("%s rows · source: our_analysis/results/%s",
                     format(r$rows, big.mark = ","), r$source)) })
  output$lk_tab <- renderDT(enr_dt(lk_df(), scroll = "560px"))

  TABLE_DL <- list(
    # CM object-mode diagnostic
    list(id = "objtest_tab", base = "cm_objectmode_comparison", df = function() objtest_df()),
    # PC-dimension sweep
    list(id = "pcd_tab", base = function() paste0("pcdims_", input$pcd_obj %||% "cm", "_comparison"),
         df = function() pcd_df()),
    # Gene-set provenance
    list(id = "gsp_tab", base = "gene_set_provenance", df = function() gsp_df()),
    # Clustering variants
    list(id = "clu_mk", base = function() paste0("markers_", input$clu_var %||% "variant"),
         df = function() clu_mk_df()),
    list(id = "clu_de", base = function() paste0("DE_", input$clu_var %||% "variant", "_", input$clu_cl %||% ""),
         df = function() clu_de_df()),
    list(id = "gsp_bench", base = "gene_set_benchmark", df = function() gsp_bench_df()),
    # Precomputed upstream results
    list(id = "lk_tab", base = function() paste0("linked_", input$lk_table %||% "result"),
         df = function() lk_df()),
    # Differential expression
    list(id = "ct_table", base = function() paste0("DE_", input$ct_tp, "_", input$ct_sel),
         df = function() drop_conf(ct_d(), input$ct_hideconf)),
    list(id = "deg_table", base = "DEG_subset",
         df = function() drop_conf(deg_res(), input$deg_hideconf)),
    # Cardiomyocyte deep-dive
    list(id = "cm_detab", base = function() { ct <- cm_ct()
           if (is.null(ct)) paste0("DE_KOvsWT_pooled_", input$cm_sub)
           else paste0("DE_", ct$key, "_", input$cm_sub, "_", input$cm_stratum %||% "G1",
                       if (identical(input$cm_grid, "de2")) "_curated" else "") },
         df = function() drop_conf(cm_d(), input$cm_hideconf)),
    list(id = "cm_summary", base = "cm_subcluster_summary_res0.2", df = function() cm_summary_df()),
    list(id = "cm_topmarkers", base = "cm_top_markers_res0.2", df = function() cm_topmarkers_df()),
    list(id = "cm_sub_idgo_tab", base = function() paste0("identity_GO_", input$cm_enr_sub),
         df = function() cm_sub_idgo_df()),
    list(id = "cm_sub_kogo_tab", base = function() cm_enr_dl_base("GO", "_up"),
         df = function() cm_sub_kogo_df()),
    list(id = "cm_sub_kodn_tab", base = function() cm_enr_dl_base("GO", "_down"),
         df = function() cm_sub_kodn_df()),
    list(id = "cm_sub_gsea_tab", base = function() cm_enr_dl_base("GSEA", ""),
         df = function() cm_sub_gsea_df()),
    # Four-group
    list(id = "fg_counts_tab", base = function() paste0("fourgroup_counts_res", FG$built$res),
         df = function() fg_counts_wide()),
    list(id = "fg_detab", base = function() paste0("fourgroup_DE_", input$fg_cluster, "_",
                                                   input$fg_contrast, "_", input$fg_stratum,
                                                   if (identical(input$fg_grid, "de2")) "_curated" else ""),
         df = function() fg_d()),
    list(id = "fg_summary_tab", base = function() paste0("g1_maturation_summary_",
                                                         input$fg_score_stratum %||% "all"),
         df = function() fg_summary_df(input$fg_score %||% "sig_maturation_nocc",
                                       input$fg_score_stratum %||% "all")),
    # Four-group enrichment (build_fourgroup_enrichment.R)
    list(id = "fg_enr_up_tab", base = function() paste0("GO_", input$fg_enr_ont %||% "BP", "_",
                                                        input$fg_enr_contrast, "_", input$fg_enr_cluster, "_up"),
         df = function() fg_enr_up_df()),
    list(id = "fg_enr_dn_tab", base = function() paste0("GO_", input$fg_enr_ont %||% "BP", "_",
                                                        input$fg_enr_contrast, "_", input$fg_enr_cluster, "_down"),
         df = function() fg_enr_dn_df()),
    list(id = "fg_enr_gsea_tab", base = function() paste0("GSEA_", input$fg_enr_contrast, "_",
                                                          input$fg_enr_cluster, "_", input$fg_enr_stratum),
         df = function() fg_enr_gsea_dat()),
    list(id = "fg_enr_audit_tab", base = function() paste0("enrichment_audit_", input$fg_enr_contrast),
         df = function() fg_enr_audit_dat()),
    # Maturation n P7 KO
    list(id = "mi_table", base = function() paste0("maturation_intersect_", input$mi_cluster),
         df = function() mi_tab_df()),
    list(id = "mi_spec_tab", base = "candidate_P7_specificity", df = function() mi_spec_df()),
    # WT programs x KO clusters
    list(id = "xc_wt_tab", base = function() paste0("WTp0p7_by_category_",
           (xc_p())$wt_cluster, "_", (xc_p())$stratum),
         df = function() xc_wt_tab_df()),
    list(id = "xc_ko_tab", base = function() paste0("P7KOvsWT_by_cluster_", (xc_p())$stratum),
         df = function() xc_ko_tab_df()),
    list(id = "xc_ko_pivot", base = function() paste0("P7KOvsWT_gene_x_cluster_", (xc_p())$stratum),
         df = function() xc_ko_pivot_dat()),
    list(id = "xc_stats", base = function() paste0("crossings_overlap_stats_", (xc_p())$stratum),
         df = function() xc_stats_df()),
    list(id = "xc_genes", base = function() paste0("crossings_genes_", (xc_p())$stratum),
         df = function() xc_genes_df()),
    # Gene-set Venn
    list(id = "vn_stats", base = function() paste0("venn_overlap_stats_", input$vn_cluster %||% "AllCM"),
         df = function() vn_stats(vn_v())),
    list(id = "vn_regions", base = function() paste0("venn_regions_", input$vn_cluster %||% "AllCM"),
         df = function() vn_region_df(vn_v())),
    # Pathways & enrichment (cell-type level, KO vs WT pooled over P0+P7)
    list(id = "enr_gsea_tab", base = function() paste0("GSEA_", input$enr_ct, "_", input$enr_tp),
         df = function() enr_gsea_table_df(input$enr_ct, input$enr_tp)),
    list(id = "enr_go_tab", base = function() paste0("GO_", input$enr_ct, "_", input$enr_tp),
         df = function() enr_go_table_df(input$enr_ct, input$enr_tp)),
    # Maturation & metabolism gene map
    list(id = "gm_table", base = function() paste0("gene_map_", input$gm_panel %||% "avg"),
         df = function() gm_d()),
    # Dev / Help
    list(id = "cc_tab", base = function() paste0("cell_communication_", input$cc_tp),
         df = function() commun_table_df(input$cc_pathway, input$cc_tp)),
    list(id = "ann_tab", base = "annotation_check", df = function() refmap_table_df()),
    list(id = "doublet_tab", base = "qc_doublets", df = function() tabs$doublet),
    list(id = "score_meta_tab", base = "module_score_definitions", df = function() score_meta_df()),
    list(id = "score_sets_tab", base = "module_score_gene_sets", df = function() score_sets_df())
  )
  for (.t in TABLE_DL) local({
    t <- .t
    register_dl(output, t$id, t$df, t$base)
  })
}

shinyApp(ui, server)
