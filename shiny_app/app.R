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
de_annot <- function(d, pos = "up in KO", neg = "up in WT") {
  d$neglogp <- -log10(pmax(d$pvalue, 1e-300))
  d$class <- ifelse(d$confounder, "sex/construct",
              ifelse(abs(d$log2FoldChange) >= 1, ifelse(d$log2FoldChange > 0, pos, neg), "n.s."))
  d$class <- factor(d$class, levels = c(pos, neg, "n.s.", "sex/construct"))
  d
}
# DE volcano from a trimmed DE data frame (static ggplot — kept for reference)
de_volcano <- function(d, ttl) {
  validate(need(!is.null(d) && nrow(d), "No DE results for this selection (cluster too small / unbalanced)."))
  d <- de_annot(d)
  ggplot(d, aes(log2FoldChange, neglogp, color = class)) +
    geom_point(size = 1.1, alpha = .6) + scale_color_manual(values = VOLC_PAL) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey60") +
    theme_minimal(base_size = 13) +
    labs(x = "log2 fold change (KO / WT)", y = "-log10 p (ranking only, n=1)", color = NULL, title = ttl)
}
# interactive plotly volcano: hover shows gene/stats, click emits the gene via
# customdata (captured by event_data(source = source_id)) to drive the DE table.
de_volcano_ly <- function(d, ttl, source_id, pos = "up in KO", neg = "up in WT",
                           xlab = "log2 fold change (KO / WT)", highlight = NULL) {
  validate(need(!is.null(d) && nrow(d), "No DE results for this selection (cluster too small / unbalanced)."))
  d <- de_annot(d, pos, neg)
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
      shapes = lapply(c(-1, 1), function(v) list(type = "line", x0 = v, x1 = v,
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
gsea_barplot_gg <- function(d, ttl, topn = 20) {
  validate(need(!is.null(d) && nrow(d), "No GSEA results for this selection."))
  d <- head(d[order(-abs(d$NES)), ], topn)
  d$dir <- ifelse(d$NES > 0, "up in KO", "up in WT")
  d$pathway <- factor(d$pathway, levels = d$pathway[order(d$NES)])
  ggplot(d, aes(NES, pathway, fill = dir,
        text = paste0(pathway, "<br>NES: ", NES, "<br>padj: ", padj, "<br>size: ", size))) +
    geom_col() + geom_vline(xintercept = 0, color = "grey60") +
    scale_fill_manual(values = c("up in KO" = "#c62828", "up in WT" = "#1565c0")) +
    theme_minimal(base_size = 12) +
    labs(x = "NES (>0 enriched toward KO-up)", y = NULL, fill = NULL, title = ttl)
}
gsea_barplot_df <- function(d, ttl, topn = 20)
  ggplotly(gsea_barplot_gg(d, ttl, topn), tooltip = "text") |> layout(margin = list(l = 0, t = 40))
go_dotplot_gg <- function(d, ttl, topn = 20) {
  validate(need(!is.null(d) && nrow(d), "No GO BP results for this selection."))
  d <- head(d[order(d$p.adjust), ], topn)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  ggplot(d, aes(FoldEnrichment, Description, size = Count, color = p.adjust,
        text = paste0(Description, "<br>fold: ", FoldEnrichment, "<br>padj: ", p.adjust, "<br>genes: ", Count))) +
    geom_point() + scale_color_viridis_c(option = "magma", direction = -1) +
    theme_minimal(base_size = 11) +
    labs(x = "fold enrichment", y = NULL, color = "padj", size = "genes", title = ttl)
}
go_dotplot_df <- function(d, ttl, topn = 20)
  ggplotly(go_dotplot_gg(d, ttl, topn), tooltip = "text") |> layout(margin = list(l = 0, t = 40))
# static "all subclusters at once" overviews (faceted; renderPlot, low-memory) --
# "All clusters" enrichment view: one full-size plot per subcluster, two per row
# (server registers a renderPlot per subcluster id "cm_<kind>_all_<CMn>").
cm_enr_grid <- function(kind) {
  subs <- cm_subs("0.2")
  outs <- lapply(subs, function(cl) plotOutput(paste0("cm_", kind, "_all_", cl), height = "340px"))
  rows <- lapply(seq(1, length(outs), by = 2), function(i)
    fluidRow(column(6, outs[[i]]),
             if (i + 1L <= length(outs)) column(6, outs[[i + 1L]])))
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
enr_gsea_plot <- function(ct, tp, topn = 20)
  gsea_barplot_df(enr_gsea(ct, tp), paste0("GSEA — ", ct, " ", tp), topn)
enr_gsea_table <- function(ct, tp) {
  d <- enr_gsea(ct, tp)
  enr_dt(d[order(d$padj), intersect(c("pathway","NES","padj","size","leadingEdge"), names(d))])
}
enr_go <- function(ct, tp) { d <- ENR$go
  validate(need(!is.null(d), "GO results are not in this data build."))
  d[d$celltype == ct & d$timepoint == tp, , drop = FALSE] }
enr_go_plot <- function(ct, tp, topn = 20)
  go_dotplot_df(enr_go(ct, tp), paste0("GO BP enriched in KO-up genes — ", ct, " ", tp), topn)
enr_go_table <- function(ct, tp) {
  d <- enr_go(ct, tp)
  enr_dt(d[order(d$p.adjust), intersect(c("ID","Description","FoldEnrichment","p.adjust","Count","geneID"), names(d))])
}
# E2F-family regulon activity (KO - WT) across cell type x timepoint
enr_e2f_heat <- function() {
  d <- tabs$e2f_regulon
  validate(need(!is.null(d) && nrow(d), "No E2F regulon activity table."))
  d$col <- paste0(d$celltype, " ", d$timepoint)
  p <- ggplot(d, aes(col, source, fill = KO_minus_WT,
        text = paste0(source, "<br>", col, "<br>KO-WT: ", round(KO_minus_WT, 3)))) +
    geom_tile(color = "grey92") +
    scale_fill_gradient2(low = "#1565c0", mid = "white", high = "#c62828", midpoint = 0) +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "KO - WT", title = "E2F-family regulon activity (KO - WT)")
  ggplotly(p, tooltip = "text") |> layout(margin = list(t = 40))
}
# top TFs by |KO - WT| activity for a cell type (from ENR$tf)
enr_tf_top <- function(ct, topn = 20) {
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
    scale_fill_manual(values = c("TRUE" = "#c62828", "FALSE" = "#1565c0"), guide = "none") +
    theme_minimal(base_size = 11) +
    labs(x = "KO - WT activity", y = NULL, title = paste0("Top TFs by |KO-WT| — ", ct))
  ggplotly(p, tooltip = "text") |> layout(margin = list(l = 0, t = 40))
}
# log2FC heatmap: top genes (by max |LFC| across groups) x groups, fill = KO/WT log2FC
lfc_heat <- function(de_list, topn = 22, ttl = NULL) {
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
    labs(x = NULL, y = NULL, fill = "log2FC\n(KO/WT)", title = ttl)
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

# ---- publication "Figure options": reusable control block + export wiring -------
# emits a namespaced (prefix_*) control set; `export` picks the download UI (vector
# ggsave for ggplot panels, camera-button PNG for the WebGL UMAP).
figure_controls <- function(prefix, export = c("ggplot","umap"), base = TRUE, palette = TRUE,
                            rename = TRUE, labels = FALSE, default_base = 13) {
  export <- match.arg(export); p <- function(s) paste0(prefix, "_", s)
  tagList(
    textInput(p("title"), "Title (blank = default)", ""),
    if (labels) checkboxInput(p("labels"), "Show cluster labels", TRUE),
    if (labels) sliderInput(p("labelsize"), "Label font size", 8, 28, 12, 1),
    if (base)   sliderInput(p("basesize"), "Base font size", 8, 24, default_base, 1),
    checkboxInput(p("legend"), "Show legend", TRUE),
    if (palette) selectInput(p("palette"), "Palette", names(PALETTES_DISCRETE)),
    if (rename) tagList(tags$small("Rename categories (double-click a label cell):"),
                        DTOutput(p("renametab"))),
    hr(),
    checkboxInput(p("export_on"), "Custom export options", FALSE),
    conditionalPanel(sprintf("input.%s_export_on", prefix),
      if (export == "ggplot") tagList(
        div(style = "display:flex;gap:6px",
          numericInput(p("w"), "W (in)", 7, 1, 20, 0.5),
          numericInput(p("h"), "H (in)", 5, 1, 20, 0.5),
          numericInput(p("dpi"), "DPI", 300, 72, 600, 1)),
        div(style = "display:flex;gap:6px;margin-top:4px",
          downloadButton(p("dl_pdf"), "PDF", class = "btn-sm btn-outline-secondary"),
          downloadButton(p("dl_svg"), "SVG", class = "btn-sm btn-outline-secondary"),
          downloadButton(p("dl_png"), "PNG", class = "btn-sm btn-outline-secondary")))
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
dl_ggplot <- function(prefix, plot_reactive, input, fmt) {
  downloadHandler(
    filename = function() paste0(prefix, "_", Sys.Date(), ".", fmt),
    content  = function(file) {
      g <- function(s) input[[paste0(prefix, "_", s)]]
      ggsave(file, apply_fig_opts(plot_reactive(), prefix, input), device = fmt,
             width = g("w") %||% 7, height = g("h") %||% 5, dpi = g("dpi") %||% 300, units = "in")
    })
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
commun_heat <- function(pathway, tp, metric) {
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
  ggplotly(p, tooltip = "text") |> layout(margin = list(t = 40))
}
commun_table <- function(pathway, tp) {
  validate(need(!is.null(COMMUN$scores), "No communication table in this build."))
  d <- COMMUN$scores[COMMUN$scores$timepoint == tp, , drop = FALSE]
  if (pathway != "All") d <- d[d$pathway == pathway, , drop = FALSE]
  d <- d[order(-abs(d$delta)), ]
  enr_dt(d[, intersect(c("pathway","ligand","receptor","sender","receiver","WT","KO","delta"), names(d))])
}
# ---- reference-marker annotation check (REFMAP$confusion) --------------------
refmap_heat <- function() {
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
  ggplotly(p, tooltip = "text") |> layout(margin = list(t = 40))
}
refmap_table <- function() {
  validate(need(!is.null(REFMAP$confusion), "No annotation-check table in this build."))
  d <- REFMAP$confusion
  enr_dt(d[order(d$celltype, -d$prop), ])
}

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
fg_skip_msg <- function(cluster, contrast, stratum) {
  s <- FG$skipped
  if (!is.null(s)) {
    r <- s[s$cluster == cluster & s$contrast == contrast & s$stratum == stratum, , drop = FALSE]
    if (nrow(r)) return(sprintf(
      "Not computed — %s. Group A: %d cells, group B: %d cells (floor is %d).%s",
      r$reason[1], r$n_A[1], r$n_B[1], FG$built$min_cells,
      if (stratum == "G1") " Try the “All cells” stratum." else ""))
  }
  "No DE table for this selection."
}
fg_de <- function(cluster, contrast, stratum) {
  fg_ok(); req(cluster, contrast, stratum)
  d <- FG$de[[cluster]][[paste0(contrast, "__", stratum)]]
  validate(need(!is.null(d) && nrow(d), fg_skip_msg(cluster, contrast, stratum)))
  d
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
fg_specificity_df <- function(genes, stratum) {
  fg_ok(); validate(need(length(genes), "Pick at least one gene."))
  pull <- function(cl, key, g) {
    d <- FG$de[[cl]][[paste0(key, "__", stratum)]]
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
  d[order(d$gene, factor(d$cluster, levels = FG_CLUSTERS)), , drop = FALSE]
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

# `panel` selects which pair of coordinates to plot: the timepoint-averaged axes, or
# P0 / P7 on their own. Quadrant and distance are recomputed against that panel's own
# centre rather than reused from the averaged one.
gm_df <- function(panel = "avg", quadrants = NULL, min_dist = 0,
                  hide_sets = TRUE, geneset = "__all__") {
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
  # genes inside the scoring sets sit at the extremes of their own axis by
  # construction; hidden by default so the map isn't just recovering its own inputs
  if (isTRUE(hide_sets)) d <- d[is.na(d$in_score_set), , drop = FALSE]
  if (!is.null(geneset) && geneset != "__all__")
    d <- d[d$gene %in% genes_for_set(geneset), , drop = FALSE]
  if (!is.null(quadrants) && length(quadrants)) d <- d[d$quadrant %in% quadrants, , drop = FALSE]
  if (!is.null(min_dist)) d <- d[d$distance >= min_dist, , drop = FALSE]
  validate(need(nrow(d), "No genes pass these filters — lower the distance cut or re-enable a quadrant."))
  d[order(-d$distance), , drop = FALSE]
}
# table view: the columns worth reading, plus the per-timepoint AUCs so a gene's
# movement between P0 and P7 is visible without switching panels
gm_table <- function(d) {
  tps <- setdiff(GM_PANELS, "avg")
  cols <- intersect(c("gene","quadrant","distance","x","y",
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
    card(full_screen = TRUE, card_header(textOutput("umap_title")), plotlyOutput("umap", height = "640px")))),

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

  nav_menu("Differential expression",
  nav_panel("DE by cell type", layout_sidebar(
    sidebar = sidebar(width = 300,
      radioButtons("ct_tp", "Timepoint", c("P0","P7"), inline = TRUE),
      selectInput("ct_sel", "Cell type", choices = NULL),
      textInput("ct_search", "Filter genes (substring)", ""),
      checkboxInput("ct_hideconf", "Hide sex/construct genes", FALSE),
      hr(), helpText("KO-vs-WT differential expression within each cell type.",
                     br(), strong("p-axis ranks candidates only — not valid at n = 1."))),
    navset_card_tab(
      nav_panel("Volcano + table",
        helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it and show its info below."),
        layout_columns(col_widths = c(6, 6),
          plotlyOutput("ct_volcano", height = "470px"),
          div(uiOutput("ct_pick_ui"), DTOutput("ct_table", height = "440px"))),
        uiOutput("ct_geneinfo")),
      nav_panel("Heatmap (top genes × cell types)", plotlyOutput("ct_heat", height = "620px"))))),

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
      actionButton("deg_run", "Compute DEGs", class = "btn-primary"),
      hr(), helpText("Descriptive Wilcoxon (presto) on log-norm expression of the ",
                     "filtered live cells. Hypothesis-generating only (n = 1); ",
                     "for rigorous KO-vs-WT use the precomputed DE tabs.")),
    div(textOutput("deg_n"), style = "font-size:13px;margin-bottom:4px"),
    helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it and show its info below."),
    layout_columns(col_widths = c(6, 6),
      plotlyOutput("deg_volcano", height = "470px"),
      div(uiOutput("deg_pick_ui"), DTOutput("deg_table", height = "440px"))),
    uiOutput("deg_geneinfo"))),

  nav_panel("Pathways & enrichment", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("enr_tp", "Timepoint", c("P0","P7"), selected = "P7"),
      selectInput("enr_ct", "Cell type", choices = NULL),
      hr(), helpText("Pre-computed pathway/GO/TF enrichment of the KO-vs-WT signal ",
                     "(fgsea Hallmark/KEGG/E2F, GO biological process, decoupleR TF activity). ",
                     strong("Descriptive only — n = 1."))),
    navset_card_tab(
      nav_panel("GSEA pathways",
        plotlyOutput("enr_gsea_plot", height = "440px"),
        DTOutput("enr_gsea_tab", height = "360px")),
      nav_panel("GO biological process",
        plotlyOutput("enr_go_plot", height = "440px"),
        DTOutput("enr_go_tab", height = "360px")),
      nav_panel("TF / regulon activity",
        helpText("E2F-family regulon activity across cell types (KO − WT), then the top TFs for the selected cell type."),
        plotlyOutput("enr_e2f_heat", height = "380px"),
        plotlyOutput("enr_tf_top", height = "460px")))))),

  nav_menu("Cardiomyocytes",
  nav_panel("Cardiomyocyte deep-dive", layout_sidebar(
    sidebar = sidebar(width = 320,
      conditionalPanel("input.cm_tabs == 'de'",
        selectInput("cm_sub", "Subcluster (for DE)", choices = NULL),
        checkboxInput("cm_hideconf", "Hide sex/construct genes (DE)", FALSE)),
      selectInput("cm_mapcolor", "Map: colour by",
                  c("Subcluster" = "subcluster", "Cell-cycle phase" = "Phase", "Cycling" = "cycling",
                    "Genotype" = "genotype", "Timepoint" = "timepoint", "Gene" = "gene")),
      conditionalPanel("input.cm_mapcolor == 'gene'",
        selectInput("cm_geneset", "Gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("cm_gene", "Gene", choices = NULL, options = list(maxOptions = 50L))),
      selectInput("cm_map_split", "Split map by (res 0.2)",
                  c("(none)" = "none", "Genotype (WT|KO)" = "genotype",
                    "Timepoint (P0|P7)" = "timepoint", "Genotype × Timepoint" = "both")),
      conditionalPanel("input.cm_tabs == 'bars'",
        radioButtons("cm_bar_mode", "Y axis", c("Proportion" = "prop", "Count" = "count"), inline = TRUE)),
      hr(), helpText("True re-clustering of cardiomyocytes. Explore subgroup identity,",
                     "KO-vs-WT differences per subgroup, and cell-cycle state.",
                     br(), "Split the map by timepoint to see whether two cycling subclusters",
                     "separate by P0/P7 or by S vs G2/M phase."),
      accordion(open = FALSE, multiple = TRUE,
        accordion_panel("Cell-cycle figure options",
          figure_controls("cmphase", palette = TRUE, rename = TRUE)),
        accordion_panel("Composition figure options",
          figure_controls("cmbar", palette = TRUE, rename = FALSE)))),
    navset_card_tab(id = "cm_tabs",
      nav_panel("Subcluster map",
        helpText("Hover any cell to highlight all cells of its subcluster; move off to restore the full map.",
                 br(), "When split (res 0.2) is on, panels are coloured by subcluster and “colour by” is ignored."),
        plotlyOutput("cm_map", height = "600px")),
      nav_panel("Identity (marker heatmap)", plotlyOutput("cm_markerheat", height = "660px")),
      nav_panel("KO-vs-WT DE (per subgroup)", value = "de",
        helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it and show its info below."),
        fluidRow(
          column(6, plotlyOutput("cm_volcano", height = "440px")),
          column(6, uiOutput("cm_pick_ui"), DTOutput("cm_detab"))),
        uiOutput("cm_geneinfo"),
        div(class = "mt-3",
            h5("DE heatmap — top genes × subclusters (log2FC KO/WT)"),
            plotlyOutput("cm_lfcheat", height = "560px"))),
      nav_panel("Cell cycle", plotOutput("cm_phase", height = "560px")),
      nav_panel("Composition (stacked bars)", value = "bars",
        helpText("Composition of each res-0.2 subcluster, broken down four ways — genotype (WT/KO), ",
                 "timepoint (P0/P7), cell-cycle phase, and cycling status. Each is its own plot; scroll to see all four."),
        plotOutput("cm_bar_geno",  height = "300px"),
        plotOutput("cm_bar_tp",    height = "300px"),
        plotOutput("cm_bar_phase", height = "300px"),
        plotOutput("cm_bar_cyc",   height = "300px")),
      nav_panel("Per-cluster summary", value = "summary",
        helpText("One row per res-0.2 subcluster: identity, composition, top marker & KO-vs-WT genes, and top pathways — scan all clusters without clicking through tabs. Sizes/percentages are from the displayed (sampled) cells."),
        div(downloadButton("cm_summary_dl", "Download CSV", class = "btn-sm btn-outline-secondary"),
            style = "margin-bottom:8px"),
        DTOutput("cm_summary")),
      nav_panel("Top markers", value = "topmarkers",
        helpText("One row per res-0.2 subcluster: top identity markers (by z-scored mean expression) plus ",
                 "each cluster's top identity GO term, top KO-vs-WT GSEA pathway, and top KO-up / KO-down genes — ",
                 "a quick read on what each subcluster is doing biologically."),
        DTOutput("cm_topmarkers")),
      nav_panel("Subcluster enrichment", value = "subenr",
        helpText("Per res-0.2 subcluster: identity markers and the KO-vs-WT signal, enriched. ",
                 strong("Descriptive only — n = 1.")),
        div(class = "d-flex align-items-center gap-3 mb-2",
            radioButtons("cm_enr_mode", NULL, c("Single cluster" = "one", "All clusters" = "all"),
                         selected = "one", inline = TRUE),
            conditionalPanel("input.cm_enr_mode == 'one'",
              selectInput("cm_enr_sub", NULL, choices = NULL, width = "260px"))),
        conditionalPanel("input.cm_enr_mode == 'one'",
          navset_card_tab(
            nav_panel("Identity GO",
              plotlyOutput("cm_sub_idgo_plot", height = "440px"),
              DTOutput("cm_sub_idgo_tab", height = "320px")),
            nav_panel("KO-vs-WT GO",
              plotlyOutput("cm_sub_kogo_plot", height = "440px"),
              DTOutput("cm_sub_kogo_tab", height = "320px")),
            nav_panel("KO-vs-WT GSEA",
              plotlyOutput("cm_sub_gsea_plot", height = "440px"),
              DTOutput("cm_sub_gsea_tab", height = "320px")))),
        conditionalPanel("input.cm_enr_mode == 'all'",
          navset_card_tab(
            nav_panel("Identity GO", cm_enr_grid("idgo")),
            nav_panel("KO-vs-WT GO", cm_enr_grid("kogo")),
            nav_panel("KO-vs-WT GSEA", cm_enr_grid("gsea")))))))),

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
        radioButtons("fg_stratum", "Cells used",
                     c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                     selected = "G1"),
        checkboxInput("fg_hideconf", "Hide sex/construct genes", FALSE),
        div(downloadButton("fg_de_dl", "Download DEG table (CSV)",
                           class = "btn-sm btn-outline-secondary"), style = "margin-bottom:8px")),
      conditionalPanel("input.fg_tabs == 'counts'",
        radioButtons("fg_count_mode", "Y axis",
                     c("% of subcluster" = "prop", "Cell count" = "count"), inline = TRUE),
        div(downloadButton("fg_counts_dl", "Download counts (CSV)",
                           class = "btn-sm btn-outline-secondary"), style = "margin-bottom:8px")),
      conditionalPanel("input.fg_tabs == 'g1'",
        selectizeInput("fg_g1_clusters", "Subclusters", choices = NULL, multiple = TRUE),
        selectInput("fg_score", "Maturation / state score", choices = NULL),
        radioButtons("fg_score_stratum", "Cells used",
                     c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                     selected = "all"),
        div(downloadButton("fg_scores_dl", "Download score summary (CSV)",
                           class = "btn-sm btn-outline-secondary"), style = "margin-bottom:8px")),
      hr(),
      helpText(strong("Sort caveat. "), FG_SORT_NOTE),
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
        DTOutput("fg_counts_tab")),
      nav_panel("Four-group DE", value = "de",
        helpText("Hover a point for the gene & stats; click a point — or a table row — to highlight it."),
        uiOutput("fg_de_note"),
        fluidRow(
          column(6, plotlyOutput("fg_volcano", height = "440px")),
          column(6, uiOutput("fg_pick_ui"), DTOutput("fg_detab"))),
        uiOutput("fg_geneinfo")),
      nav_panel("G1 & maturation", value = "g1",
        helpText("Top: cell-cycle phase composition per group (G1 % printed on each bar). ",
                 "Bottom: per-cell maturation / state scores. ",
                 "The question is whether P7 KO cardiomyocytes sit at a less mature score than P7 WT ",
                 "— compare within the G1 stratum to hold cycling composition fixed."),
        plotOutput("fg_phase_plot", height = "420px"),
        plotOutput("fg_score_plot", height = "440px"))))),

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
        checkboxInput("mi_hideconf", "Hide sex/construct genes", TRUE),
        div(downloadButton("mi_dl", "Download intersection (CSV)",
                           class = "btn-sm btn-outline-secondary"), style = "margin-bottom:8px")),
      conditionalPanel("input.mi_tabs == 'candidates'",
        selectInput("mi_geneset", "Gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("mi_genes", "Genes", choices = NULL, multiple = TRUE,
                       options = list(maxOptions = 50L)),
        actionLink("mi_reset_genes", "reset to the shortlist"),
        selectizeInput("mi_cand_clusters", "Subclusters", choices = NULL, multiple = TRUE),
        radioButtons("mi_cand_stratum", "Cells used",
                     c("G1 only (phase-matched)" = "G1", "All cells (raw)" = "all"),
                     selected = "all"),
        div(class = "mt-2",
          downloadButton("mi_cand_dl", "Download expression grid (CSV)",
                         class = "btn-sm btn-outline-secondary"),
          downloadButton("mi_spec_dl", "Download P7-specificity (CSV)",
                         class = "btn-sm btn-outline-secondary"))),
      hr(),
      helpText(strong("Maturation axis. "),
               "Genes are ranked by comparing the most- vs least-mature cardiomyocytes, ",
               "within each timepoint and then averaged, so the axis is maturation and not P0-vs-P7. ",
               "It uses the ", strong("cycle-free"), " maturation score (Mki67 / Top2a / Ccnd1 removed) ",
               "— otherwise “less mature ⇒ more cycling” would be partly circular."),
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
        DTOutput("mi_table")),
      nav_panel("Candidate genes", value = "candidates",
        helpText("Any gene, across CM subclusters × the four groups. ",
                 "Size = % of cells expressing, colour = mean expression. ",
                 "The table below reads the KO effect off the precomputed contrasts: ",
                 strong("P7_specificity"), " > 0 means the KO effect is larger at P7 than at P0, and ",
                 strong("priority_concentration"), " > 0 means it is larger inside CM2/CM4/CM5 than outside."),
        plotOutput("mi_candidates", height = "480px"),
        DTOutput("mi_spec_tab")))))),

  nav_menu("Dev",
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
        plotlyOutput("cc_heat", height = "540px")),
      nav_panel("Interaction table", DTOutput("cc_tab", height = "520px"))))),

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
        plotOutput("cyc_violins", height = "600px")),
      nav_panel("Cycling vs cytokinesis", plotOutput("cyc_scatter", height = "600px"))))),

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
        checkboxInput("gm_hidesets", "Hide genes from the scoring sets", TRUE),
        selectInput("gm_geneset", "Restrict to gene set", choices = GENE_SET_CHOICES, selected = "__all__"),
        selectizeInput("gm_gene", "Find a gene", choices = NULL, options = list(maxOptions = 50L)),
        numericInput("gm_labeln", "Label top N by distance", 20, 0, 200, 5),
        div(downloadButton("gm_dl", "Download gene map (CSV)",
                           class = "btn-sm btn-outline-secondary"), style = "margin-bottom:8px")),
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
                 "their own axis by construction, so leaving them in would partly just recover the inputs.")),
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
      nav_panel("By four groups", value = "violin", plotOutput("mat_violin", height = "560px")),
      nav_panel("Maturation vs metabolism (cells)", value = "cells",
        plotOutput("mat_scatter", height = "600px")),
      nav_panel("Gene map", value = "genemap",
        helpText("Hover a point for the gene and its coordinates; click a point — or a table row — ",
                 "to ring it and show its details."),
        uiOutput("gm_note"),
        plotlyOutput("gm_scatter", height = "560px"),
        uiOutput("gm_pick_ui"),
        DTOutput("gm_table"),
        uiOutput("gm_geneinfo")))))),

  nav_spacer(),
  nav_menu("Help",
  nav_panel("QC & normalization", div(style = "max-width:1000px;padding:8px 4px",
    uiOutput("qcfigs"),
    h5("Doublet rate by lane (numbers)"),
    div(style = "overflow:auto", tableOutput("doublet_tab")))),

  nav_panel("Annotation check", div(style = "max-width:1000px;padding:8px 4px",
    helpText("Each cell is scored against published-style developmental mouse-heart marker panels; ",
             "the argmax lineage (\"predicted\") is cross-tabulated against the existing cell-type label. ",
             "The diagonal should dominate — off-diagonal mass flags populations worth double-checking. ",
             strong("Concordance check, not probabilistic label transfer.")),
    plotlyOutput("ann_heat", height = "520px"),
    h5("Confusion (row-normalised)", style = "margin-top:12px"),
    DTOutput("ann_tab", height = "360px"))),

  nav_panel("About / caveats", div(style = "max-width:820px;padding:8px 4px", htmlOutput("about"))))
)

# -------------------------------------------------------------- SERVER --------
server <- function(input, output, session) {
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
    updateSelectInput(session, "cm_sub", choices = setNames(subs, vapply(subs, function(s) sub_label("0.2", s), "")),
                      selected = subs[1])
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
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("vln_dl_", f)]] <- dl_ggplot("vln", vln_plot, input, f) })

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
  for (.f in c("pdf","svg","png")) local({ f <- .f; output[[paste0("dot_dl_", f)]] <- dl_ggplot("dot", dot_plot, input, f) })
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
                         highlight = ct_pick()))
  outputOptions(output, "ct_volcano", suspendWhenHidden = FALSE)  # render at startup so plotly_click source registers before its click observer fires
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
  output$ct_heat <- renderPlotly({
    req(input$ct_tp)
    keys <- grep(paste0("^", input$ct_tp, "_"), names(ctDE), value = TRUE)
    dl <- setNames(ctDE[keys], gsub("_", " ", sub(paste0("^", input$ct_tp, "_"), "", keys)))
    ggheat(lfc_heat(dl, 24, paste0("Top KO-vs-WT genes across cell types — ", input$ct_tp)))
  })

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
  output$cm_markerheat <- renderPlotly({
    h <- heat[["res0.2"]]; validate(need(!is.null(h), "No marker heatmap for this resolution."))
    long <- h$long; long$gene <- factor(long$gene, levels = rev(h$genes)); long$cluster <- factor(long$cluster, levels = h$clusters)
    p <- ggplot(long, aes(cluster, gene, fill = z)) + geom_tile() +
      scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
      theme_minimal(base_size = 11) + theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = "subcluster", y = "marker gene", fill = "z-score\nmean expr",
           title = "Subcluster identity markers — res 0.2")
    ggheat(p)
  })
  cm_d    <- reactive({ req(input$cm_sub); subDE[["res0.2"]][[input$cm_sub]] })
  cm_pick <- reactiveVal(NULL)
  cm_tab  <- reactive(de_table(drop_conf(cm_d(), input$cm_hideconf)))   # full table (selection never filters it)
  cm_dt_proxy <- DT::dataTableProxy("cm_detab")
  output$cm_volcano <- renderPlotly(de_volcano_ly(drop_conf(cm_d(), input$cm_hideconf),
                         paste0(input$cm_sub, " — ", sub_label("0.2", input$cm_sub)), "cm_volcano",
                         highlight = cm_pick()))
  outputOptions(output, "cm_volcano", suspendWhenHidden = FALSE)  # render at startup so plotly_click source registers before its click observer fires
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
  observeEvent(input$cm_sub, cm_pick(NULL), ignoreInit = TRUE)
  output$cm_lfcheat <- renderPlotly(
    ggheat(lfc_heat(subDE[["res0.2"]], 22, "KO-vs-WT log2FC across CM subclusters — res 0.2")))
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
  output$cm_bar_geno  <- renderPlot(apply_fig_opts(cm_bar_one("genotype",  "By genotype (WT / KO)"), "cmbar", input))
  output$cm_bar_tp    <- renderPlot(apply_fig_opts(cm_bar_one("timepoint", "By timepoint (P0 / P7)"), "cmbar", input))
  output$cm_bar_phase <- renderPlot(apply_fig_opts(cm_bar_one("Phase",     "By cell-cycle phase"), "cmbar", input))
  output$cm_bar_cyc   <- renderPlot(apply_fig_opts(cm_bar_one("cycling",   "By cycling status"), "cmbar", input))
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
  output$cm_summary_dl <- downloadHandler(
    filename = function() paste0("cm_subcluster_summary_res0.2_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(cm_summary_df(), f, row.names = FALSE))

  # ---- CM top markers + biology per subcluster (stylized; res 0.2) ----
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

  # ---- CM per-subcluster enrichment (precomputed in ENR$sub; res 0.2) ----
  local({
    avail <- if (!is.null(ENR$sub)) sort(unique(c(ENR$sub$identity_go$subcluster,
                                                   ENR$sub$go$subcluster, ENR$sub$gsea$subcluster))) else character(0)
    subs <- cm_subs("0.2"); subs <- if (length(avail)) subs[subs %in% avail] else subs
    updateSelectInput(session, "cm_enr_sub",
                      choices = setNames(subs, vapply(subs, function(s) sub_label("0.2", s), "")),
                      selected = subs[1])
  })
  output$cm_sub_idgo_plot <- renderPlotly({ req(input$cm_enr_sub)
    go_dotplot_df(enr_sub_df("identity_go", input$cm_enr_sub), paste0("Identity GO BP — ", input$cm_enr_sub)) })
  output$cm_sub_idgo_tab <- renderDT({ req(input$cm_enr_sub)
    d <- enr_sub_df("identity_go", input$cm_enr_sub)
    enr_dt(d[order(d$p.adjust), intersect(c("ID","Description","FoldEnrichment","p.adjust","Count","geneID"), names(d))]) })
  output$cm_sub_kogo_plot <- renderPlotly({ req(input$cm_enr_sub)
    d <- enr_sub_df("go", input$cm_enr_sub); d <- d[d$direction == "KO_up", , drop = FALSE]
    go_dotplot_df(d, paste0("GO BP enriched in KO-up genes — ", input$cm_enr_sub)) })
  output$cm_sub_kogo_tab <- renderDT({ req(input$cm_enr_sub)
    d <- enr_sub_df("go", input$cm_enr_sub)
    enr_dt(d[order(d$p.adjust), intersect(c("Description","direction","FoldEnrichment","p.adjust","Count","geneID"), names(d))]) })
  output$cm_sub_gsea_plot <- renderPlotly({ req(input$cm_enr_sub)
    gsea_barplot_df(enr_sub_df("gsea", input$cm_enr_sub), paste0("GSEA — ", input$cm_enr_sub)) })
  output$cm_sub_gsea_tab <- renderDT({ req(input$cm_enr_sub)
    d <- enr_sub_df("gsea", input$cm_enr_sub)
    enr_dt(d[order(d$padj), intersect(c("pathway","NES","padj","size","leadingEdge"), names(d))]) })
  # "All clusters" view: one full-size plot per subcluster (two per row on screen)
  local({
    for (.cl in cm_subs("0.2")) local({
      cl <- .cl; ttl <- sub_label("0.2", cl)
      output[[paste0("cm_idgo_all_", cl)]] <- renderPlot(
        go_dotplot_gg(enr_sub_df("identity_go", cl), ttl, topn = 8))
      output[[paste0("cm_kogo_all_", cl)]] <- renderPlot({
        d <- enr_sub_df("go", cl); go_dotplot_gg(d[d$direction == "KO_up", , drop = FALSE], ttl, topn = 8) })
      output[[paste0("cm_gsea_all_", cl)]] <- renderPlot(
        gsea_barplot_gg(enr_sub_df("gsea", cl), ttl, topn = 10))
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
  observe(updateSelectizeInput(session, "mi_genes", choices = genes_for_set(input$mi_geneset),
                               selected = FG_SHORTLIST, server = TRUE))
  observeEvent(input$mi_reset_genes,
    updateSelectizeInput(session, "mi_genes", selected = FG_SHORTLIST, server = TRUE))

  # group sizes
  fg_counts_p <- reactive(fg_counts_plot(input$fg_count_mode %||% "prop",
                                         input$fgcount_basesize %||% 13))
  output$fg_counts_plot <- renderPlot(apply_fig_opts(fg_counts_p(), "fgcount", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("fgcount_dl_", f)]] <- dl_ggplot("fgcount", fg_counts_p, input, f) })
  output$fg_counts_tab <- renderDT(DT::datatable(fg_counts_wide(), rownames = FALSE,
    options = list(pageLength = 14, scrollX = TRUE, dom = "ft"),
    class = "compact stripe hover"))
  output$fg_counts_dl <- downloadHandler(
    filename = function() paste0("fourgroup_counts_res", FG$built$res, "_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(FG$counts, f, row.names = FALSE))

  # four-group DE: volcano <-> table <-> gene card (same wiring as the other DE tabs)
  fg_d    <- reactive(drop_conf(fg_de(input$fg_cluster, input$fg_contrast, input$fg_stratum),
                                input$fg_hideconf))
  fg_tab  <- reactive(de_table(fg_d()))
  fg_pick <- reactiveVal(NULL)
  fg_dt_proxy <- DT::dataTableProxy("fg_detab")
  output$fg_de_note <- renderUI({
    ct <- fg_ct(input$fg_contrast); req(ct)
    n <- try(fg_d(), silent = TRUE)
    bad <- !inherits(n, "try-error")
    nA <- if (bad) n$n_A[1] else NA; nB <- if (bad) n$n_B[1] else NA
    txt <- if (!bad) "" else sprintf(" — %d vs %d cells, %d genes past the gate", nA, nB, nrow(n))
    thin <- FG$built$thin_cells %||% 50
    div(style = "font-size:13px;margin-bottom:6px",
        HTML(paste0("<b>", ct$label, "</b> in ",
                    if (input$fg_cluster == "AllCM") "all cardiomyocytes" else input$fg_cluster,
                    if (input$fg_stratum == "G1") ", G1 cells only" else ", all cells", txt)),
        # An arm can clear the 10-cell floor and still be far too thin to trust —
        # say so here rather than letting the table look like any other.
        if (bad && min(nA, nB) < thin)
          div(style = "color:#c62828;font-weight:600",
              HTML(sprintf("&#9888; Smallest group is %d cells — treat this contrast as unreliable%s.",
                min(nA, nB),
                if (max(nA, nB) / max(min(nA, nB), 1) >= 10)
                  sprintf(" (%.0f× imbalance)", max(nA, nB) / max(min(nA, nB), 1)) else ""))),
        if (input$fg_stratum == "all" && grepl("P0_vs_P7", input$fg_contrast))
          span(style = "color:#c62828", HTML(" &middot; <b>sort-confounded</b> — see the caveat in the sidebar.")))
  })
  output$fg_volcano <- renderPlotly({
    ct <- fg_ct(input$fg_contrast); req(ct)
    de_volcano_ly(fg_d(),
      paste0(ct$label, " — ", if (input$fg_cluster == "AllCM") "all CM" else input$fg_cluster,
             if (input$fg_stratum == "G1") " (G1)" else ""),
      "fg_volcano", pos = ct$pos, neg = ct$neg, xlab = ct$xlab, highlight = fg_pick())
  })
  outputOptions(output, "fg_volcano", suspendWhenHidden = FALSE)  # register the click source at startup
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
  observeEvent(list(input$fg_cluster, input$fg_contrast, input$fg_stratum),
               fg_pick(NULL), ignoreInit = TRUE)
  output$fg_de_dl <- downloadHandler(
    filename = function() paste0("fourgroup_DE_", input$fg_cluster, "_", input$fg_contrast,
                                 "_", input$fg_stratum, "_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(fg_d(), f, row.names = FALSE))

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
  output$fg_scores_dl <- downloadHandler(
    filename = function() paste0("fourgroup_scores_", Sys.Date(), ".csv"),
    content  = function(f) {
      validate(need(!is.null(FG$scores), "No score summary in this data build."))
      write.csv(FG$scores, f, row.names = FALSE) })

  # ---- Maturation intersection + candidate genes ----
  mi_quad_p <- reactive(fg_quadrant_plot(input$mi_cluster, input$mi_hideconf,
                                         bs = input$miquad_basesize %||% 13))
  output$mi_quadrant <- renderPlot(apply_fig_opts(mi_quad_p(), "miquad", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("miquad_dl_", f)]] <- dl_ggplot("miquad", mi_quad_p, input, f) })
  mi_tab_df <- reactive(fg_intersect_df(input$mi_cluster, input$mi_quad, input$mi_hideconf))
  output$mi_table <- renderDT(DT::datatable(mi_tab_df(), rownames = FALSE,
    options = list(pageLength = 25, scrollX = TRUE, scrollY = "420px",
                   scrollCollapse = TRUE, dom = "ftip"),
    class = "compact stripe hover") |>
    DT::formatSignif(intersect(c("mat_log2FC","mat_auc","p7ko_log2FC","p7ko_padj"), names(mi_tab_df())), 3))
  output$mi_dl <- downloadHandler(
    filename = function() paste0("maturation_intersect_", input$mi_cluster, "_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(mi_tab_df(), f, row.names = FALSE))

  mi_cand_p <- reactive(fg_candidate_plot(input$mi_genes, input$mi_cand_clusters,
                                          input$mi_cand_stratum %||% "all",
                                          input$micand_basesize %||% 13))
  output$mi_candidates <- renderPlot(apply_fig_opts(mi_cand_p(), "micand", input))
  for (.f in c("pdf","svg","png")) local({ f <- .f;
    output[[paste0("micand_dl_", f)]] <- dl_ggplot("micand", mi_cand_p, input, f) })
  mi_spec_df <- reactive(fg_specificity_df(input$mi_genes, input$mi_cand_stratum %||% "all"))
  output$mi_spec_tab <- renderDT(DT::datatable(mi_spec_df(), rownames = FALSE,
    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"),
    class = "compact stripe hover") |>
    DT::formatSignif(intersect(c("P7_KO_vs_WT","P0_KO_vs_WT","P7_specificity",
                                 "CM2_4_5_mean_absLFC","other_clusters_mean_absLFC",
                                 "priority_concentration"), names(mi_spec_df())), 3))
  output$mi_cand_dl <- downloadHandler(
    filename = function() paste0("candidate_expression_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(
      fg_candidate_df(input$mi_genes, input$mi_cand_clusters, input$mi_cand_stratum %||% "all"),
      f, row.names = FALSE))
  output$mi_spec_dl <- downloadHandler(
    filename = function() paste0("candidate_P7_specificity_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(mi_spec_df(), f, row.names = FALSE))

  # ---- Subset & DEGs (interactive descriptive DE) ----
  observeEvent(input$deg_by, {
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
                  xlab = "logFC  (A / B)", highlight = deg_pick())
  })
  outputOptions(output, "deg_volcano", suspendWhenHidden = FALSE)  # render at startup so plotly_click source registers before its click observer fires
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
  output$enr_gsea_plot <- renderPlotly({ req(input$enr_ct); enr_gsea_plot(input$enr_ct, input$enr_tp) })
  output$enr_gsea_tab  <- renderDT({ req(input$enr_ct); enr_gsea_table(input$enr_ct, input$enr_tp) })
  output$enr_go_plot   <- renderPlotly({ req(input$enr_ct); enr_go_plot(input$enr_ct, input$enr_tp) })
  output$enr_go_tab    <- renderDT({ req(input$enr_ct); enr_go_table(input$enr_ct, input$enr_tp) })
  output$enr_e2f_heat  <- renderPlotly(enr_e2f_heat())
  output$enr_tf_top    <- renderPlotly({ req(input$enr_ct); enr_tf_top(input$enr_ct) })

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
  output$cyc_scatter <- renderPlot(ploidy_scatter(input$cyc_ct, input$cyc_basesize %||% 13))

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
                            input$gm_hidesets %||% TRUE, input$gm_geneset %||% "__all__"))
  gm_tab  <- reactive(gm_table(gm_d()))
  gm_pick <- reactiveVal(NULL)
  gm_dt_proxy <- DT::dataTableProxy("gm_table")
  output$gm_scatter <- renderPlotly(gm_plot_ly(gm_d(), input$gm_labeln %||% 20, gm_pick(),
                                                input$gm_panel %||% "avg"))
  outputOptions(output, "gm_scatter", suspendWhenHidden = FALSE)  # register the click source at startup
  output$gm_note <- renderUI({
    d <- try(gm_d(), silent = TRUE)
    n <- if (inherits(d, "try-error")) 0 else nrow(d)
    tagList(
      div(style = "font-size:13px;margin-bottom:4px",
          sprintf("%d genes shown of %d on the map — %s.", n, if (is.null(GM)) 0L else nrow(GM),
                  if (identical(input$gm_panel %||% "avg", "avg")) "P0 and P7 averaged"
                  else paste0(input$gm_panel, " cells only")),
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
  output$gm_dl <- downloadHandler(
    filename = function() paste0("gene_map_", input$gm_panel %||% "avg", "_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(gm_d(), f, row.names = FALSE))

  # ---- Cell-cell signalling (curated L-R; build_communication.R) ----
  output$cc_heat <- renderPlotly(commun_heat(input$cc_pathway, input$cc_tp, input$cc_metric))
  output$cc_tab  <- renderDT(commun_table(input$cc_pathway, input$cc_tp))

  # ---- Annotation check (reference-marker concordance; build_refmap.R) ----
  output$ann_heat <- renderPlotly(refmap_heat())
  output$ann_tab  <- renderDT(refmap_table())
}

shinyApp(ui, server)
