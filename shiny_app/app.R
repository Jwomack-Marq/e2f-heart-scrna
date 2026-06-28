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

app   <- readRDS("app_data.rds")
meta  <- app$meta; expr <- app$expr; genes <- app$genes
cmm   <- app$cm$meta; RES <- app$cm$res
heat  <- app$heat; tabs <- app$tables; CONF <- app$confound
ctDE  <- tabs$ct_DE; subDE <- tabs$sub_DE; subSum <- tabs$sub_summary; subType <- tabs$sub_subtype
figs  <- app$figs
GI    <- app$geneinfo          # per-gene info table (NULL on an un-enriched build)
ENR   <- app$enrich            # list(gsea, go, tf) of precomputed enrichment (may be NULL)
EXPR  <- app$deg_expr          # broad-gene log-norm matrix (genes x downsampled cells) for the DEG explorer
DMETA <- if (!is.null(app$deg_meta)) app$deg_meta else meta   # metadata aligned to EXPR cols (fallback so UI builds)
GENES_FULL <- app$deg_genes
# searchable union: curated panel (full 30k cells) + broad genes (deg_expr, ~8k cells).
# Lets the UMAP/Gene-detail views show ANY gene that shows up in a volcano, falling
# back to the broad matrix when a gene is outside the curated panel (expr_vec below).
ALL_GENES <- sort(unique(c(genes, GENES_FULL)))
in_panel  <- function(g) !is.null(g) && nzchar(g) && g %in% rownames(expr)

has <- function(col, df = meta) col %in% names(df)
CAT_COLS  <- Filter(has, c("celltype","genotype","timepoint","Phase","cycling","cm_subtype","seurat_clusters"))
CONT_COLS <- Filter(has, c("pseudotime","S.Score","G2M.Score"))

nice <- c(gene = "Gene expression", celltype = "Cell type", genotype = "Genotype (KO/WT)",
          timepoint = "Timepoint (P0/P7)", Phase = "Cell-cycle phase", cycling = "Cycling (S/G2M)",
          cm_subtype = "CM subtype", seurat_clusters = "Cluster", pseudotime = "Pseudotime",
          S.Score = "S-phase score", G2M.Score = "G2/M score", subcluster = "Subcluster")
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
                           xlab = "log2 fold change (KO / WT)") {
  validate(need(!is.null(d) && nrow(d), "No DE results for this selection (cluster too small / unbalanced)."))
  d <- de_annot(d, pos, neg)
  d$hover <- sprintf(
    "<b>%s</b><br>logFC: %.2f<br>-log10 p: %.2f<br>padj: %.2g<br>%s",
    d$gene, d$log2FoldChange, d$neglogp, d$padj, as.character(d$class))
  plot_ly(d, x = ~log2FoldChange, y = ~neglogp, color = ~class, colors = volc_pal(pos, neg),
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
      margin = list(t = 34)) |>
    event_register("plotly_click")
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
de_datatable <- function(df) {
  disp <- ifelse(names(df) %in% names(DE_DISP), DE_DISP[names(df)], names(df))
  DT::datatable(df, rownames = FALSE, selection = "single", colnames = unname(disp),
    options = list(pageLength = DE_PAGELEN, scrollY = "380px", scrollCollapse = TRUE,
                   dom = "ftip", order = list()),
    class = "compact stripe hover") |>
    DT::formatSignif(intersect(c("log2FoldChange","neglog10p","baseMean","pvalue","padj"), names(df)), 3)
}
# optionally drop the sex/construct confounder genes (Xist, Y-genes, ROSA26) from a
# DE frame before it reaches a volcano/table — n=1 makes these dominate the contrast.
drop_conf <- function(d, hide) {
  if (isTRUE(hide) && !is.null(d) && "confounder" %in% names(d)) d <- d[!d$confounder, , drop = FALSE]
  d
}
# small "Showing only <gene>" banner + a link to clear the volcano-click filter
pick_banner <- function(gene, clear_id) {
  if (is.null(gene) || !nzchar(gene)) return(NULL)
  div(style = "margin-bottom:6px; font-size:13px",
      span(HTML(paste0("Showing only <b>", gene, "</b> &middot; "))),
      actionLink(clear_id, "show all genes"))
}
# info card for a picked gene, from the bundled app$geneinfo table (GI). All data
# is precomputed (no runtime network calls); external links go to the full records.
gene_info_card <- function(gene) {
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
    card_header(HTML(paste0("<b>", gene, "</b> &mdash; ", name))),
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
enr_gsea <- function(ct, tp) { d <- ENR$gsea
  validate(need(!is.null(d), "GSEA results are not in this data build."))
  d[d$celltype == ct & d$timepoint == tp, , drop = FALSE] }
enr_gsea_plot <- function(ct, tp, topn = 20) {
  d <- enr_gsea(ct, tp)
  validate(need(nrow(d), "No GSEA results for this cell type / timepoint."))
  d <- head(d[order(-abs(d$NES)), ], topn)
  d$dir <- ifelse(d$NES > 0, "up in KO", "up in WT")
  d$pathway <- factor(d$pathway, levels = d$pathway[order(d$NES)])
  p <- ggplot(d, aes(NES, pathway, fill = dir,
        text = paste0(pathway, "<br>NES: ", NES, "<br>padj: ", padj, "<br>size: ", size))) +
    geom_col() + geom_vline(xintercept = 0, color = "grey60") +
    scale_fill_manual(values = c("up in KO" = "#c62828", "up in WT" = "#1565c0")) +
    theme_minimal(base_size = 12) +
    labs(x = "NES (>0 enriched toward KO-up)", y = NULL, fill = NULL,
         title = paste0("GSEA — ", ct, " ", tp))
  ggplotly(p, tooltip = "text") |> layout(margin = list(l = 0, t = 40))
}
enr_gsea_table <- function(ct, tp) {
  d <- enr_gsea(ct, tp)
  enr_dt(d[order(d$padj), intersect(c("pathway","NES","padj","size","leadingEdge"), names(d))])
}
enr_go <- function(ct, tp) { d <- ENR$go
  validate(need(!is.null(d), "GO results are not in this data build."))
  d[d$celltype == ct & d$timepoint == tp, , drop = FALSE] }
enr_go_plot <- function(ct, tp, topn = 20) {
  d <- enr_go(ct, tp)
  validate(need(nrow(d), "No GO BP results for this cell type / timepoint."))
  d <- head(d[order(d$p.adjust), ], topn)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  p <- ggplot(d, aes(FoldEnrichment, Description, size = Count, color = p.adjust,
        text = paste0(Description, "<br>fold: ", FoldEnrichment, "<br>padj: ", p.adjust, "<br>genes: ", Count))) +
    geom_point() + scale_color_viridis_c(option = "magma", direction = -1) +
    theme_minimal(base_size = 11) +
    labs(x = "fold enrichment", y = NULL, color = "padj", size = "genes",
         title = paste0("GO BP enriched in KO-up genes — ", ct, " ", tp))
  ggplotly(p, tooltip = "text") |> layout(margin = list(l = 0, t = 40))
}
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
.ax       <- list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "")
.exprsc   <- list(c(0,"#eeeeee"), c(0.45,"#fec44f"), c(0.75,"#fc4e2a"), c(1,"#800026"))
.umap_mar <- list(l = 0, r = 0, b = 0, t = 26)
centroids <- function(df) { a <- aggregate(cbind(UMAP1, UMAP2) ~ val, df, median)
  lapply(seq_len(nrow(a)), function(i) list(x = a$UMAP1[i], y = a$UMAP2[i], text = as.character(a$val[i]),
    showarrow = FALSE, font = list(size = 12, color = "#111"), bgcolor = "rgba(255,255,255,0.55)")) }

# single-panel categorical: hover a cell -> highlight ALL cells of that cluster (crosstalk),
# double-click a legend entry to isolate; centroid labels make clusters easy to read.
umap_cat <- function(df, colvar, ttl = NULL, psize = 4.5, labels = TRUE) {
  df$val <- factor(df[[colvar]])
  grp <- paste0("umaphl_", make.names(colvar))                 # named crosstalk group to clear from JS
  sd <- crosstalk::SharedData$new(df, key = ~val, group = grp)
  p <- plot_ly(sd, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
               color = ~val, colors = pal_for(levels(df$val)),
               marker = list(size = psize, line = list(width = 0)), text = ~val, hoverinfo = "text") |>
    layout(xaxis = .ax, yaxis = .ax, margin = .umap_mar, title = list(text = ttl, font = list(size = 13)),
           legend = list(itemsizing = "constant"),
           annotations = if (labels) centroids(df) else NULL)
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
umap_cont <- function(df, val, ttl = NULL, psize = 4.5) {
  df$val <- as.numeric(val); df <- df[order(df$val, na.last = FALSE), ]
  plot_ly(df, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
          marker = list(size = psize, color = ~val, colorscale = .exprsc, showscale = TRUE,
                        colorbar = list(title = ""), line = list(width = 0)),
          text = ~round(val, 2), hoverinfo = "text") |>
    layout(xaxis = .ax, yaxis = .ax, margin = .umap_mar, title = list(text = ttl, font = list(size = 13)))
}
# side-by-side split panels (legend / colourbar shown once)
umap_split <- function(df, colvar, splitvar, gene = NULL, continuous = FALSE, psize = 4) {
  levs_sp  <- sort(unique(as.character(df[[splitvar]])))
  val_levs <- if (!continuous) levels(factor(df[[colvar]])) else NULL
  cols     <- if (!continuous) pal_for(val_levs) else NULL
  plts <- lapply(seq_along(levs_sp), function(j) {
    d <- df[as.character(df[[splitvar]]) == levs_sp[j], ]
    if (continuous) {
      d$val <- if (!is.null(gene)) expr_vec(gene, d$cell) else d[[colvar]]; d <- d[order(d$val, na.last = FALSE), ]
      plot_ly(d, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
              marker = list(size = psize, color = ~val, colorscale = .exprsc, showscale = (j == 1),
                            colorbar = list(title = ""), line = list(width = 0)),
              text = ~round(val, 2), hoverinfo = "text") |> layout(xaxis = .ax, yaxis = .ax)
    } else {
      d$val <- factor(d[[colvar]], levels = val_levs); p <- plot_ly()
      for (i in seq_along(val_levs)) { di <- d[d$val == val_levs[i], ]
        p <- add_trace(p, data = di, x = ~UMAP1, y = ~UMAP2, type = "scattergl", mode = "markers",
                       name = val_levs[i], legendgroup = val_levs[i], showlegend = (j == 1),
                       marker = list(size = psize, color = cols[i], line = list(width = 0)),
                       text = val_levs[i], hoverinfo = "text") }
      p |> layout(xaxis = .ax, yaxis = .ax)
    }
  })
  anns <- lapply(seq_along(levs_sp), function(j) list(text = paste0(labof(splitvar), ": ", levs_sp[j]),
    x = (j - 0.5) / length(levs_sp), y = 1.0, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 13)))
  subplot(plts, nrows = 1, shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE) |>
    layout(margin = .umap_mar, legend = list(itemsizing = "constant"), annotations = anns)
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

# ---------------------------------------------------------------- UI ----------
ui <- page_navbar(
  title = "E2F7/8 heart scRNA-seq", theme = bs_theme(version = 5, bootswatch = "flatly"),

  nav_panel("UMAP explorer", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("color_by", "Colour cells by",
                  choices = setNames(c("gene", CAT_COLS, CONT_COLS),
                                     c(labof("gene"), labof(CAT_COLS), labof(CONT_COLS))),
                  selected = "celltype"),
      conditionalPanel("input.color_by == 'gene'",
        selectizeInput("gene", "Gene", choices = NULL, options = list(maxOptions = 50L))),
      selectInput("split", "Split panels by", c("(none)" = "none",
                  setNames(intersect(c("genotype","timepoint"), CAT_COLS),
                           labof(intersect(c("genotype","timepoint"), CAT_COLS))))),
      sliderInput("ptsize", "Point size", 1, 9, 4.5, 0.5),
      hr(), helpText("Hover a cell to highlight its whole cluster; double-click a legend entry to isolate one; ",
                     "single-click to toggle. Colour by a gene, then split by Genotype to compare KO vs WT.")),
    card(full_screen = TRUE, card_header(textOutput("umap_title")), plotlyOutput("umap", height = "640px")))),

  nav_panel("Gene detail", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectizeInput("g2", "Gene", choices = NULL, options = list(maxOptions = 50L)),
      selectInput("grp", "Group by", setNames(CAT_COLS, labof(CAT_COLS)),
                  selected = if (has("celltype")) "celltype" else CAT_COLS[1]),
      selectInput("sp2", "Split by", c("(none)" = "none",
                  setNames(intersect(c("genotype","timepoint"), CAT_COLS),
                           labof(intersect(c("genotype","timepoint"), CAT_COLS)))), selected = "genotype")),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Expression distribution (violin)"), plotOutput("vln", height = "460px")),
      card(card_header("% expressing & mean expression"), plotOutput("dot", height = "460px"))))),

  nav_panel("Composition", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("comp_fill", "Show fractions of", setNames(CAT_COLS, labof(CAT_COLS)),
                  selected = if (has("celltype")) "celltype" else CAT_COLS[1]),
      selectInput("comp_x", "Across groups", setNames(
                  intersect(c("orig.ident","genotype","timepoint"), names(meta)),
                  labof(intersect(c("orig.ident","genotype","timepoint"), names(meta)))),
                  selected = if (has("orig.ident")) "orig.ident" else "genotype")),
    card(full_screen = TRUE, card_header("Cell-type / state proportions"), plotOutput("comp", height = "560px")))),

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
        helpText("Hover a point for the gene & stats; click a point to show just that gene in the table."),
        layout_columns(col_widths = c(6, 6),
        plotlyOutput("ct_volcano", height = "470px"),
        div(uiOutput("ct_pick_ui"), DTOutput("ct_table", height = "440px"))),
        uiOutput("ct_geneinfo")),
      nav_panel("Heatmap (top genes × cell types)", plotlyOutput("ct_heat", height = "620px"))))),

  nav_panel("Cardiomyocyte deep-dive", layout_sidebar(
    sidebar = sidebar(width = 320,
      selectInput("cm_res", "Subcluster resolution",
                  setNames(RES, paste0("res ", RES)), selected = RES[length(RES)]),
      selectInput("cm_sub", "Subcluster (for DE)", choices = NULL),
      checkboxInput("cm_hideconf", "Hide sex/construct genes (DE)", FALSE),
      selectInput("cm_mapcolor", "Map: colour by",
                  c("Subcluster" = "subcluster", "Cell-cycle phase" = "Phase", "Cycling" = "cycling",
                    "Genotype" = "genotype", "Timepoint" = "timepoint", "Gene" = "gene")),
      conditionalPanel("input.cm_mapcolor == 'gene'",
        selectizeInput("cm_gene", "Gene", choices = NULL, options = list(maxOptions = 50L))),
      selectInput("cm_phase_split", "Cell-cycle plot: split by",
                  c("Genotype" = "genotype", "Timepoint" = "timepoint",
                    "Genotype × Timepoint" = "both"), selected = "both"),
      hr(), helpText("True re-clustering of cardiomyocytes. Explore subgroup identity,",
                     "KO-vs-WT differences per subgroup, and cell-cycle state.",
                     br(), "Split the cell-cycle plot by timepoint to see whether two cycling",
                     "subclusters separate by P0/P7 or by S vs G2/M phase.")),
    navset_card_tab(
      nav_panel("Subcluster map",
        helpText("Hover any cell to highlight all cells of its subcluster; move off to restore the full map."),
        plotlyOutput("cm_map", height = "600px")),
      nav_panel("Identity (marker heatmap)", plotlyOutput("cm_markerheat", height = "660px")),
      nav_panel("KO-vs-WT DE (per subgroup)",
        helpText("Hover a point for the gene & stats; click a point to show just that gene in the table."),
        layout_columns(col_widths = c(6, 6),
          plotlyOutput("cm_volcano", height = "440px"),
          div(uiOutput("cm_pick_ui"), DTOutput("cm_detab", height = "410px"))),
        uiOutput("cm_geneinfo"),
        card(card_header("DE heatmap — top genes × subclusters (log2FC KO/WT)"),
             plotlyOutput("cm_lfcheat", height = "560px"))),
      nav_panel("Cell cycle", plotOutput("cm_phase", height = "560px"))))),

  nav_panel("E2F focus", layout_sidebar(
    sidebar = sidebar(width = 320,
      selectInput("e2f_ct", "Cell type",
                  choices = c("All cells" = "All",
                              setNames(sort(unique(as.character(meta$celltype))),
                                       gsub("_", " ", sort(unique(as.character(meta$celltype)))))),
                  selected = if ("Cardiomyocyte" %in% meta$celltype) "Cardiomyocyte" else "All"),
      hr(),
      helpText(strong("Reading the E2F8 'up in KO' signal."),
               br(), "An apparent rise in E2f8 transcript in the E2F7/8 KO most often reflects ",
               "the conditional allele: reads outside the excised exons are still quantified, and ",
               "loss of E2F8 auto-repression can elevate the residual transcript. It does NOT ",
               "establish that functional E2F8 protein is up. Confirm which exons the flox removes ",
               "vs. where reads map before interpreting direction.")),
    card(full_screen = TRUE,
      card_header("E2f7 / E2f8 log-norm expression by genotype × timepoint"),
      plotOutput("e2f_expr", height = "560px")))),

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
    helpText("Hover a point for the gene & stats; click to show just that gene below."),
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
        plotlyOutput("enr_tf_top", height = "460px"))))),

  nav_panel("QC & normalization", div(style = "max-width:1000px;padding:8px 4px",
    uiOutput("qcfigs"),
    h5("Doublet rate by lane (numbers)"),
    div(style = "overflow:auto", tableOutput("doublet_tab")))),

  nav_panel("About / caveats", div(style = "max-width:820px;padding:8px 4px", htmlOutput("about")))
)

# -------------------------------------------------------------- SERVER --------
server <- function(input, output, session) {
  for (id in c("gene","g2","cm_gene")) {
    sel <- if ("Gabbr2" %in% ALL_GENES) "Gabbr2" else ALL_GENES[1]
    updateSelectizeInput(session, id, choices = ALL_GENES, selected = sel, server = TRUE)
  }
  # cell-type choices depend on timepoint
  observeEvent(input$ct_tp, {
    if (!length(input$ct_tp) || !nzchar(input$ct_tp)) return()
    cts <- sub(paste0("^", input$ct_tp, "_"), "", grep(paste0("^", input$ct_tp, "_"), names(ctDE), value = TRUE))
    updateSelectInput(session, "ct_sel", choices = setNames(cts, gsub("_", " ", cts)),
                      selected = if ("Cardiomyocyte" %in% cts) "Cardiomyocyte" else cts[1])
  }, ignoreNULL = FALSE)
  # subcluster choices depend on resolution
  observeEvent(input$cm_res, {
    if (!length(input$cm_res) || !nzchar(input$cm_res)) return()
    subs <- names(subDE[[paste0("res", input$cm_res)]])
    subs <- subs[order(as.integer(sub("CM", "", subs)))]
    updateSelectInput(session, "cm_sub", choices = setNames(subs, vapply(subs, function(s) sub_label(input$cm_res, s), "")),
                      selected = subs[1])
  }, ignoreNULL = FALSE)

  # ---- UMAP explorer ----
  output$umap_title <- renderText({
    if (input$color_by != "gene") return(paste0("Coloured by ", labof(input$color_by)))
    note <- if (!in_panel(input$gene)) "  (broad matrix — shown on the ~8k-cell subset)" else ""
    paste0("Expression of ", input$gene, note)
  })
  output$umap <- renderPlotly({
    cb <- input$color_by; cont <- cb %in% c("gene", CONT_COLS); ps <- input$ptsize
    if (input$split == "none") {
      if (cont) umap_cont(meta, if (cb == "gene") expr_vec(input$gene, meta$cell) else meta[[cb]],
                          if (cb == "gene") input$gene else labof(cb), psize = ps)
      else umap_cat(meta, cb, psize = ps)
    } else umap_split(meta, cb, input$split, gene = if (cb == "gene") input$gene else NULL,
                      continuous = cont, psize = max(2, ps - 1))
  })

  # ---- Gene detail ----
  output$vln <- renderPlot({
    df <- meta; df$expr <- expr_vec(input$g2, df$cell); df$grp <- factor(df[[input$grp]])
    base <- theme_minimal(base_size = 13) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    if (input$sp2 != "none") {
      df$splitv <- factor(df[[input$sp2]])
      ggplot(df, aes(grp, expr, fill = splitv)) +
        geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2, position = position_dodge(.9)) +
        base + labs(x = labof(input$grp), y = paste0(input$g2, " (log-norm)"), fill = labof(input$sp2))
    } else {
      ggplot(df, aes(grp, expr, fill = grp)) + geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2) +
        base + guides(fill = "none") + labs(x = labof(input$grp), y = paste0(input$g2, " (log-norm)"))
    }
  })
  output$dot <- renderPlot({
    df <- meta; df$expr <- expr_vec(input$g2, df$cell); df$grp <- factor(df[[input$grp]]); sp <- input$sp2
    key <- if (sp != "none") interaction(df$grp, df[[sp]], sep = " · ", drop = TRUE) else df$grp
    agg <- do.call(rbind, lapply(split(df, key), function(s) data.frame(
      grp = s$grp[1], split = if (sp != "none") as.character(s[[sp]][1]) else "all",
      pct = 100 * mean(s$expr > 0), mean = mean(s$expr))))
    ggplot(agg, aes(grp, split, size = pct, color = mean)) + geom_point() +
      scale_color_viridis_c(option = "magma", direction = -1) + scale_size_area(max_size = 12) +
      theme_minimal(base_size = 13) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
      labs(x = labof(input$grp), y = if (sp != "none") labof(sp) else "", size = "% expr", color = "mean", title = input$g2)
  })

  # ---- Composition ----
  output$comp <- renderPlot({
    df <- meta; x <- input$comp_x; f <- input$comp_fill
    tb <- as.data.frame(prop.table(table(df[[x]], df[[f]]), margin = 1)); names(tb) <- c("x","fill","prop")
    ggplot(tb, aes(x, prop, fill = fill)) + geom_col() + theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1)) + labs(x = labof(x), y = "proportion", fill = labof(f))
  })

  # ---- DE by cell type ----
  ct_d    <- reactive({ req(input$ct_tp, input$ct_sel); ctDE[[paste(input$ct_tp, input$ct_sel, sep = "_")]] })
  ct_pick <- reactiveVal(NULL)                                     # gene picked by clicking the volcano
  ct_tab  <- reactive({                                            # exactly what the table shows
    d <- de_table(drop_conf(ct_d(), input$ct_hideconf), input$ct_search); g <- ct_pick()
    if (!is.null(g) && g %in% d$gene) d <- d[d$gene == g, , drop = FALSE]
    d
  })
  output$ct_volcano <- renderPlotly(de_volcano_ly(drop_conf(ct_d(), input$ct_hideconf),
                         paste0(input$ct_tp, " ", gsub("_", " ", input$ct_sel)), "ct_volcano"))
  output$ct_table   <- renderDT(de_datatable(ct_tab()))
  output$ct_pick_ui  <- renderUI(pick_banner(ct_pick(), "ct_clear"))
  output$ct_geneinfo <- renderUI(gene_info_card(ct_pick()))
  observeEvent(event_data("plotly_click", source = "ct_volcano"),
    ct_pick(event_data("plotly_click", source = "ct_volcano")$customdata))
  observeEvent(input$ct_clear, ct_pick(NULL))
  observeEvent(list(input$ct_tp, input$ct_sel, input$ct_search), ct_pick(NULL), ignoreInit = TRUE)
  output$ct_heat <- renderPlotly({
    req(input$ct_tp)
    keys <- grep(paste0("^", input$ct_tp, "_"), names(ctDE), value = TRUE)
    dl <- setNames(ctDE[keys], gsub("_", " ", sub(paste0("^", input$ct_tp, "_"), "", keys)))
    ggheat(lfc_heat(dl, 24, paste0("Top KO-vs-WT genes across cell types — ", input$ct_tp)))
  })

  # ---- CM deep-dive ----
  output$cm_map <- renderPlotly({
    req(input$cm_res, input$cm_mapcolor); cb <- input$cm_mapcolor; df <- cmm
    if (cb == "gene") return(umap_cont(df, expr_vec(input$cm_gene, df$cell), input$cm_gene, psize = 5))
    df$mapval <- if (cb == "subcluster")
      factor(paste0("CM", df[[cm_subcol(input$cm_res)]]), levels = cm_subs(input$cm_res)) else factor(df[[cb]])
    umap_cat(df, "mapval", ttl = paste0(labof(cb), if (cb == "subcluster") paste0(" — res ", input$cm_res) else ""),
             psize = 5, labels = (cb == "subcluster"))
  })
  output$cm_markerheat <- renderPlotly({
    req(input$cm_res)
    h <- heat[[paste0("res", input$cm_res)]]; validate(need(!is.null(h), "No marker heatmap for this resolution."))
    long <- h$long; long$gene <- factor(long$gene, levels = rev(h$genes)); long$cluster <- factor(long$cluster, levels = h$clusters)
    p <- ggplot(long, aes(cluster, gene, fill = z)) + geom_tile() +
      scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
      theme_minimal(base_size = 11) + theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = "subcluster", y = "marker gene", fill = "z-score\nmean expr",
           title = paste0("Subcluster identity markers — res ", input$cm_res))
    ggheat(p)
  })
  cm_d    <- reactive({ req(input$cm_res, input$cm_sub); subDE[[paste0("res", input$cm_res)]][[input$cm_sub]] })
  cm_pick <- reactiveVal(NULL)
  cm_tab  <- reactive({
    d <- de_table(drop_conf(cm_d(), input$cm_hideconf)); g <- cm_pick()
    if (!is.null(g) && g %in% d$gene) d <- d[d$gene == g, , drop = FALSE]
    d
  })
  output$cm_volcano <- renderPlotly(de_volcano_ly(drop_conf(cm_d(), input$cm_hideconf),
                         paste0(input$cm_sub, " — ", sub_label(input$cm_res, input$cm_sub)), "cm_volcano"))
  output$cm_detab   <- renderDT(de_datatable(cm_tab()))
  output$cm_pick_ui  <- renderUI(pick_banner(cm_pick(), "cm_clear"))
  output$cm_geneinfo <- renderUI(gene_info_card(cm_pick()))
  observeEvent(event_data("plotly_click", source = "cm_volcano"),
    cm_pick(event_data("plotly_click", source = "cm_volcano")$customdata))
  observeEvent(input$cm_clear, cm_pick(NULL))
  observeEvent(list(input$cm_res, input$cm_sub), cm_pick(NULL), ignoreInit = TRUE)
  output$cm_lfcheat <- renderPlotly({ req(input$cm_res)
    ggheat(lfc_heat(subDE[[paste0("res", input$cm_res)]], 22,
             paste0("KO-vs-WT log2FC across CM subclusters — res ", input$cm_res))) })
  output$cm_phase <- renderPlot({
    req(input$cm_res)
    df <- cmm; df$sub <- factor(paste0("CM", df[[cm_subcol(input$cm_res)]]), levels = cm_subs(input$cm_res))
    validate(need("Phase" %in% names(df), "No cell-cycle phase data."))
    sp <- input$cm_phase_split %||% "both"
    has_tp <- "timepoint" %in% names(df)
    split_cols <- switch(sp,
      genotype = "genotype",
      timepoint = if (has_tp) "timepoint" else "genotype",
      both = if (has_tp) c("genotype","timepoint") else "genotype")
    # fraction of each Phase within sub × (split groups): margin = every dim but Phase (dim 2)
    tab <- table(df[c("sub","Phase", split_cols)])
    tb <- as.data.frame(prop.table(tab, setdiff(seq_along(dim(tab)), 2L)))
    names(tb)[match("Freq", names(tb))] <- "prop"
    tb$prop[is.nan(tb$prop)] <- 0     # empty subcluster × group cells -> 0, not NaN
    fw <- if (length(split_cols) == 2) facet_grid(reformulate(split_cols[2], split_cols[1])) else facet_wrap(reformulate(split_cols))
    ggplot(tb, aes(sub, prop, fill = Phase)) + geom_col() + fw +
      scale_fill_manual(values = c(G1 = "#bdbdbd", S = "#1565c0", G2M = "#c62828"), na.value = "grey90") +
      theme_minimal(base_size = 13) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      labs(x = "subcluster", y = "fraction of cells",
           title = paste0("Cell-cycle phase by subcluster — res ", input$cm_res))
  })

  # ---- E2F focus (E2f7 / E2f8 expression by genotype x timepoint) ----
  output$e2f_expr <- renderPlot({
    req(input$e2f_ct)
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
    long$genotype <- factor(long$genotype)
    ggplot(long, aes(genotype, expr, fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, alpha = .55, linewidth = .2) +
      stat_summary(fun = mean, geom = "point", size = 2.4, color = "black") +
      facet_grid(gene ~ timepoint, scales = "free_y") +
      scale_fill_manual(values = c(KO = "#c62828", WT = "#1565c0"), na.value = "grey70") +
      theme_minimal(base_size = 13) + guides(fill = "none") +
      labs(x = "genotype", y = "log-norm expression",
           title = paste0("E2f7 / E2f8 — ", gsub("_", " ", input$e2f_ct),
                          " (black dot = mean; descriptive, n = 1)"))
  })

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
  deg_tab <- reactive({
    d <- de_table(drop_conf(deg_res(), input$deg_hideconf)); g <- deg_pick()
    if (!is.null(g) && g %in% d$gene) d <- d[d$gene == g, , drop = FALSE]
    d
  })
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
                  xlab = "logFC  (A / B)")
  })
  output$deg_table   <- renderDT(de_datatable(deg_tab()))
  output$deg_pick_ui <- renderUI(pick_banner(deg_pick(), "deg_clear"))
  output$deg_geneinfo <- renderUI(gene_info_card(deg_pick()))
  observeEvent(event_data("plotly_click", source = "deg_volcano"),
    deg_pick(event_data("plotly_click", source = "deg_volcano")$customdata))
  observeEvent(input$deg_clear, deg_pick(NULL))

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
}

shinyApp(ui, server)
