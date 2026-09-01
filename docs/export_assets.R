#!/usr/bin/env Rscript
# export_assets.R
# ---------------------------------------------------------------------------
# Produce every figure and generated table the methods book (docs/) embeds.
#
#   docs/export.sh                     # all of it, in e2f-enrich:latest
#   docs/export.sh --only=fourgroup    # just the ids matching a regex
#   docs/export.sh --list              # print the export ids and exit
#
# Two rules this script exists to enforce:
#
#   1. A FIGURE IN THE BOOK IS THE FIGURE IN THE APP. It is not redrawn from a
#      description. The app's plotting functions are loaded out of app.R itself --
#      everything above the UI marker -- and called directly, so `go_dotplot_gg()`
#      here and `go_dotplot_gg()` in the running app are the same closure over the
#      same bundle. Where a panel is a server reactive rather than a callable
#      function (there is no way to call `comp_plot` from outside the server) the
#      logic is mirrored in the MIRRORS section below, each one carrying the app.R
#      line range it mirrors so the two can be diffed by eye.
#
#   2. NO NUMBER IN THE BOOK IS TYPED BY HAND. Load-bearing tables are written to
#      docs/_generated/*.md as markdown fragments and included by the chapters, so
#      a rebuilt bundle changes the book. Same discipline build_deck.R follows, for
#      the same reason: prose drifts from data, and a methods document is where
#      that drift becomes a wrong claim in front of an audience.
#
# READ-ONLY on the bundle. Unlike shiny_app/build_*.R this never calls saveRDS --
# the book documents the analysis, it does not redo it.
# ---------------------------------------------------------------------------

suppressMessages({
  library(Matrix); library(ggplot2); library(shiny); library(bslib)
  library(plotly);  library(DT);      library(svglite); library(ragg)
})

ARGS   <- commandArgs(trailingOnly = TRUE)
ONLY   <- sub("^--only=", "", grep("^--only=", ARGS, value = TRUE))
LIST   <- "--list" %in% ARGS
BUNDLE <- if (file.exists("/in/app_data.rds")) "/in/app_data.rds" else "shiny_app/app_data.rds"
REPO   <- if (dir.exists("/repo")) "/repo" else "."
OUT_A  <- if (dir.exists("/out/assets")) "/out/assets" else "docs/assets"
OUT_G  <- if (dir.exists("/out/generated")) "/out/generated" else "docs/_generated"
DELIV  <- file.path(REPO, "deliverables/2026-08-21")
# PowerPoint cannot be trusted with SVG through pandoc's pptx writer, and no
# SVG->PNG converter exists on this machine. So when a slide-asset dir is mounted,
# every vector figure also gets a raster twin from the SAME plot object -- not a
# converted file, the same ggplot rendered through a second device. One plotting
# source, two output formats, no conversion step to drift.
SLIDES <- Sys.getenv("SLIDE_PNG_DIR", "")
if (nzchar(SLIDES) && !dir.exists(SLIDES)) dir.create(SLIDES, recursive = TRUE)

stopifnot(file.exists(BUNDLE))
dir.create(OUT_A, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_G, showWarnings = FALSE, recursive = TRUE)

# ---- load the app's own functions ------------------------------------------
# app.R sources download_helpers.R and readRDS()s app_data.rds relative to its own
# directory, so give it a staging cwd where both resolve. Symlinks, not copies: the
# repo and the bundle are mounted read-only and must stay that way.
STAGE <- file.path(tempdir(), "appstage")
dir.create(STAGE, showWarnings = FALSE)
# Link ALL of shiny_app/*.R, never a hand-listed subset: app.R source()s its helpers
# at startup and reads build_signature_scores.R for the curated gene lists, so the set
# of runtime files grows whenever app.R grows a dependency. A hand-maintained list here
# breaks silently the first time that happens -- which is exactly how the dev container
# failed to boot in 1c1b1f7, and how this script broke when studio_helpers.R appeared.
for (f in list.files(file.path(REPO, "shiny_app"), pattern = "\\.R$", full.names = TRUE)) {
  tgt <- file.path(STAGE, basename(f))
  if (!file.exists(tgt)) file.symlink(f, tgt)
}
if (!file.exists(file.path(STAGE, "app_data.rds")))
  file.symlink(BUNDLE, file.path(STAGE, "app_data.rds"))
APP_R <- file.path(REPO, "shiny_app/app.R")
src   <- readLines(APP_R, warn = FALSE)
ui_at <- grep("^# -{5,} UI -{5,}$", src)
stopifnot("app.R's UI marker is not unique -- the prefix boundary moved" = length(ui_at) == 1L)
cat(sprintf("Loading app.R functions (lines 1-%d of %d) ...\n", ui_at - 1L, length(src)))
owd <- setwd(STAGE)
eval(parse(text = paste(src[seq_len(ui_at - 1L)], collapse = "\n")), envir = globalenv())
setwd(owd)
cat(sprintf("  bundle built %s | %s cells | %s CM | %s panel genes | %s broad genes\n",
            as.character(app$built), format(nrow(meta), big.mark = ","),
            format(nrow(cmm), big.mark = ","), format(length(genes), big.mark = ","),
            format(length(GENES_FULL), big.mark = ",")))

# ---- bookkeeping ------------------------------------------------------------
MANIFEST <- list()
FAILED   <- character(0)
SKIPPED  <- character(0)
GAPS     <- list()

# A slot that was never built is not an export failure -- it is a documented gap,
# and the book has to say so out loud. Anything else lets a chapter imply a result
# that does not exist in this bundle. Recorded here, rendered into _gaps.md, and
# included by the chapter that would otherwise have shown it.
gap <- function(id, chapter, what, why) {
  GAPS[[length(GAPS) + 1L]] <<- data.frame(
    chapter = chapter, missing = what, why = why, stringsAsFactors = FALSE)
  SKIPPED <<- c(SKIPPED, id)
  cat(sprintf("  %-34s GAP - %s\n", id, why))
}

register <- function(id, chapter, file, what) {
  MANIFEST[[length(MANIFEST) + 1L]] <<- data.frame(
    id = id, chapter = chapter, file = file, what = what, stringsAsFactors = FALSE)
}
wanted <- function(id) !length(ONLY) || grepl(ONLY, id)

# `expr` is lazily evaluated inside the handler, so a failure in one export is
# reported and the rest still run -- a missing bundle slot should not cost the
# other 30 figures.
step <- function(id, chapter, what, expr) {
  if (!wanted(id)) { SKIPPED <<- c(SKIPPED, id); return(invisible()) }
  if (LIST) { cat(sprintf("  %-34s %s\n", id, what)); return(invisible()) }
  cat(sprintf("  %-34s ", id))
  tryCatch({
    f <- force(expr)
    register(id, chapter, f, what)
    cat("ok\n")
  }, error = function(e) {
    FAILED <<- c(FAILED, sprintf("%s: %s", id, conditionMessage(e)))
    cat("FAIL -", conditionMessage(e), "\n")
  })
}

# ---- writers ---------------------------------------------------------------
# SVG, not PNG: these are vector plots, the files are an order of magnitude smaller,
# they stay sharp when a boss zooms in on a projector, and the repo already tracks
# SVG figures (model/figures/*.svg).
# `raster = TRUE` for point clouds. A 30,000-cell UMAP as SVG is 4.4 MB of
# individual <circle> elements -- it bloats the repo, and a browser stutters
# scrolling past it. Rasterising those at 200 dpi costs nothing a reader can see
# and buys a factor of ~20. Everything else (bars, tiles, dot plots, lollipops)
# stays vector, where the text stays crisp under a projector's zoom.
fig <- function(id, plot, w = 8, h = 5, raster = FALSE, dpi = 200) {
  ext  <- if (raster) ".png" else ".svg"
  path <- file.path(OUT_A, paste0(id, ext))
  tmp  <- paste0(path, ".part")
  # write to a scratch name and move on success, so a plot that errors mid-render
  # cannot leave a 0-byte file behind that later looks like a successful export.
  ok <- FALSE
  on.exit(if (!ok) unlink(tmp), add = TRUE)
  if (raster) ggsave(tmp, plot, device = ragg::agg_png, width = w, height = h,
                     dpi = dpi, units = "in", bg = "white", limitsize = FALSE)
  else        ggsave(tmp, plot, device = svglite::svglite, width = w, height = h,
                     bg = "white", limitsize = FALSE)
  ok <- TRUE
  file.rename(tmp, path)
  # raster twin for the slides, at print dpi
  if (nzchar(SLIDES)) {
    png2 <- file.path(SLIDES, paste0(id, ".png"))
    if (raster) file.copy(path, png2, overwrite = TRUE)
    else ggsave(png2, plot, device = ragg::agg_png, width = w, height = h,
                dpi = 200, units = "in", bg = "white", limitsize = FALSE)
  }
  paste0("assets/", basename(path))
}
# The five upstream QC figures live in the bundle as base64 data URIs (app$figs) --
# they were rendered by the upstream pipeline, which is not in this repo, so they
# are the only figures here that cannot be regenerated from code.
figs_png <- function(id, slot) {
  d <- app$figs[[slot]]
  if (is.null(d) || !nzchar(d)) stop("app$figs$", slot, " is empty")
  b64  <- sub("^data:image/[a-z]+;base64,", "", d)
  path <- file.path(OUT_A, paste0(id, ".png"))
  writeBin(base64enc_decode(b64), path)
  if (nzchar(SLIDES)) file.copy(path, file.path(SLIDES, basename(path)), overwrite = TRUE)
  paste0("assets/", basename(path))
}
# base64 decoder in base R -- avoids depending on base64enc/openssl being present.
base64enc_decode <- function(x) {
  x <- gsub("[^A-Za-z0-9+/=]", "", x)
  chars <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1]]
  v <- match(strsplit(sub("=+$", "", x), "")[[1]], chars) - 1L
  bits <- unlist(lapply(v, function(i) as.integer(intToBits(i))[6:1]))
  n <- length(bits) %/% 8L
  bits <- bits[seq_len(n * 8L)]
  as.raw(colSums(matrix(bits, nrow = 8L) * 2^(7:0)))
}

# markdown pipe table from a data.frame -- Quarto renders it as a table and it stays
# diffable in git, which a rendered image of a table would not be.
md_table <- function(df, digits = 3, max_rows = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!is.null(max_rows) && nrow(df) > max_rows) df <- head(df, max_rows)
  fmt <- function(x) {
    if (is.numeric(x)) {
      blank <- function(v) { v[is.na(x)] <- ""; v }   # NA reads as absent, not as "NA"
      if (all(is.na(x) | x == round(x)) && max(abs(x), na.rm = TRUE) >= 1000)
        return(blank(formatC(x, format = "d", big.mark = ",")))
      if (all(is.na(x) | x == round(x))) return(blank(as.character(x)))
      # Hypergeometric and Wilcoxon p-values here run to 1e-280. "fg" would print
      # that as 280 leading zeros and destroy the table, so any column holding a
      # value outside a readable range switches to scientific notation wholesale
      # -- per column, not per cell, so a column stays visually comparable.
      nz <- x[is.finite(x) & x != 0]
      sci <- length(nz) && (min(abs(nz)) < 1e-4 || max(abs(nz)) >= 1e6)
      if (sci) return(ifelse(is.na(x), "", formatC(x, format = "g", digits = digits)))
      return(formatC(signif(x, digits), format = "fg", digits = digits, flag = "#"))
    }
    x <- as.character(x); x[is.na(x)] <- ""
    # a literal newline inside a cell ends the row and breaks the whole table;
    # the app's multi-line labels (set names carry their contrast on a second
    # line) hit this. Collapse to a separator that reads the same inline.
    x <- gsub("[\r\n]+", " · ", x)
    gsub("|", "\\|", x, fixed = TRUE)
  }
  cols <- lapply(df, fmt)
  hdr  <- gsub("_", " ", names(df))
  algn <- vapply(df, function(x) if (is.numeric(x)) "---:" else ":---", "")
  rows <- vapply(seq_len(nrow(df)), function(i)
    paste0("| ", paste(vapply(cols, function(cl) cl[i], ""), collapse = " | "), " |"), "")
  c(paste0("| ", paste(hdr, collapse = " | "), " |"),
    paste0("| ", paste(algn, collapse = " | "), " |"), rows)
}
frag <- function(id, df, caption = NULL, digits = 3, max_rows = NULL, note = NULL) {
  # An empty frame renders as a two-line table with no columns -- valid markdown,
  # silently wrong, and exactly the kind of thing that survives review. Fail here.
  if (is.null(df) || !NROW(df) || !NCOL(df))
    stop("no rows -- the bundle slot behind this table is empty or absent")
  path <- file.path(OUT_G, paste0(id, ".md"))
  body <- md_table(df, digits, max_rows)
  if (!is.null(caption)) body <- c(body, "", paste0(": ", caption, " {tbl-colwidths=\"auto\"}"))
  if (!is.null(note))    body <- c(body, "", paste0("::: {.callout-note appearance=\"minimal\"}\n",
                                                    note, "\n:::"))
  writeLines(c(paste0("<!-- generated by docs/export_assets.R -- do not edit -->"), "", body), path)
  paste0("_generated/", basename(path))
}

# ---- MIRRORS ---------------------------------------------------------------
# Panels the app builds inside server reactives (not callable from here). Each is a
# line-for-line mirror of the cited block; if a mirror and the app ever disagree,
# the app is right and this is a bug.
theme_pub <- theme_minimal(base_size = 13)

# mirrors umap_cat()/umap_cont() (app.R:435-472) onto a static device. The app's
# UMAP is WebGL plotly; there is no vector export of it, so this is the one figure
# whose *rendering* differs from the app while the data and palette do not.
umap_gg <- function(df, colvar, ttl = NULL, psize = .22, labels = TRUE, continuous = FALSE) {
  d <- df
  if (continuous) {
    d$val <- as.numeric(df[[colvar]]); d <- d[order(d$val, na.last = FALSE), ]
    ggplot(d, aes(UMAP1, UMAP2, colour = val)) + geom_point(size = psize, alpha = .7) +
      scale_colour_gradientn(colours = c("#eeeeee","#fec44f","#fc4e2a","#800026"),
                             values = c(0, .45, .75, 1), na.value = "grey92") +
      theme_umap + labs(title = ttl, colour = NULL)
  } else {
    d$val <- factor(d[[colvar]])
    cen <- aggregate(cbind(UMAP1, UMAP2) ~ val, d, median)
    p <- ggplot(d, aes(UMAP1, UMAP2, colour = val)) + geom_point(size = psize, alpha = .55) +
      scale_colour_manual(values = disc_pal(levels(d$val), "Default")) +
      guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
      theme_umap + labs(title = ttl, colour = labof(colvar))
    if (labels) p <- p + geom_text(data = cen, aes(UMAP1, UMAP2, label = val),
                                   inherit.aes = FALSE, size = 3, colour = "grey10")
    p
  }
}
# mirrors comp_plot (app.R:2905-2914)
comp_gg <- function(x, f) {
  d <- meta
  # the app offers genotype x timepoint as a composition axis; it is derived at
  # display time rather than stored as a column, so derive it the same way here
  d$splitgrp <- factor(paste(d$genotype, d$timepoint, sep = "-"), levels = FG_GROUPS)
  tb <- as.data.frame(prop.table(table(d[[x]], d[[f]]), margin = 1))
  names(tb) <- c("x", "fill", "prop")
  ggplot(tb, aes(x, prop, fill = fill)) + geom_col() + theme_pub +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(x = labof(x), y = "proportion", fill = labof(f))
}
# mirrors vln_plot (app.R:2856-2872) and dot_plot (app.R:2877-2890)
vln_gg <- function(gene, grp, split = NULL) {
  df <- meta; df$expr <- expr_vec(gene, df$cell); df$grp <- factor(df[[grp]])
  base <- theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  if (!is.null(split)) {
    df$splitv <- factor(df[[split]])
    ggplot(df, aes(grp, expr, fill = splitv)) +
      geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2,
                  position = position_dodge(.9)) +
      base + labs(x = labof(grp), y = paste0(gene, " (log-norm)"), fill = labof(split))
  } else {
    ggplot(df, aes(grp, expr, fill = grp)) +
      geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2) +
      base + guides(fill = "none") + labs(x = labof(grp), y = paste0(gene, " (log-norm)"))
  }
}
dot_gg <- function(gene, grp, split) {
  df <- meta; df$expr <- expr_vec(gene, df$cell); df$grp <- factor(df[[grp]])
  key <- interaction(df$grp, df[[split]], sep = " · ", drop = TRUE)
  agg <- do.call(rbind, lapply(split(df, key), function(s) data.frame(
    grp = s$grp[1], split = as.character(s[[split]][1]),
    pct = 100 * mean(s$expr > 0, na.rm = TRUE), mean = mean(s$expr, na.rm = TRUE))))
  agg$split <- factor(agg$split, levels = levels(factor(df[[split]])))
  ggplot(agg, aes(grp, split, size = pct, color = mean)) + geom_point() +
    scale_color_viridis_c(option = "magma", direction = -1) + scale_size_area(max_size = 12) +
    theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(x = labof(grp), y = labof(split), size = "% expr", color = "mean", title = gene)
}
# mirrors cm_markerheat_p (app.R:2973-2983)
markerheat_gg <- function() {
  h <- heat[["res0.2"]]; long <- h$long
  long$gene <- factor(long$gene, levels = rev(h$genes))
  long$cluster <- factor(long$cluster, levels = h$clusters)
  ggplot(long, aes(cluster, gene, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "subcluster", y = "marker gene", fill = "z-score\nmean expr",
         title = "Subcluster identity markers — res 0.2")
}
# mirrors cm_phase_plot (app.R:3070-3090)
cmphase_gg <- function() {
  df <- cmm; df$sub <- factor(paste0("CM", df[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
  tab <- table(df[c("sub", "Phase", "genotype", "timepoint")])
  tb  <- as.data.frame(prop.table(tab, setdiff(seq_along(dim(tab)), 2L)))
  names(tb)[match("Freq", names(tb))] <- "prop"
  tb$prop[is.nan(tb$prop)] <- 0
  tb$Phase <- factor(tb$Phase, levels = c("G1", "S", "G2M"))
  ggplot(tb, aes(sub, prop, fill = Phase)) + geom_col() +
    facet_grid(genotype ~ timepoint) +
    scale_fill_manual(values = setNames(c("#bdbdbd", "#1565c0", "#c62828"), c("G1","S","G2M"))) +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(x = "subcluster", y = "fraction of cells",
         title = "Cell-cycle phase by subcluster — res 0.2")
}
# mirrors e2f_plot (app.R:3362-3391)
e2f_gg <- function(ct = "Cardiomyocyte") {
  df <- meta[as.character(meta$celltype) == ct, ]
  eg <- intersect(c("E2f7", "E2f8"), c(rownames(expr), rownames(EXPR)))
  long <- do.call(rbind, lapply(eg, function(g) data.frame(
    gene = g, expr = expr_vec(g, df$cell), genotype = df$genotype, timepoint = df$timepoint)))
  long <- long[is.finite(long$expr), ]
  long$genotype <- factor(long$genotype, levels = levels(factor(meta$genotype)))
  ggplot(long, aes(genotype, expr, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, alpha = .55, linewidth = .2) +
    stat_summary(fun = mean, geom = "point", size = 2.4, color = "black") +
    facet_grid(gene ~ timepoint, scales = "free_y") +
    scale_fill_manual(values = setNames(c("#c62828", "#1565c0"), c("KO", "WT"))) +
    theme_pub + guides(fill = "none") +
    labs(x = "genotype", y = "log-norm expression",
         title = paste0("E2f7 / E2f8 — ", ct, " (black dot = mean; descriptive, n = 1)"))
}
# mirrors e2f_fc_plot (app.R:3401-3424)
e2ffc_gg <- function(set = "E2F targets", ct = "Cardiomyocyte") {
  tg <- intersect(GENE_SETS[[set]], rownames(expr))
  rows <- do.call(rbind, lapply(c("P0", "P7"), function(tp) {
    d <- ctDE[[paste(tp, ct, sep = "_")]]; if (is.null(d)) return(NULL)
    i <- match(tg, d$gene); i <- i[!is.na(i)]
    data.frame(gene = d$gene[i], timepoint = tp, lfc = d$log2FoldChange[i], padj = d$padj[i])
  }))
  rows$dir  <- ifelse(rows$lfc >= 0, "up in KO", "up in WT")
  rows$gene <- factor(rows$gene,
                      levels = unique(rows$gene[order(ave(rows$lfc, rows$gene, FUN = mean))]))
  ggplot(rows, aes(lfc, gene, color = dir)) +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_segment(aes(x = 0, xend = lfc, yend = gene), linewidth = .5) + geom_point(size = 2) +
    facet_wrap(~ timepoint) + scale_color_manual(values = VOLC_PAL, name = NULL) +
    theme_pub + labs(x = "log2 fold change (KO / WT)", y = NULL,
                     title = paste0(set, " — KO vs WT in ", ct))
}
# mirrors cyc_violins_plot (app.R:4205-4225)
cyc_gg <- function(ct = "Cardiomyocyte") {
  cols <- intersect(c("sig_prolif","sig_cytokinesis","sig_ccexit","sig_ploidy"), names(meta))
  df <- meta[as.character(meta$celltype) == ct, ]
  long <- do.call(rbind, lapply(cols, function(c) {
    d <- df[!is.na(df[[c]]), , drop = FALSE]
    data.frame(genotype = d$genotype, timepoint = as.character(d$timepoint),
               score = labof(c), value = d[[c]])
  }))
  long$score <- factor(long$score, levels = labof(cols))
  ggplot(long, aes(genotype, value, fill = timepoint)) +
    geom_violin(scale = "width", trim = TRUE, alpha = .85, linewidth = .2,
                position = position_dodge(.9)) +
    facet_wrap(~ score, scales = "free_y") + theme_pub +
    labs(x = NULL, y = "module score", fill = "timepoint",
         title = paste0("Cycle-exit / ploidy scores — ", ct))
}
# mirrors gm_plot_ly (app.R:1960-2018) as a static panel
genemap_gg <- function(panel = "avg", label_n = 24) {
  d <- gm_df(panel, hide_sets = TRUE)
  ctr <- gm_centre(panel)
  d$quadrant <- factor(d$quadrant, levels = GM_QUADS)
  lab <- head(d[order(-d$distance), ], label_n)
  p <- ggplot(d, aes(x, y, colour = quadrant)) +
    geom_hline(yintercept = ctr[["met"]], colour = "grey70") +
    geom_vline(xintercept = ctr[["mat"]], colour = "grey70") +
    geom_point(size = .7, alpha = .5) +
    scale_colour_manual(values = GM_QUAD_PAL, name = NULL) +
    theme_pub + labs(x = "maturation AUC (mature →)", y = "metabolic AUC (oxidative →)",
      title = sprintf("Gene map — %s panel (%s genes, set genes hidden)", panel,
                      format(nrow(d), big.mark = ",")),
      caption = sprintf("Axes split at each one's own median (mat %.3f, met %.3f), not 0.5.",
                        ctr[["mat"]], ctr[["met"]]))
  # seed= is load-bearing, not decoration: ggrepel starts each label from a random
  # offset and iterates, so without it this figure's pixels differ on every export and
  # genemap.png shows up dirty in git after a no-op re-run. That churn is worse than
  # cosmetic -- it makes "did this figure actually change?" unanswerable from the diff.
  if (requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.7,
                                      max.overlaps = 40, show.legend = FALSE, seed = 42)
  p
}

# ===========================================================================
# 01 -- upstream: the PIPseeker run itself
# ===========================================================================
# original_Han_analysis/ holds the upstream run we received: PIPseeker's per-lane
# run_config.csv, its metrics, and the STAR logs. These are the primary record of how
# the count matrices were made, so the chapter cites them rather than prose. The heavy
# parts of that folder are git-ignored; these text files are not.
cat("\n== 01 upstream provenance (original_Han_analysis) ==\n")
ORIG  <- file.path(REPO, "original_Han_analysis")
RAWPS <- file.path(ORIG, "raw_pipseeker_outputs")
LANES <- if (dir.exists(RAWPS)) sort(list.dirs(RAWPS, recursive = FALSE, full.names = FALSE)) else character(0)

kv <- function(path) {                       # PIPseeker writes key,value with no header
  if (!file.exists(path)) return(NULL)
  d <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("k", "v"))
  setNames(as.list(d$v), d$k)
}
metrics <- function(lane, sens) kv(file.path(RAWPS, lane, "metrics",
                                             paste0("sensitivity_", sens), "metrics_summary.csv"))
runcfg  <- function(lane) kv(file.path(RAWPS, lane, "run_config.csv"))

if (!length(LANES)) {
  gap("tbl-upstream-run", "01", "the PIPseeker run parameters and per-lane metrics",
      "original_Han_analysis/raw_pipseeker_outputs is not present in this working tree.")
  gap("tbl-sensitivity", "01", "the cell-calling sensitivity comparison", "same")
  gap("tbl-star", "01", "STAR alignment statistics", "same")
  gap("tbl-mito-source", "01", "mitochondrial fraction as PIPseeker measured it", "same")
} else {

  step("tbl-upstream-run", "01", "the PIPseeker run, per lane",
       frag("tbl-upstream-run", local({
         do.call(rbind, lapply(LANES, function(l) {
           c0 <- runcfg(l); m <- metrics(l, 5)
           data.frame(lane = l,
                      chemistry = c0$chemistry %||% NA,
                      pipseeker = c0$pipeline_version %||% NA,
                      # lanes whose run_config.csv was not included in the copy we
                      # received carry the command-line values from
                      # preprocessing_scripts/pipseeker_scripts.txt instead
                      # one lane's run_config.csv was not included in the copy we
                      # received; for it, report the range PIPseeker actually computed,
                      # which is what the metrics folders record
                      sensitivity = if (is.null(c0)) {
                          have <- gsub("sensitivity_", "", list.files(file.path(RAWPS, l, "metrics"),
                                                                      pattern = "^sensitivity_"))
                          if (length(have)) sprintf("%s-%s*", min(have), max(have)) else NA
                        } else sprintf("%s-%s", c0$min_sensitivity %||% "?",
                                                c0$max_sensitivity %||% "?"),
                      cells = as.integer(m$num_cell_barcodes),
                      mean_reads_per_cell = as.integer(m$mean_reads_per_cell),
                      median_genes = as.integer(m$median_genes_in_cells),
                      saturation_pct = as.numeric(m$sequencing_saturation),
                      stringsAsFactors = FALSE)
         }))
       }), caption = "The upstream run as PIPseeker recorded it. `sensitivity` is the range computed; the `filtered_matrix/sensitivity_5/` matrices are the ones everything downstream reads."))

  step("tbl-sensitivity", "01", "cells called at each sensitivity",
       frag("tbl-sensitivity", local({
         do.call(rbind, lapply(LANES, function(l) {
           r <- data.frame(lane = l, stringsAsFactors = FALSE)
           for (s in 3:5) {
             m <- metrics(l, s)
             r[[paste0("sens_", s)]] <- if (is.null(m)) NA_integer_ else as.integer(m$num_cell_barcodes)
           }
           r$used <- "sensitivity 5"
           r
         }))
       }), caption = "Cell barcodes called at each sensitivity PIPseeker was asked to compute. Blank means that level was not requested for that lane.",
          note = "Sensitivity is the single most consequential upstream choice: on P0KO_lane6, the only lane computed at all three levels, it moves the cell count from 8,271 to 12,045 -- a 46 % swing. The most permissive level available was the one used."))

  step("tbl-star", "01", "STAR alignment, per lane",
       frag("tbl-star", local({
         do.call(rbind, lapply(LANES, function(l) {
           f <- file.path(RAWPS, l, "star", "Log.final.out")
           if (!file.exists(f)) return(NULL)
           tx <- readLines(f, warn = FALSE)
           g <- function(pat) {
             h <- grep(pat, tx, value = TRUE, fixed = TRUE)
             if (!length(h)) return(NA_character_)
             trimws(sub(".*\\|", "", h[1]))
           }
           m <- metrics(l, 5)
           data.frame(lane = l,
                      input_reads = as.numeric(gsub("[^0-9]", "", g("Number of input reads"))),
                      uniquely_mapped_pct = as.numeric(gsub("%", "", g("Uniquely mapped reads %"))),
                      multi_loci_pct = as.numeric(gsub("%", "", g("% of reads mapped to multiple loci"))),
                      mapped_to_txome_pct = as.numeric(m$mapping_pct_txome),
                      stringsAsFactors = FALSE)
         }))
       }), caption = "STAR 2.7.6a, run inside PIPseeker, against Ensembl GRCm38 primary assembly with the release-99 GTF."))

  step("tbl-mito-source", "01", "mitochondrial fraction as measured upstream",
       frag("tbl-mito-source", local({
         d <- do.call(rbind, lapply(LANES, function(l) {
           m <- metrics(l, 5)
           data.frame(sample = sub("_lane[0-9]+$", "", l), lane = l,
                      pct_mito_in_cells = as.numeric(m$pct_mito_in_cells),
                      stringsAsFactors = FALSE)
         }))
         agg <- aggregate(pct_mito_in_cells ~ sample, d, mean)
         agg <- agg[match(c("P0WT","P0KO","P7WT","P7KO"), agg$sample), , drop = FALSE]
         agg <- agg[!is.na(agg$sample), , drop = FALSE]
         names(agg)[2] <- "pct_mito_mean_of_2_lanes"
         merge(agg, setNames(aggregate(pct_mito_in_cells ~ sample, d,
                 function(x) paste(sprintf("%.2f", x), collapse = ", ")),
                 c("sample", "per_lane")), by = "sample", sort = FALSE)
       }), digits = 3,
          caption = "PIPseeker's own per-lane mitochondrial estimate, before any of our filtering.",
          note = "The KO is *lower* than WT at P0 and *higher* at P7 -- the same direction and shape as the downstream mitochondrial confound, measured independently and upstream of everything we did. The two lanes of each sample agree to within 0.03 pp, which is what technical replicates should look like."))
}

# ===========================================================================
# 01 -- what our re-analysis changed, measured against the original run
# ===========================================================================
# The hand-off README says the original R scripts were edited for "QC correctness fixes --
# doublet removal, mouse mito QC, ..." and points at a METHODS_COMPARISON.md that is not in
# this tree. Rather than take that on trust, this compares the original run's own merged
# objects against the cells that survive into our bundle, cell by cell. Every claim about
# what changed is then a measurement.
cat("\n== 01 original vs ours ==\n")
ORIGPROC <- if (dir.exists("/orig")) "/orig/processing" else file.path(ORIG, "processing")
HAVE_SEURAT <- requireNamespace("SeuratObject", quietly = TRUE) &&
               dir.exists(ORIGPROC) && length(list.files(ORIGPROC, pattern = "\\.rds$"))

if (!HAVE_SEURAT) {
  gap("tbl-changed", "01", "the original-vs-ours comparison",
      paste("needs SeuratObject and original_Han_analysis/processing/*.rds.",
            "The heavy parts of that folder are git-ignored, so a fresh clone will not have them."))
} else {
  ORIG_CMP <- local({
    suppressMessages(library(SeuratObject))
    m <- meta; m$bc <- sub("^[^_]+_", "", m$cell)
    do.call(rbind, lapply(c("P0KO","P0WT","P7KO","P7WT"), function(sm) {
      f <- file.path(ORIGPROC, sprintf("merge.lanes.%s.rds", sm))
      if (!file.exists(f)) return(NULL)
      o  <- readRDS(f); md <- o@meta.data
      mt <- grep("^mt-", rownames(o), value = TRUE)
      # per layer: the two layers are the two lanes, and their cells are disjoint
      pm <- unlist(lapply(SeuratObject::Layers(o), function(lay) {
        x <- SeuratObject::LayerData(o, layer = lay)
        setNames(100 * Matrix::colSums(x[mt, , drop = FALSE]) /
                 pmax(Matrix::colSums(x), 1), colnames(x))
      }))
      data.frame(sample = sm, cell = rownames(md),
                 nFeature = md$nFeature_RNA, pct_mt = as.numeric(pm[rownames(md)]),
                 doublet = md$doublet %in% c(TRUE, "True", "TRUE"),
                 kept = rownames(md) %in% m$bc[m$orig.ident == sm],
                 stringsAsFactors = FALSE)
    }))
  })

  step("tbl-changed", "01", "what the re-analysis changed, per sample",
       frag("tbl-changed", local({
         d <- ORIG_CMP
         do.call(rbind, lapply(unique(d$sample), function(sm) {
           x <- d[d$sample == sm, ]
           elig <- sum(!x$doublet & x$pct_mt <= 20)
           # Columns are laid out so the arithmetic CLOSES left to right. An earlier
           # version showed only original / removals / bundle and left the downsample
           # implicit in a percentage, which read as ~8,500 cells going missing.
           data.frame(sample = sm,
                      original_cells = nrow(x),
                      minus_doublets = -sum(x$doublet),
                      minus_mito_over_20 = -sum(x$pct_mt > 20 & !x$doublet),
                      eligible = elig,
                      minus_not_bundled = -(elig - sum(x$kept)),
                      in_our_bundle = sum(x$kept),
                      stringsAsFactors = FALSE)
         }))
       }), caption = "The original run's merged objects against the cells that reach our bundle. Each row subtracts left to right and closes on the last column.",
          note = "Two different things are happening in this table and they should not be confused. The two negative QC columns are **filters** -- those cells are judged unusable. The `minus_not_bundled` column is the **browser downsample**: a stratified ~50 % sample taken so the app loads in a browser, not a quality judgement. Those cells still exist upstream, and the precomputed differential-expression tables were built before it. Every one of our 30,030 cells is present in the original objects, so both runs start from the identical PIPseeker matrices -- nothing was re-counted."))

  step("tbl-changed-detail", "01", "retention by doublet flag and mito bin",
       frag("tbl-changed-detail", local({
         d <- ORIG_CMP
         rows <- list()
         add <- function(group, level, i) rows[[length(rows) + 1L]] <<- data.frame(
           group = group, level = level, cells = length(i),
           retained = sum(d$kept[i]), retained_pct = round(100 * mean(d$kept[i]), 1),
           stringsAsFactors = FALSE)
         add("Scrublet doublet", "not a doublet", which(!d$doublet))
         add("Scrublet doublet", "flagged doublet", which(d$doublet))
         br <- c(-Inf, 5, 10, 15, 20, Inf)
         lb <- c("< 5 %", "5-10 %", "10-15 %", "15-20 %", "> 20 %")
         b  <- cut(d$pct_mt, br, labels = lb)
         for (lv in lb) { i <- which(b == lv & !d$doublet); if (length(i)) add("Mitochondrial %", lv, i) }
         cl <- !d$doublet & d$pct_mt <= 20
         fb <- cut(d$nFeature, c(-Inf, 7000, 8000, 9000, 10000, Inf),
                   labels = c("< 7,000", "7-8,000", "8-9,000", "9-10,000", "> 10,000"))
         for (lv in levels(fb)) { i <- which(fb == lv & cl); if (length(i)) add("Genes detected", lv, i) }
         do.call(rbind, rows)
       }), caption = "Retention of the original run's cells in our bundle, by the property being tested. A category retained at ~50 % was subject only to the downsample; one retained at 0 % was filtered out.",
          note = "Two filters are categorical and one is graded. Doublets and cells above 20 % mitochondrial are removed outright. Cells with very high gene counts are progressively depleted rather than cut at a fixed number -- consistent with a quantile- or MAD-based upper bound rather than a hard threshold; the highest-complexity cell we keep has 9,307 genes."))
}

# ===========================================================================
# 01 -- upstream: QC, normalization, annotation
# ===========================================================================
cat("\n== 01 upstream (from the bundle) ==\n")
step("qc-filtering",  "01", "upstream QC filtering figure",  figs_png("qc-filtering", "filtering"))
step("qc-violins",    "01", "per-lane QC violins",           figs_png("qc-violins",   "qc_violins"))
step("qc-doublets",   "01", "doublet detection figure",       figs_png("qc-doublets",  "doublet"))
step("qc-hvg",        "01", "HVG / SCTransform figure",       figs_png("qc-hvg",       "hvg"))
step("qc-harmony",    "01", "Harmony before/after",           figs_png("qc-harmony",   "harmony"))
step("tbl-doublet",   "01", "per-lane doublet rates",
     frag("tbl-doublet", tabs$doublet,
          caption = "Per-lane doublet calls carried in the bundle (`app$tables$doublet`)."))
step("refmap-heat",   "01", "annotation concordance heatmap", fig("refmap-heat", refmap_heat_gg(), 8, 5))
step("tbl-refmap-cov", "01", "reference panel coverage",
     frag("tbl-refmap-cov", REFMAP$coverage,
          caption = "Marker panels scored by `build_refmap.R`, and how many of each panel's genes exist in the curated matrix."))
step("tbl-refmap-conf", "01", "annotation confusion table",
     frag("tbl-refmap-conf", refmap_table_df(),
          caption = "Predicted (argmax over marker panels) vs assigned cell type, row-normalised."))
step("tbl-celltypes", "01", "cells per cell type x group",
     frag("tbl-celltypes",
          local({
            t <- as.data.frame.matrix(table(meta$celltype, paste(meta$genotype, meta$timepoint, sep = "-")))
            data.frame(celltype = rownames(t), t[, intersect(FG_GROUPS, names(t)), drop = FALSE],
                       total = rowSums(t), check.names = FALSE)
          }),
          caption = "Cells per cell type and four-group, as bundled (a stratified downsample of the upstream object)."))

# ===========================================================================
# 02 -- the cell atlas
# ===========================================================================
cat("\n== 02 atlas ==\n")
step("umap-celltype", "02", "UMAP coloured by cell type",
     fig("umap-celltype", umap_gg(meta, "celltype", "Cell type"), 7.5, 5.5, raster = TRUE))
step("umap-phase",    "02", "UMAP coloured by cell-cycle phase",
     fig("umap-phase", umap_gg(meta, "Phase", "Cell-cycle phase", labels = FALSE), 7.5, 5.5, raster = TRUE))
step("umap-gene",     "02", "UMAP coloured by Myh6 expression",
     fig("umap-gene", local({ d <- meta; d$v <- expr_vec("Myh6", d$cell)
                              umap_gg(d, "v", "Myh6 (log-norm)", continuous = TRUE) }), 7.5, 5.5, raster = TRUE))
step("comp-celltype", "02", "cell-type composition by group",
     fig("comp-celltype", comp_gg("splitgrp", "celltype"), 7, 4.6))
step("gene-violin",   "02", "gene-detail violin",
     fig("gene-violin", vln_gg("Myh6", "celltype", "genotype"), 8, 4.6))
step("gene-dot",      "02", "gene-detail dot plot",
     fig("gene-dot", dot_gg("Myh6", "celltype", "genotype"), 7.5, 3.6))
step("tbl-matrices",  "02", "the two expression matrices",
     frag("tbl-matrices", data.frame(
       matrix = c("`app$expr` (curated panel)", "`app$deg_expr` (broad)"),
       genes  = c(nrow(expr), nrow(EXPR)),
       cells  = c(ncol(expr), ncol(EXPR)),
       used_for = c("every live plot: UMAP colour, violins, dot plots, module scores",
                    "on-the-fly DE, and any gene absent from the curated panel")),
       caption = "The two matrices in the bundle, and which questions each answers."))

# ===========================================================================
# 03 -- differential expression
# ===========================================================================
cat("\n== 03 differential expression ==\n")
CT_KEY <- "P7_Cardiomyocyte"
step("de-volcano",  "03", "pooled KO-vs-WT volcano (static form of the app's plotly panel)",
     fig("de-volcano", de_volcano(ctDE[[CT_KEY]], "P7 Cardiomyocyte — KO vs WT"), 7, 5, raster = TRUE))
step("de-lfcheat",  "03", "log2FC heatmap, top genes x cell type",
     fig("de-lfcheat", lfc_heat(ctDE[grep("^P7_", names(ctDE))], topn = 22,
                                ttl = "P7 — top genes by |log2FC| across cell types"), 7.5, 6.5))
step("tbl-de-head", "03", "top rows of a DE table, as the app shows them",
     frag("tbl-de-head", head(de_table(ctDE[[CT_KEY]]), 12), max_rows = 12,
          caption = "First rows of `app$tables$ct_DE[[\"P7_Cardiomyocyte\"]]` after `de_table()` formatting."))
step("tbl-de-sizes", "03", "how many DE tables, and how big",
     frag("tbl-de-sizes", local({
       k <- names(ctDE)
       data.frame(table = k,
                  genes_kept = vapply(ctDE, nrow, 0L),
                  confounder_rows = vapply(ctDE, function(d) sum(d$confounder), 0L),
                  row.names = NULL)[order(k), ]
     }), caption = "Every precomputed cell-type DE table in the bundle."))

# ===========================================================================
# 04 -- pathway enrichment
# ===========================================================================
cat("\n== 04 enrichment ==\n")
step("enr-go-dot",  "04", "GO BP dot plot, cell-type level",
     fig("enr-go-dot", go_dotplot_gg(enr_go("Cardiomyocyte", "P7"),
           "GO BP enriched in KO-up genes — Cardiomyocyte P7"), 8, 5.5))
step("enr-gsea-bar", "04", "GSEA bar plot, cell-type level",
     fig("enr-gsea-bar", gsea_barplot_gg(enr_gsea("Cardiomyocyte", "P7"),
           "GSEA — Cardiomyocyte P7"), 8, 5.5))
step("enr-e2f-heat", "04", "E2F-family regulon activity heatmap",
     fig("enr-e2f-heat", enr_e2f_heat_gg(), 8, 4.5))
step("enr-tf-top",   "04", "top TFs by |KO-WT| activity",
     fig("enr-tf-top", enr_tf_top_gg("Cardiomyocyte"), 7, 5.5))
step("tbl-enr-regimes", "04", "the two enrichment threshold regimes",
     frag("tbl-enr-regimes", data.frame(
       builder = c("build_subcluster_enrichment.R", "build_fourgroup_enrichment.R",
                   "analysis/2026-08-21_email/02_enrich.R"),
       question = c("KO vs WT, pooled over P0+P7, per CM subcluster",
                    "each four-group contrast, per subcluster and stratum",
                    "the two contrasts sent to the collaborator"),
       input_gate = c("padj < 0.05 & |log2FC| >= 1", "padj < 0.05 & |log2FC| >= 0.25",
                      "padj < 0.05 & |log2FC| >= 0.25"),
       GO_p_q = c("0.2 / 0.2", "0.05 / 0.2", "0.05 / 0.2"),
       ontologies = c("BP", "BP, MF, CC", "BP, MF, CC"),
       universe = c("genes in that cluster's DE table",
                    "genes detected in >= 5% of one arm (~11k)",
                    "genes with max(pct) >= 5% in that table")),
       caption = "Three enrichment runs with deliberately different thresholds. Conflating them is the easiest way to misread the app.",
       note = "The permissive 0.2/0.2 regime is an exploratory screen over small per-subcluster lists; the 0.05/0.2 regime is what the contrast-specific results and the collaborator deliverable use."))

# ===========================================================================
# 05 -- cardiomyocyte deep-dive
# ===========================================================================
cat("\n== 05 CM deep-dive ==\n")
step("cm-umap",      "05", "CM subcluster UMAP, res 0.2",
     fig("cm-umap", local({
       d <- cmm; d$sub <- factor(paste0("CM", d[["SCT_snn_res.0.2"]]), levels = cm_subs("0.2"))
       umap_gg(d, "sub", "CM subclusters — res 0.2", psize = .3)
     }), 7.5, 5.5, raster = TRUE))
step("cm-markerheat", "05", "subcluster identity marker heatmap", fig("cm-markerheat", markerheat_gg(), 7.5, 8))
step("cm-phase",      "05", "phase composition by subcluster",    fig("cm-phase", cmphase_gg(), 8.5, 5))
step("cm-lfcheat",    "05", "pooled KO-vs-WT log2FC across subclusters",
     fig("cm-lfcheat", lfc_heat(subDE[["res0.2"]], topn = 22,
           ttl = "Pooled KO vs WT — top genes across CM subclusters"), 7.5, 6.5))
step("tbl-cm-summary", "05", "per-subcluster summary sheet",
     frag("tbl-cm-summary", subSum[["res0.2"]], digits = 3,
          caption = "`app$tables$sub_summary[[\"res0.2\"]]` — one row per CM subcluster."))
step("tbl-cm-subtype", "05", "nearest CM subtype per subcluster",
     frag("tbl-cm-subtype", subType[["res0.2"]],
          caption = "Each res-0.2 subcluster's nearest curated CM subtype."))

# ===========================================================================
# 06 -- E2F focus
# ===========================================================================
cat("\n== 06 E2F ==\n")
step("e2f-violin", "06", "E2f7/E2f8 expression by genotype x timepoint",
     fig("e2f-violin", e2f_gg(), 7, 5))
step("e2f-targets", "06", "E2F target log2FC lollipop",
     fig("e2f-targets", e2ffc_gg(), 7.5, 6.5))
step("tbl-e2f", "06", "E2f7/E2f8 detection by group",
     frag("tbl-e2f", local({
       do.call(rbind, lapply(c("E2f7", "E2f8"), function(g) {
         v <- expr_vec(g, meta$cell)
         do.call(rbind, lapply(FG_GROUPS, function(gr) {
           i <- paste(meta$genotype, meta$timepoint, sep = "-") == gr & !is.na(v)
           data.frame(gene = g, group = gr, n_cells = sum(i),
                      pct_detected = 100 * mean(v[i] > 0), mean_lognorm = mean(v[i]))
         }))
       }))
     }), caption = "The knockout is not visible in the transcript: detection and mean expression of E2f7/E2f8 by group.",
        note = "This is the single most important negative control in the dataset. See the chapter text for why a conditional allele plus a 3'-biased assay can produce it."))

# ===========================================================================
# 07 -- four-group design
# ===========================================================================
cat("\n== 07 four-group ==\n")
step("fg-counts",  "07", "four-group sizes per subcluster",
     fig("fg-counts", fg_counts_plot("prop"), 8.5, 5))
step("fg-phase",   "07", "phase composition, four groups",
     fig("fg-phase", fg_phase_plot(c("AllCM", "CM1", "CM2", "CM3", "CM4", "CM5")), 8.5, 5))
step("fg-score",   "07", "maturation score by four group",
     fig("fg-score", fg_score_plot("sig_maturation_nocc",
           c("AllCM", "CM1", "CM2", "CM3", "CM7", "CM8"), "G1"), 8.5, 5))
step("fg-summary", "07", "KO-WT maturation gap summary",
     fig("fg-summary", fg_summary_plot("sig_maturation_nocc", "G1"), 8, 5))
step("fg-volcano", "07", "P7 KO-vs-WT volcano, CM2, G1 stratum",
     fig("fg-volcano", local({
       ct <- fg_ct("P7_KO_vs_WT")
       d  <- fg_de("CM2", "P7_KO_vs_WT", "G1")
       de_volcano(d, "CM2 — P7 KO vs WT (G1 stratum)")
     }), 7, 5, raster = TRUE))
step("tbl-fg-counts", "07", "four-group counts table",
     frag("tbl-fg-counts", fg_counts_wide(),
          caption = "Cells per subcluster x group, with the underpowered flag the app shows."))
step("tbl-fg-contrasts", "07", "the four contrasts, as defined",
     frag("tbl-fg-contrasts", FG$built$contrasts,
          caption = "`app$fourgroup$built$contrasts` — every contrast the four-group grids carry, with its direction labels."))
step("tbl-fg-params", "07", "build_fourgroup.R parameters as built",
     frag("tbl-fg-params", local({
       b <- FG$built
       flat <- function(x) paste(as.character(x), collapse = ", ")
       keys <- intersect(c("groups","resolution","min_cells","thin_cells","n_genes","n_cells_de",
                           "n_genes2","n_cells_total","matrix","mat_auc","pct_min","lfc_min"), names(b))
       data.frame(parameter = keys, value = vapply(keys, function(k) flat(b[[k]]), ""))
     }), caption = "Recorded in the bundle by the builder, not typed here."))
step("tbl-fg-coverage", "07", "contrast coverage across the two grids",
     frag("tbl-fg-coverage", local({
       cl <- FG_CLUSTERS
       rows <- do.call(rbind, lapply(cl, function(c) {
         k1 <- names(FG$de[[c]]); k2 <- names(FG$de2[[c]])
         data.frame(cluster = c, de_broad = length(k1), de2_curated = length(k2),
                    only_in_curated = length(setdiff(k2, k1)),
                    only_in_broad = length(setdiff(k1, k2)))
       }))
       rows
     }), caption = "Why the app offers two DE matrices: neither grid covers every contrast."))
step("tbl-fg-skipped", "07", "contrasts that were not computed, and why",
     frag("tbl-fg-skipped", FG$skipped, max_rows = 40,
          caption = "`app$fourgroup$skipped` — each arm that fell below the floor, with its cell counts."))
HAVE_FGE <- !is.null(FGE) && !is.null(FGE$go)
FGE_WHY  <- paste("app$enrich$fourgroup is not in this bundle -- build_fourgroup_enrichment.R",
                  "has not been run against it (~2 h, ~460 enrichGO calls). The app shows its",
                  "\"run the builder\" message on that panel for the same reason.")
if (HAVE_FGE) {
  step("fg-enr-go", "07", "four-group GO, P7 KO vs WT",
       fig("fg-enr-go", go_dotplot_gg(
         # "A_up", not "up": the builder writes A_up/B_up and the app translates
         # its own up/down into those (app.R:1019-1029). Passing "up" here matches
         # no rows and yields an empty panel rather than an error.
         fg_enr_df("go", "CM2", "P7_KO_vs_WT", "G1", "BP", "A_up"),
         "GO BP — up in P7 KO, CM2 (G1)"), 8, 5.5))
  step("fg-enr-gsea", "07", "four-group GSEA, P7 KO vs WT",
       fig("fg-enr-gsea", local({
         ct <- fg_ct("P7_KO_vs_WT")
         gsea_barplot_gg(fg_enr_df("gsea", "CM2", "P7_KO_vs_WT", "G1"),
                         "GSEA — P7 KO vs WT, CM2 (G1)",
                         up_lab = ct$pos, down_lab = ct$neg)
       }), 8, 5.5))
  step("tbl-fg-audit", "07", "enrichment coverage audit",
       frag("tbl-fg-audit", local({
         a <- FGE$audit
         a[a$cluster %in% c("AllCM","CM1","CM2") & a$ontology == "BP", , drop = FALSE]
       }), max_rows = 30,
          caption = "The audit row behind every GO panel: genes in, universe size, and which fallback rule selected the list."))
} else {
  gap("fg-enr-go",   "07", "GO/GSEA figures for the four-group contrasts", FGE_WHY)
  gap("tbl-fg-audit", "07", "the enrichment coverage audit table", FGE_WHY)
}

# ===========================================================================
# 08 -- module scores and signature axes
# ===========================================================================
cat("\n== 08 scores ==\n")
step("score-coverage", "08", "module-score definitions and coverage",
     frag("score-coverage", SCOREMETA,
          caption = "`app$score_meta` — every score, the sets behind it, which matrix it was computed on, and how many of its genes were found."))
step("cyc-violins", "08", "cycle-exit / ploidy score violins", fig("cyc-violins", cyc_gg(), 8, 5.5))
step("ploidy-scatter", "08", "proliferation vs cytokinesis",
     fig("ploidy-scatter", ploidy_scatter("Cardiomyocyte"), 8, 4.6, raster = TRUE))
step("mat-scatter", "08", "transcriptional vs metabolic maturation",
     fig("mat-scatter", mat_scatter("Cardiomyocyte", stratum = "all"), 8.5, 5, raster = TRUE))
step("mat-violin", "08", "maturation score by four group",
     fig("mat-violin", score_violin("sig_maturation_nocc", "Cardiomyocyte", stratum = "G1"), 7, 5))
step("genemap", "08", "gene map, averaged panel", fig("genemap", genemap_gg("avg"), 8, 6, raster = TRUE))
step("tbl-genemap-quads", "08", "gene-map quadrant occupancy",
     frag("tbl-genemap-quads", local({
       d <- gm_df("avg", hide_sets = TRUE)
       t <- as.data.frame(table(quadrant = d$quadrant), stringsAsFactors = FALSE)
       t$pct <- round(100 * t$Freq / sum(t$Freq), 1)
       names(t)[2] <- "n_genes"
       t$diagonal <- ifelse(t$quadrant %in% c("mature+oxidative", "immature+glycolytic"),
                            "expected coupling", "uncoupled")
       t
     }), caption = "Quadrant occupancy with score-set genes hidden. The two diagonal quadrants summing above 50% is the coupling result."))
step("tbl-genemap-centres", "08", "the per-panel axis centres",
     frag("tbl-genemap-centres", local({
       do.call(rbind, lapply(GM_PANELS, function(p) {
         c0 <- gm_centre(p)
         data.frame(panel = p, mat_centre = c0[["mat"]], met_centre = c0[["met"]],
                    offset_from_0.5_mat = c0[["mat"]] - 0.5,
                    offset_from_0.5_met = c0[["met"]] - 0.5)
       }))
     }), digits = 4,
        caption = "Each panel's own median, and how far it sits from a naive 0.5 split."))

# ===========================================================================
# 09 -- set overlaps
# ===========================================================================
cat("\n== 09 overlaps ==\n")
step("fg-quadrant", "09", "maturation axis x P7 KO log2FC quadrant map",
     fig("fg-quadrant", fg_quadrant_plot("CM2"), 8, 5.5, raster = TRUE))
step("tbl-intersect", "09", "the intersection table",
     frag("tbl-intersect", head(fg_intersect_df("CM2"), 20), max_rows = 20,
          caption = "Top rows of the CM2 intersection: gene-level maturation association crossed with its P7 KO response."))
# AllCM rather than a single subcluster: at these thresholds a per-subcluster
# down-list is a handful of genes, and a Venn over one gene teaches nothing about
# the method. The overlap statistics are the point of the tab, and they need sets
# large enough for an expectation to be meaningful.
# The tab's own defaults (app.R:2533-2541): the two temporal contrasts crossed
# against the canonical cell-cycle panel, all cardiomyocytes, phase-matched,
# |log2FC| >= 0.25. Documenting the default view is the point -- it is what a
# reader sees before touching a control.
VN_P   <- list(cluster = "AllCM", stratum = "G1", grid = "de", padj = 0.05, lfc = 0.25)
VN_IDS <- c("de:WT_P0_vs_P7:both", "de:KO_P0_vs_P7:both", "cur:__canonical__")
step("venn", "09", "gene-set Venn",
     fig("venn", local({
       vs <- vn_sets(VN_IDS, VN_P$cluster, VN_P$stratum, VN_P$grid, VN_P$padj, VN_P$lfc)
       vn_plot(vs, ttl = "All cardiomyocytes, G1 stratum")
     }), 6.5, 5.5))
step("tbl-venn-stats", "09", "overlap statistics for the Venn",
     frag("tbl-venn-stats", local({
       vs <- vn_sets(VN_IDS, VN_P$cluster, VN_P$stratum, VN_P$grid, VN_P$padj, VN_P$lfc)
       vn_stats(vs)
     }), caption = "Every pairwise overlap against its hypergeometric expectation. The picture cannot carry this; the table has to."))
XC_P <- list(wt_cluster = "AllCM", mat_clusters = XC_MAT_CLUSTERS, cyc_clusters = XC_CYC_CLUSTERS,
             mat_set = "CM maturation", cyc_set = XC_CANON, stratum = "all", grid = "de",
             padj = 0.05, measure = "auc", eff = 0.60, eff_auc = 0.60, eff_lfc = 0.25,
             minc = 1, hide_mt = TRUE)
step("xc-venn", "09", "one of the four WT-program x KO-cluster crossings",
     fig("xc-venn", local({
       cmp <- XC_COMPARISONS[[1]]
       vn_plot(xc_comparison(cmp, XC_P), ttl = xc_label(cmp, XC_P))
     }), 6.5, 5))
step("tbl-xc-audit", "09", "AUC vs log2FC audit for the WT lists",
     frag("tbl-xc-audit", xc_measure_audit(XC_P),
          caption = "The same four WT P0->P7 lists under both effect-size measures. A symmetric |log2FC| cut cannot classify a cell-cycle gene on this data, which is why the tab defaults to AUC."))
step("tbl-xc-stats", "09", "overlap statistics for the four crossings",
     frag("tbl-xc-stats", local({
       do.call(rbind, lapply(XC_COMPARISONS, function(cmp) {
         s <- vn_stats(xc_comparison(cmp, XC_P))
         cbind(comparison = cmp$key, s)
       }))
     }), caption = "All four crossings, with fold enrichment and hypergeometric p."))

# ===========================================================================
# 10 -- cell-cell signalling
# ===========================================================================
cat("\n== 10 signalling ==\n")
step("commun-heat", "10", "sender x receiver signalling heatmap",
     fig("commun-heat", commun_heat_gg("VEGF", "P7", "delta"), 7.5, 5))
step("tbl-commun-pairs", "10", "the curated ligand-receptor pairs",
     frag("tbl-commun-pairs", COMMUN$pairs,
          caption = "Every ligand->receptor pair scored by `build_communication.R`. Curated, not database-derived."))
step("tbl-commun-scores", "10", "top signalling deltas",
     frag("tbl-commun-scores", head(commun_table_df("VEGF", "P7"), 15), max_rows = 15,
          caption = "VEGF at P7: sender, receiver, WT and KO scores, and the KO-WT delta."))

# ===========================================================================
# 11 -- confounds, sensitivity, reproducibility
# ===========================================================================
cat("\n== 11 confounds and reproducibility ==\n")
step("tbl-confound", "11", "the flagged sex/construct genes",
     frag("tbl-confound", data.frame(gene = CONF,
            role = "flagged in every table, excluded from every enrichment input"),
          caption = "`app$confound` — the genes that separate the two animals rather than the two genotypes."))
step("tbl-mito", "11", "mitochondrial read-fraction sensitivity",
     frag("tbl-mito", local({
       f <- file.path(DELIV, "csv/_mito_fraction.csv")
       if (!file.exists(f)) stop("run analysis/2026-08-21_email/run.sh first (", f, " missing)")
       read.csv(f, stringsAsFactors = FALSE)
     }), caption = "Per-arm mitochondrial share of log-normalised signal, from `analysis/2026-08-21_email/05_mito_sensitivity.R`."))
step("tbl-mito-go", "11", "GO terms surviving mt- removal",
     frag("tbl-mito-go", local({
       f <- file.path(DELIV, "csv/_mito_go_comparison.csv")
       if (!file.exists(f)) stop("missing ", f)
       read.csv(f, stringsAsFactors = FALSE)
     }), max_rows = 30,
        caption = "GO BP term counts with and without mitochondrial genes, per cluster and direction."))
step("tbl-stamp", "11", "provenance stamp",
     frag("tbl-stamp", data.frame(
       item = c("bundle built", "bundle file", "cells bundled", "CM cells bundled",
                "curated panel genes", "broad matrix genes", "broad matrix cells",
                "cell types", "CM subclusters (res 0.2)", "four-group DE tables (broad)",
                "four-group DE tables (curated)", "gene-map genes", "exported"),
       value = c(as.character(app$built), basename(BUNDLE),
                 format(nrow(meta), big.mark = ","), format(nrow(cmm), big.mark = ","),
                 format(length(genes), big.mark = ","),
                 format(length(GENES_FULL), big.mark = ","), format(ncol(EXPR), big.mark = ","),
                 length(unique(as.character(meta$celltype))), length(cm_subs("0.2")),
                 sum(vapply(FG$de, length, 0L)), sum(vapply(FG$de2, length, 0L)),
                 format(nrow(GM), big.mark = ","), format(Sys.Date()))),
       caption = "What the figures and tables in this book were generated from."))

# ---- manifest and exit -----------------------------------------------------
if (LIST) quit(save = "no", status = 0)

# _manifest.md and _gaps.md describe the WHOLE book, so only a full run can write a
# correct one. They used to be overwritten on every run including --only=, which meant a
# partial export silently truncated the manifest to just the steps that ran (it was found
# at 6 rows against 71 in git) and rewrote _gaps.md to claim "nothing missing" on the
# strength of a run that had checked almost nothing. A wrong index is worse than a stale
# one, because the stale one at least looks stale. Leave both alone unless this was a
# full run.
PARTIAL <- length(ONLY) > 0

G <- do.call(rbind, GAPS)
if (PARTIAL) {
  cat("\n-- partial run (--only=", ONLY, "): _manifest.md and _gaps.md left untouched;",
      "\n   re-run docs/export.sh with no --only to refresh them.\n", sep = "")
} else {
  if (!is.null(G)) {
    writeLines(c("<!-- generated by docs/export_assets.R -- do not edit -->", "",
                 md_table(G)), file.path(OUT_G, "_gaps.md"))
  } else {
    writeLines(c("<!-- generated by docs/export_assets.R -- do not edit -->", "",
                 "*Nothing missing: every figure and table in this book was exported",
                 "from the current bundle.*"), file.path(OUT_G, "_gaps.md"))
  }
}

M <- do.call(rbind, MANIFEST)
if (!is.null(M) && !PARTIAL) {
  writeLines(c("<!-- generated by docs/export_assets.R -- do not edit -->", "",
               md_table(M[order(M$chapter, M$id), c("chapter","id","what","file")])),
             file.path(OUT_G, "_manifest.md"))
}

cat(sprintf("\n== %d exported, %d failed, %d gap/skipped ==\n",
            if (is.null(M)) 0L else nrow(M), length(FAILED), length(SKIPPED)))
if (!is.null(G)) cat(sprintf("%d documented gap(s) written to _generated/_gaps.md\n", nrow(G)))
if (length(FAILED)) {
  cat("\nFAILURES:\n"); cat(paste0("  - ", FAILED, collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}
cat("Done.\n")
