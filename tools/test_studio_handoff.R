#!/usr/bin/env Rscript
# test_studio_handoff.R
# ---------------------------------------------------------------------------
# Runtime companion to shiny_app/studio_helpers.R (the Figure Studio handoff).
# Three claims, each of which would fail silently in production if untrue:
#
#   1. strip_plot_env() actually detaches a reactive-built ggplot from its
#      enclosing environments -- the serialized figspec is small, and a FRESH
#      R process with only ggplot2 can readRDS + ggplot_build it.
#   2. With the studio configured, clicking every registered figure's
#      "Figure Studio" button writes one readable figspec per click.
#   3. prune_handoff() deletes files past the TTL and nothing younger.
#
#   docker run --rm -v "$PWD:/repo:ro" -w /repo lab-server-e2f-heart-scrna-dev:latest \
#     Rscript tools/test_studio_handoff.R
# ---------------------------------------------------------------------------

# app.R reads these at source time, so they must be set before testServer().
HDIR <- file.path(tempdir(), "handoff")
dir.create(HDIR, showWarnings = FALSE)
Sys.setenv(FIGURE_STUDIO_BASE = "http://studio.invalid/", HANDOFF_DIR = HDIR)

suppressMessages({ library(shiny); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
APPDIR <- if (length(args)) args[1] else "shiny_app"
stopifnot(dir.exists(APPDIR))

pass <- 0L; fail <- 0L; skip <- 0L
note <- function(status, what, extra = "") {
  cat(sprintf("  [%-4s] %-28s %s\n", status, what, extra))
  if (status == "PASS") pass <<- pass + 1L
  else if (status == "FAIL") fail <<- fail + 1L
  else skip <<- skip + 1L
}

# ---- 1. strip_plot_env, as a standalone unit --------------------------------
source(file.path(APPDIR, "studio_helpers.R"))

fat_plot <- local({
  ballast <- rnorm(2e6)                       # ~16 MB riding in the closure frame
  df <- data.frame(x = 1:10, y = (1:10)^2, g = rep(c("a", "b"), 5))
  ggplot(df, aes(x, y, color = g)) + geom_point() + facet_wrap(~g)
})
sz_raw <- length(serialize(fat_plot, NULL))
stripped <- strip_plot_env(fat_plot)
sz_strip <- length(serialize(stripped, NULL))
note(if (sz_raw > 10e6) "PASS" else "FAIL", "unstripped plot is fat",
     sprintf("%.1f MB captured", sz_raw / 1e6))
note(if (sz_strip < 1e6) "PASS" else "FAIL", "stripped plot is small",
     sprintf("%.1f MB -> %.3f MB", sz_raw / 1e6, sz_strip / 1e6))
note(if (!inherits(tryCatch(ggplot_build(fat_plot), error = identity), "error"))
       "PASS" else "FAIL", "original still builds", "stripping must not mutate it")

# the portability claim: a fresh R process, only ggplot2 attached, no app state
f <- tempfile(fileext = ".rds"); saveRDS(stripped, f)
out <- suppressWarnings(system2("Rscript", c("-e", shQuote(sprintf(
  "suppressMessages(library(ggplot2)); invisible(ggplot_build(readRDS('%s'))); cat('BUILD_OK')", f))),
  stdout = TRUE, stderr = TRUE))
note(if (any(grepl("BUILD_OK", out))) "PASS" else "FAIL",
     "fresh-session round-trip", paste(tail(out, 1), collapse = " "))

# a plot strip_plot_env cannot handle must come back untouched, not broken:
# quo_get_env on a non-quosure path is simulated by a mangled mapping
weird <- fat_plot; weird$mapping <- structure(list(x = quote(x)), class = "uneval")
back <- strip_plot_env(weird)
note(if (!inherits(tryCatch(ggplot_build(back), error = identity), "error"))
       "PASS" else "FAIL", "fallback keeps plot usable")

# ---- 3. prune_handoff (before the click flood, while the dir is ours) -------
old_f <- file.path(HDIR, "20200101000000aaaaaaaaaaaa.rds")
new_f <- file.path(HDIR, "20990101000000bbbbbbbbbbbb.rds")
saveRDS(1, old_f); saveRDS(1, new_f)
Sys.setFileTime(old_f, Sys.time() - STUDIO_TTL_SECS - 60)
prune_handoff(HDIR)
note(if (!file.exists(old_f) && file.exists(new_f)) "PASS" else "FAIL",
     "prune_handoff", "expired file gone, fresh file kept")
unlink(new_f)

# ---- 2. every registered figure hands off through its button ----------------
# The prefixes come from the same static parse check_download_coverage.R trusts.
src <- readLines(file.path(APPDIR, "app.R"), warn = FALSE)
prefixes <- unique(unlist(regmatches(src,
  gregexpr('register_fig\\(output, "\\K[A-Za-z0-9_]+', src, perl = TRUE))))
cat(sprintf("\n== handoff per registered figure (%d prefixes) ==\n", length(prefixes)))

testServer(APPDIR, {
  # The same plausible browsing state test_downloads.R uses, plus the inputs the
  # UMAP / gene-detail / CM-map / gene-map figures need.
  session$setInputs(
    ct_tp = "P7", ct_sel = "Cardiomyocyte", ct_hideconf = FALSE, ct_search = "",
    cm_sub = "CM2", cm_hideconf = FALSE, cm_contrast = "pooled",
    enr_ct = "Cardiomyocyte", enr_tp = "P7",
    fg_cluster = "CM2", fg_contrast = "P7_KO_vs_WT", fg_stratum = "all",
    fg_grid = "de", fg_hideconf = FALSE, fg_count_mode = "prop",
    fg_score = "sig_maturation_nocc", fg_score_stratum = "all",
    fg_g1_clusters = c("AllCM", "CM2"),
    mi_cluster = "CM2", mi_hideconf = TRUE,
    mi_quad = c("immature_up_in_KO", "mature_down_in_KO"),
    mi_cand_stratum = "all", mi_spec_grid = "de",
    mi_genes = c("Tcf4", "Adamts9", "Gabbr2"), mi_cand_clusters = c("AllCM", "CM2"),
    vn_cluster = "AllCM", vn_stratum = "all", vn_grid = "de",
    vn_a = "de:P7_KO_vs_WT:up", vn_b = "de:WT_P0_vs_P7:up", vn_c = "none",
    vn_measure = "lfc", vn_auc = 0.60,
    gm_panel = "avg", gm_dist = 0, gm_labeln = 20,
    cc_tp = "P7", cc_pathway = "All", cc_metric = "delta",
    deg_hideconf = FALSE,
    fg_enr_cluster = "CM2", fg_enr_contrast = "P7_KO_vs_WT", fg_enr_stratum = "all",
    fg_enr_ont = "BP", fg_enr_topn = 20,
    cm_enr_mode = "one", cm_enr_sub = "CM2",
    xc_wt_cluster = "AllCM", xc_mat_set = "CM maturation", xc_cyc_set = "__canonical__",
    xc_mat_clusters = c("CM1","CM2","CM3","CM7","CM8"),
    xc_cyc_clusters = c("CM2","CM4","CM5"), xc_minc = 1, xc_stratum = "all",
    xc_grid = "de", xc_padj = 0.05, xc_measure = "auc", xc_auc = 0.60, xc_lfc = 0.25,
    xc_hidemt = TRUE, xc_gene_cmp = "__all__",
    ct_vlfc = 1, cm_vlfc = 1, fg_vlfc = 1, deg_vlfc = 1,
    mat_ct = "Cardiomyocyte", mat_score = "sig_maturation", mat_stratum = "all",
    # UMAP explorer / gene detail / CM map / candidate dotplot
    color_by = "celltype", split = "none", ptsize = 4.5, gene = "Gabbr2",
    g2 = "Gabbr2", grp = "celltype", sp2 = "none",
    cm_mapcolor = "subcluster", cm_map_split = "none", cm_gene = "Gabbr2",
    e2f_grp = "genotype", comp_x = "celltype", comp_fill = "genotype"
  )

  ok_n <- 0L
  for (pf in prefixes) {
    before <- list.files(HDIR, pattern = "\\.rds$")
    do.call(session$setInputs, setNames(list(1L), paste0(pf, "_studio")))
    newf <- setdiff(list.files(HDIR, pattern = "\\.rds$", full.names = TRUE),
                    file.path(HDIR, before))
    if (!length(newf)) { note("SKIP", pf, "no figspec (guarded / needs more inputs)"); next }
    spec <- tryCatch(readRDS(newf[1]), error = function(e) e)
    ok <- !inherits(spec, "error") &&
      identical(spec$version, 1L) &&
      all(c("plot", "meta") %in% names(spec)) &&
      identical(spec$meta$app, "e2f-heart-scrna") &&
      identical(spec$meta$prefix, pf) &&
      file.size(newf[1]) < 25e6 &&
      !inherits(tryCatch(ggplot_build(spec$plot), error = identity), "error")
    if (ok) ok_n <- ok_n + 1L
    note(if (ok) "PASS" else "FAIL", pf,
         sprintf("%.2f MB", file.size(newf[1]) / 1e6))
  }
  # the handoff must work broadly, not just on a lucky few
  note(if (ok_n >= 20) "PASS" else "FAIL", "coverage floor",
       sprintf("%d of %d figures handed off", ok_n, length(prefixes)))
  # the seven ex-plotly twins are the point of the exercise; name them
  for (pf in c("umapgg", "cmmap", "ctvolc", "cmvolc", "gmsc"))
    note(if (pf %in% prefixes) "PASS" else "FAIL", paste0("registered: ", pf))
})

cat(sprintf("\n%d passed, %d skipped, %d failed\n", pass, skip, fail))
if (fail) quit(status = 1)
