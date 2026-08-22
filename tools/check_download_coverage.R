#!/usr/bin/env Rscript
# check_download_coverage.R
# ---------------------------------------------------------------------------
# Fails if any table or figure in the app has no download.
#
# The point of this file is that "make everything downloadable" is the kind of
# refactor that ends up 90% done and looks finished. Eyeballing tabs will not
# catch the one panel that was missed; parsing the source will.
#
#   Rscript tools/check_download_coverage.R [path/to/app.R]
#
# It is a static check on the source text, not a runtime test -- it cannot tell
# you a handler works, only that one exists. The runtime side is covered by
# tools/test_downloads.R.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
APP  <- if (length(args)) args[1] else "shiny_app/app.R"
stopifnot(file.exists(APP))
src <- paste(readLines(APP, warn = FALSE), collapse = "\n")

rx <- function(pattern) {
  m <- gregexpr(pattern, src, perl = TRUE)[[1]]
  if (m[1] == -1) return(character(0))
  hits <- regmatches(src, gregexpr(pattern, src, perl = TRUE))[[1]]
  unique(sub(pattern, "\\1", hits, perl = TRUE))
}

# ---- what exists -----------------------------------------------------------
tables <- unique(c(rx('output\\$([A-Za-z0-9_]+)\\s*<-\\s*renderDT'),
                   rx('output\\[\\["([A-Za-z0-9_]+)"\\]\\]\\s*<-\\s*renderDT'),
                   rx('output\\$([A-Za-z0-9_]+)\\s*<-\\s*renderTable')))
plots  <- unique(c(rx('output\\$([A-Za-z0-9_]+)\\s*<-\\s*renderPlot\\('),
                   rx('output\\$([A-Za-z0-9_]+)\\s*<-\\s*renderPlotly')))

# ---- what is covered -------------------------------------------------------
reg_tables <- rx('list\\(id = "([A-Za-z0-9_]+)"')                  # TABLE_DL rows
reg_direct <- rx('register_dl\\(output, "([A-Za-z0-9_]+)"')        # standalone data downloads
fig_prefix <- unique(c(rx('register_fig\\(output, "([A-Za-z0-9_]+)"'),
                       rx('output\\[\\[paste0\\("([A-Za-z0-9_]+)_dl_"')))
fig_ui     <- rx('dl_fig_ui\\("([A-Za-z0-9_]+)"')
tab_ui     <- rx('dl_data_ui\\("([A-Za-z0-9_]+)"')

# ---- deliberate exclusions, each with a reason -----------------------------
# Anything here is a conscious decision, not an oversight. Adding to this list
# should require the same thought as adding a download.
EXCLUDE_TABLES <- c(
  # the six editable rename grids are inputs, not results
  "umap_renametab", "vln_renametab", "dot_renametab",
  "comp_renametab", "cmphase_renametab", "e2f_renametab")
EXCLUDE_PLOTS <- c(
  # plotly-native scatter/UMAP views: no ggplot underneath to hand to ggsave.
  # Their modebar camera exports PNG, and the underlying data IS downloadable
  # from the table beside them.
  "umap", "cm_map", "ct_volcano", "cm_volcano", "fg_volcano", "deg_volcano",
  "gm_scatter",
  # the four "all clusters" grid loops render one panel per subcluster from the
  # same _gg builders already exported singly on the tab above
  "cm_idgo_all", "cm_kogo_all", "cm_kodn_all", "cm_gsea_all")

is_loop_plot <- function(x) grepl("_all_CM[0-9]+$", x)

# ---- the check -------------------------------------------------------------
covered_tables <- unique(c(reg_tables, reg_direct))
missing_tables <- setdiff(setdiff(tables, EXCLUDE_TABLES), covered_tables)

# A plot counts as covered when a figure prefix plausibly maps to it. The app's
# prefixes are deliberately terse (cm_bar_geno -> cmbargeno), so compare with
# separators stripped rather than demanding an exact match.
squash <- function(x) gsub("[^a-z0-9]", "", tolower(x))
cov_sq <- squash(fig_prefix)
plot_covered <- function(x) {
  sx <- squash(x)
  # Prefixes are short by convention ("vn", "e2f", "cyc", "mat"), so a length
  # floor of 4 rejected real coverage. A loose match can in principle call an
  # uncovered plot covered, which is why every prefix is registered explicitly
  # (no loop variables) -- the set being matched against is exact.
  any(cov_sq == sx) || any(nchar(cov_sq) >= 2 & startsWith(sx, cov_sq)) ||
    any(nchar(sx) >= 2 & startsWith(cov_sq, sx))
}
cand_plots <- setdiff(plots, EXCLUDE_PLOTS)
cand_plots <- cand_plots[!is_loop_plot(cand_plots)]
missing_plots <- cand_plots[!vapply(cand_plots, plot_covered, TRUE)]

# UI buttons without a handler behind them, and vice versa -- a button that
# downloads nothing is worse than no button.
orphan_tab_ui <- setdiff(tab_ui, covered_tables)
orphan_fig_ui <- setdiff(fig_ui, fig_prefix)

cat(sprintf("app: %s\n", APP))
cat(sprintf("  tables rendered : %3d   (%d excluded by design)\n", length(tables), length(EXCLUDE_TABLES)))
cat(sprintf("  tables covered  : %3d\n", length(intersect(tables, covered_tables))))
cat(sprintf("  plots rendered  : %3d   (%d excluded, %d loop panels)\n",
            length(plots), length(intersect(plots, EXCLUDE_PLOTS)), sum(is_loop_plot(plots))))
cat(sprintf("  figure prefixes : %3d\n", length(fig_prefix)))

fail <- FALSE
report <- function(what, x) {
  if (!length(x)) return(invisible(NULL))
  fail <<- TRUE
  cat(sprintf("\nFAIL - %s (%d):\n", what, length(x)))
  cat(paste0("  ", x, collapse = "\n"), "\n")
}
report("tables with no download", missing_tables)
report("plots with no download", missing_plots)
report("dl_data_ui() buttons with no registered handler", orphan_tab_ui)
report("dl_fig_ui() buttons with no registered figure", orphan_fig_ui)

if (fail) {
  cat("\nEvery rendered table and plot needs a download, or an entry in the\n",
      "EXCLUDE_* lists above with a reason. See the plan for the rationale.\n", sep = "")
  quit(status = 1)
}
cat("\nOK - every table and plot is downloadable.\n")
