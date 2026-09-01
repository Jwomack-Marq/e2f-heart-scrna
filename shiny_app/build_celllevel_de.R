# build_celllevel_de.R
# ---------------------------------------------------------------------------
# Carry the cell-level KO-vs-WT ranking into app_data.rds for the subclusters the
# pseudobulk pipeline cannot test.
#
# WHY. cm_subcluster_analyze.R aggregates cells into one pseudobulk sample per
# library x lane, requires 20 cells per sample and two samples per genotype, and skips a
# subcluster that cannot supply them. At res 0.2 that is CM12: 99 cells spread
# 10/10/18/21/8/8/10/14, so exactly one sample clears the floor. The subcluster then has
# no entry in app$tables$sub_DE at all -- which means it is not merely empty in the DE
# panel, it is ABSENT FROM THE DROPDOWN. A user clicking through subclusters never learns
# it exists, and the honest statement ("not testable the way the rest is tested") never
# gets made.
#
# This does not fix the testability. It carries the descriptive alternative
# (cm_subcluster_celllevel_de.R: Wilcoxon over cells, ranked by AUC) so the app can show
# the subcluster, say plainly why the usual result is missing, and offer a ranking that
# is clearly labelled as not a test.
#
#   Reads  results/tables/cm_subcluster_celllevel_res<R>_<CL>.csv
#   Writes app$cm_celllevel[["res<R>"]][["<CL>"]]
#
#   Rscript shiny_app/build_celllevel_de.R [--tables=<dir>] [--top=300]
# ---------------------------------------------------------------------------

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))
TABLES <- argval("--tables", "../our_analysis/results/tables")
TOP    <- as.integer(argval("--top", "300"))

fs <- list.files(TABLES, pattern = "^cm_subcluster_celllevel_res.*\\.csv$", full.names = TRUE)
if (!length(fs))
  stop("no cell-level tables in ", TABLES,
       "\n  run our_analysis/05_analyses/cm_subcluster_celllevel_de.R first")

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_celllevel.bak.rds\n")
saveRDS(app, "app_data.pre_celllevel.bak.rds", compress = "gzip")

out <- list()
for (f in fs) {
  key <- sub("^cm_subcluster_celllevel_(res[0-9.]+)_(CM[0-9]+)\\.csv$", "\\1|\\2", basename(f))
  res <- sub("\\|.*$", "", key); cl <- sub("^.*\\|", "", key)
  d <- read.csv(f, stringsAsFactors = FALSE)

  # The full table is every gene x every stratum (~30k rows, 5.6 MB for CM12 alone) and
  # nobody reads past the top of a ranking. Keep the top TOP per stratum by |AUC - 0.5|,
  # which is the order the script already writes in. Genes flagged as sex/construct
  # confounders stay in rather than being dropped here: the app has a "hide confounders"
  # toggle and removing them upstream would make that toggle silently do nothing.
  d <- do.call(rbind, lapply(split(d, d$stratum), function(x) {
    x <- x[order(-abs(x$auc - 0.5)), , drop = FALSE]
    head(x, TOP)
  }))
  d <- d[order(factor(d$stratum, levels = c("pooled", "P0", "P7")),
               -abs(d$auc - 0.5)), , drop = FALSE]
  rownames(d) <- NULL
  d$NOTE <- NULL          # the caveat belongs in the UI, not repeated on every row
  out[[res]][[cl]] <- d
  cat(sprintf("  %s %-6s  %d rows kept of %d  (strata: %s)\n", res, cl, nrow(d),
              nrow(read.csv(f, stringsAsFactors = FALSE)),
              paste(unique(d$stratum), collapse = ", ")))
}

app$cm_celllevel <- out

# Record WHY each of these exists, read off the pseudobulk summary rather than hardcoded,
# so the app can state the actual sample counts instead of a remembered sentence.
app$cm_celllevel_meta <- lapply(names(out), function(res) {
  s <- app$tables$sub_summary[[res]]
  cls <- names(out[[res]])
  if (is.null(s)) return(setNames(vector("list", length(cls)), cls))
  setNames(lapply(cls, function(cl) {
    r <- s[s$subcluster == cl, , drop = FALSE]
    if (!nrow(r)) return(NULL)
    list(n_cells = r$n_cells[1], n_KO_samp = r$n_KO_samp[1],
         n_WT_samp = r$n_WT_samp[1], status = r$status[1])
  }), cls)
})
names(app$cm_celllevel_meta) <- names(out)

cat("\n== Summary ==\n")
for (res in names(out)) for (cl in names(out[[res]])) {
  m <- app$cm_celllevel_meta[[res]][[cl]]
  cat(sprintf("  %s %s: %d cells, pseudobulk samples KO=%s WT=%s (%s)\n", res, cl,
              m$n_cells %||% NA, m$n_KO_samp %||% NA, m$n_WT_samp %||% NA, m$status %||% "?"))
}
cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat(sprintf("Done. Bundle is now %.1f MB.\n", file.size("app_data.rds") / 1024^2))
