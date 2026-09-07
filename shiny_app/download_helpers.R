# download_helpers.R
# ---------------------------------------------------------------------------
# Shared data-export helpers. Sourced by BOTH:
#   * shiny_app/app.R                      (every table download in the browser)
#   * analysis/<date>_<name>/03_excel.R    (the offline collaborator workbooks)
#
# They live in one file so the workbook a collaborator gets by email and the one
# they get by clicking Download are the same artefact. Everything here is pure
# data -> file; nothing depends on shiny, so the offline scripts can source it
# without pulling the app in.
#
# The workbook conventions (frozen header, autofilter, 3-significant-figure
# numerics, a README sheet carrying the caveats, the 31-character sheet-name
# guard) were written and tested for the 2026-08-21 deliverable; this is that
# code lifted out rather than reimplemented.
# ---------------------------------------------------------------------------

# ---- caveats that travel with every export --------------------------------
# A table handed over as a file loses whatever the surrounding page said about
# it, so the caveats ride along inside the workbook instead of on the tab.
DL_CAVEATS <- c(
  "n = 1 animal per genotype x timepoint. The two lanes per sample are the same library sequenced twice (technical, not biological, replicates). No contrast here has biological replication: cell-level Wilcoxon p-values are pseudoreplicated. Treat them as a RANKING, not as a hypothesis test.",
  "The KO and WT animals are different sexes. Y-linked genes (Eif2s3y, Kdm5d, Uty, Ddx3y), Xist and Tsix therefore top every KO-vs-WT list. They are flagged in the 'confounder' column and are excluded from every GO and GSEA input.",
  "E2f7 and E2f8 mRNA are NOT reduced in the KO in this dataset - most likely a conditional allele that a 3'-biased assay cannot see. Do not read the E2f7/E2f8 rows as a knockdown check.",
  "Mitochondrially-encoded (mt-) genes are up in KO in every cardiomyocyte subcluster and down in none, which is a library read-fraction difference rather than biology. They are flagged in the 'mito_encoded' column where present. It matters mainly for KO-up enrichment: those genes carry the oxidative-phosphorylation and electron-transport terms.",
  "P0 and P7 were not sorted identically. Within cardiomyocytes the cycling fraction is 16.3% at P0-WT and 25.0% at P7-WT (a ratio of 1.53x), and phase-matched vs raw log2 fold changes correlate at r = 0.99, so the effect on any P0-vs-P7 comparison is small - but the G1-matched tables are the safer read."
)

# ---- small utilities -------------------------------------------------------

# Excel caps sheet names at 31 chars and forbids []:*?/\ ; collisions are worse
# than truncation because openxlsx errors out halfway through a workbook.
dl_sheet_name <- function(x, taken = character()) {
  x <- gsub("[\\[\\]:*?/\\\\]", "-", x)
  x <- substr(x, 1, 31)
  if (!x %in% taken) return(x)
  for (i in 2:99) {
    cand <- paste0(substr(x, 1, 31 - nchar(as.character(i)) - 1), "_", i)
    if (!cand %in% taken) return(cand)
  }
  substr(paste0(x, "_", as.integer(runif(1, 100, 999))), 1, 31)
}

# 3 significant figures everywhere except counts, which must stay exact.
dl_signif <- function(df, keep = c("n_input", "n_universe", "Count", "size", "n_cells")) {
  for (nm in names(df)) {
    if (is.numeric(df[[nm]]) && !nm %in% keep && !grepl("^n_", nm))
      df[[nm]] <- signif(df[[nm]], 3)
  }
  df
}

# Column widths: explicit, not "auto". openxlsx's auto-fit scans every cell to
# size the column, which across sheets of ~20,000 rows is minutes of work for a
# cosmetic result -- and it makes long gene-list columns unusably wide anyway.
dl_widths <- function(df) {
  w <- rep(13, max(1, ncol(df)))
  wide <- c("gene", "Description", "geneID", "pathway", "leadingEdge", "sig_sets",
            "input_rule", "top_up_by_padj", "top_down_by_padj", "example_lost")
  w[names(df) %in% wide] <- 32
  w[names(df) == "Description"] <- 55
  w[names(df) %in% c("geneID", "leadingEdge")] <- 60
  w[names(df) %in% c("confounder", "log2FoldChange", "mito_encoded", "direction")] <- 15
  w
}

dl_stamp <- function() format(Sys.Date(), "%Y-%m-%d")

# ---- CSV -------------------------------------------------------------------
dl_write_csv <- function(df, file) utils::write.csv(df, file, row.names = FALSE, na = "")

# ---- workbook --------------------------------------------------------------
# sheets: a NAMED list of data frames, written in order. A character vector
# instead of a data frame is written as a notes sheet (one line per element),
# which is how the README sheet is produced.
dl_write_xlsx <- function(sheets, file, title = NULL, notes = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("openxlsx is required to write .xlsx - add it to the image.")
  stopifnot(is.list(sheets), length(sheets) > 0, !is.null(names(sheets)))

  hdr  <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#EEEEEE",
                                halign = "left", border = "bottom", borderColour = "#999999")
  wrap <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  wb   <- openxlsx::createWorkbook()

  # README first, so it is what opens.
  readme <- c(
    if (!is.null(title)) c(title, "") else NULL,
    sprintf("Exported %s from the E2F7/8 heart scRNA-seq browser.", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    if (!is.null(notes)) c(notes, "") else NULL,
    "PLEASE READ BEFORE INTERPRETING",
    paste0("  ", seq_along(DL_CAVEATS), ". ", DL_CAVEATS),
    "",
    "SHEETS IN THIS WORKBOOK",
    paste0("  ", names(sheets), "  (",
           vapply(sheets, function(x) if (is.data.frame(x)) sprintf("%d rows", nrow(x)) else "notes", ""), ")")
  )
  openxlsx::addWorksheet(wb, "README")
  openxlsx::writeData(wb, "README", data.frame(x = readme), colNames = FALSE)
  openxlsx::setColWidths(wb, "README", cols = 1, widths = 120)
  openxlsx::addStyle(wb, "README", wrap, rows = seq_along(readme), cols = 1, gridExpand = TRUE)

  taken <- "README"
  for (nm in names(sheets)) {
    d  <- sheets[[nm]]
    sn <- dl_sheet_name(nm, taken); taken <- c(taken, sn)
    openxlsx::addWorksheet(wb, sn)
    if (is.data.frame(d)) {
      if (!nrow(d)) d <- data.frame(note = "No rows for this selection.")
      d <- dl_signif(d)
      openxlsx::writeData(wb, sn, d, headerStyle = hdr, withFilter = TRUE)
      openxlsx::freezePane(wb, sn, firstActiveRow = 2)
      openxlsx::setColWidths(wb, sn, cols = seq_len(max(1, ncol(d))), widths = dl_widths(d))
    } else {
      openxlsx::writeData(wb, sn, data.frame(x = as.character(d)), colNames = FALSE)
      openxlsx::setColWidths(wb, sn, cols = 1, widths = 120)
      openxlsx::addStyle(wb, sn, wrap, rows = seq_along(d), cols = 1, gridExpand = TRUE)
    }
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(file)
}
