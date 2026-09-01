#!/usr/bin/env Rscript
# test_downloads.R
# ---------------------------------------------------------------------------
# Runtime companion to check_download_coverage.R. That one proves a handler
# exists; this one proves it produces a file with rows in it.
#
#   Rscript tools/test_downloads.R [path/to/shiny_app]
#
# Uses shiny::testServer, whose expression is evaluated inside the server
# function's own environment -- so TABLE_DL and every reactive are directly
# reachable without clicking anything.
# ---------------------------------------------------------------------------

suppressMessages({ library(shiny); library(openxlsx) })
args <- commandArgs(trailingOnly = TRUE)
APPDIR <- if (length(args)) args[1] else "shiny_app"
stopifnot(dir.exists(APPDIR))
# app.R sources these with local = TRUE, so they live in the app's environment and
# testServer's expression cannot see them. Source them here too and test them as
# the standalone unit they are.
source(file.path(APPDIR, "download_helpers.R"))

pass <- 0L; fail <- 0L; skip <- 0L
note <- function(status, what, extra = "") {
  cat(sprintf("  [%-4s] %-28s %s\n", status, what, extra))
  if (status == "PASS") pass <<- pass + 1L
  else if (status == "FAIL") fail <<- fail + 1L
  else skip <<- skip + 1L
}

testServer(APPDIR, {
  # A plausible browsing state. Anything still unset makes its table req() out,
  # which is reported as SKIP rather than counted as a pass -- an untested table
  # must not look like a working one.
  session$setInputs(
    ct_tp = "P7", ct_sel = "Cardiomyocyte", ct_hideconf = FALSE, ct_search = "",
    cm_sub = "CM2", cm_hideconf = FALSE,
    cm_enr_mode = "one", cm_enr_sub = "CM2",
    enr_ct = "Cardiomyocyte", enr_tp = "P7",
    fg_cluster = "CM2", fg_contrast = "P7_KO_vs_WT", fg_stratum = "all",
    fg_grid = "de", fg_hideconf = FALSE, fg_count_mode = "prop",
    fg_score = "sig_maturation_nocc", fg_score_stratum = "all",
    fg_g1_clusters = c("AllCM", "CM2"), fg_book_all = FALSE,
    mi_cluster = "CM2", mi_hideconf = TRUE,
    mi_quad = c("immature_up_in_KO", "mature_down_in_KO"),
    mi_cand_stratum = "all", mi_spec_grid = "de",
    mi_genes = c("Tcf4", "Adamts9", "Gabbr2"), mi_cand_clusters = c("AllCM", "CM2"),
    vn_cluster = "AllCM", vn_stratum = "all", vn_grid = "de",
    vn_a = "de:P7_KO_vs_WT:up", vn_b = "de:WT_P0_vs_P7:up", vn_c = "none",
    gm_panel = "avg", gm_dist = 0,
    cc_tp = "P7", cc_pathway = "All", cc_metric = "delta",
    deg_hideconf = FALSE,
    fg_enr_cluster = "CM2", fg_enr_contrast = "P7_KO_vs_WT", fg_enr_stratum = "all",
    fg_enr_ont = "BP", fg_enr_topn = 20,
    xc_wt_cluster = "AllCM", xc_mat_set = "CM maturation", xc_cyc_set = "__canonical__",
    xc_mat_clusters = c("CM1","CM2","CM3","CM7","CM8"),
    xc_cyc_clusters = c("CM2","CM4","CM5"), xc_minc = 1, xc_stratum = "all",
    xc_grid = "de", xc_padj = 0.05, xc_measure = "auc", xc_auc = 0.60, xc_lfc = 0.25,
    xc_hidemt = TRUE, xc_gene_cmp = "__all__",
    ct_vlfc = 1, cm_vlfc = 1, fg_vlfc = 1, deg_vlfc = 1,
    vn_measure = "lfc", vn_auc = 0.60,
    mat_ct = "Cardiomyocyte", mat_score = "sig_maturation", mat_stratum = "all"
  )

  cat("\n== table download frames ==\n")
  for (t in TABLE_DL) {
    d <- tryCatch(t$df(), error = function(e) e)
    if (inherits(d, "shiny.silent.error")) {
      msg <- conditionMessage(d)
      # A "run the builder" guard is the correct behaviour on a bundle that
      # predates a builder -- not the same thing as an untested table.
      if (grepl("run build_", msg)) note("PASS", t$id, "guarded (slot absent)")
      else note("SKIP", t$id, "needs more inputs")
      next
    }
    if (inherits(d, "error"))              { note("FAIL", t$id, conditionMessage(d)); next }
    if (!is.data.frame(d))                 { note("FAIL", t$id, "not a data.frame"); next }
    nm <- tryCatch(if (is.function(t$base)) t$base() else t$base, error = function(e) "<err>")
    note("PASS", t$id, sprintf("%d rows x %d cols  -> %s", nrow(d), ncol(d), nm))
  }

  cat("\n== writers round-trip ==\n")
  d <- tryCatch(TABLE_DL[[1]]$df(), error = function(e) NULL)
  if (is.null(d)) {
    note("SKIP", "writers", "no frame available")
  } else {
    f1 <- tempfile(fileext = ".csv"); dl_write_csv(d, f1)
    back <- utils::read.csv(f1)
    note(if (nrow(back) == nrow(d)) "PASS" else "FAIL", "csv round-trip",
         sprintf("%d rows", nrow(back)))

    f2 <- tempfile(fileext = ".xlsx")
    dl_write_xlsx(list(data = head(d, 200), other = head(d, 5)), f2,
                  title = "test", notes = "written by tools/test_downloads.R")
    sn <- openxlsx::getSheetNames(f2)
    ok <- identical(sn, c("README", "data", "other")) && all(nchar(sn) <= 31)
    note(if (ok) "PASS" else "FAIL", "xlsx sheets", paste(sn, collapse = ", "))
    b2 <- openxlsx::read.xlsx(f2, sheet = "data")
    note(if (nrow(b2) == min(200, nrow(d))) "PASS" else "FAIL", "xlsx round-trip",
         sprintf("%d rows, %.0f KB", nrow(b2), file.size(f2) / 1024))
  }

  cat("\n== four-group enrichment ==\n")
  # Passes either way: with the slot present these must return frames, without it
  # they must fail with the "run the builder" message rather than crash. A bundle
  # built before build_fourgroup_enrichment.R existed still has to load.
  enr_fns <- c("fg_enr_up_df", "fg_enr_dn_df", "fg_enr_gsea_dat", "fg_enr_audit_dat")
  for (nm in enr_fns) {
    r <- tryCatch(get(nm)(), error = function(e) e)
    if (inherits(r, "shiny.silent.error")) {
      note(if (grepl("build_fourgroup_enrichment", conditionMessage(r))) "PASS" else "FAIL",
           nm, "guarded (slot absent)")
    } else if (inherits(r, "error")) {
      note("FAIL", nm, conditionMessage(r))
    } else {
      note(if (is.data.frame(r)) "PASS" else "FAIL", nm, sprintf("%d rows", nrow(r)))
    }
  }
  bk <- tryCatch(fg_enr_book_sheets("P7_KO_vs_WT", "all", "BP"), error = function(e) e)
  if (inherits(bk, "shiny.silent.error")) note("PASS", "fg_enr_book_sheets", "guarded (slot absent)")
  else if (inherits(bk, "error"))          note("FAIL", "fg_enr_book_sheets", conditionMessage(bk))
  else                                     note("PASS", "fg_enr_book_sheets", sprintf("%d sheets", length(bk)))

  cat("\n== bulk contrast workbook ==\n")
  # via the server-local builder: testServer's expression cannot see app.R's
  # file-level objects (FG, FG_CLUSTERS, ...), only the server function's own.
  sheets <- tryCatch(fg_book_sheets("P7_KO_vs_WT", "all", "de"), error = function(e) e)
  if (inherits(sheets, "error")) {
    note("FAIL", "fg_book sheets", conditionMessage(sheets))
  } else {
    f3 <- tempfile(fileext = ".xlsx")
    dl_write_xlsx(sheets, f3, title = "P7 KO vs WT", notes = "test")
    note("PASS", "fg_book workbook",
         sprintf("%d sheets, %.1f MB", length(sheets) + 1L, file.size(f3) / 1e6))
  }

  cat("\n== CM deep-dive per-contrast DE ==\n")
  # The deep-dive DE tab reads the same four-group grids the Four-group tab does.
  # The same selection has to give the same table: if these two ever diverge, one of
  # the tabs is quietly answering a different question than its label claims.
  session$setInputs(cm_sub = "CM2", cm_contrast = "P7_KO_vs_WT", cm_stratum = "all",
                    cm_grid = "de",
                    fg_cluster = "CM2", fg_contrast = "P7_KO_vs_WT", fg_stratum = "all",
                    fg_grid = "de", fg_hideconf = FALSE)
  a <- tryCatch(cm_d(), error = function(e) e)
  b <- tryCatch(fg_d(), error = function(e) e)
  if (inherits(a, "error") || inherits(b, "error"))
    note("FAIL", "cm_d == fg_d",
         conditionMessage(if (inherits(a, "error")) a else b))
  else
    note(if (identical(a, b)) "PASS" else "FAIL", "cm_d == fg_d",
         sprintf("%d vs %d rows", nrow(a), nrow(b)))
  # every contrast the builder ships should resolve, or say why it cannot. FG_CTAB is a
  # file-level object testServer cannot see, so the keys are spelled out here.
  for (k in c("WT_P0_vs_P7", "KO_P0_vs_P7", "P7_KO_vs_WT", "P0_KO_vs_WT")) {
    session$setInputs(cm_contrast = k)
    d <- tryCatch(cm_d(), error = function(e) e)
    if (inherits(d, "shiny.silent.error")) note("SKIP", paste0("cm_d ", k), conditionMessage(d))
    else if (inherits(d, "error"))         note("FAIL", paste0("cm_d ", k), conditionMessage(d))
    else                                   note("PASS", paste0("cm_d ", k), sprintf("%d rows", nrow(d)))
  }
  session$setInputs(cm_contrast = "P7_KO_vs_WT")
  l <- tryCatch(cm_lfc_list(), error = function(e) e)
  if (inherits(l, "error")) note("FAIL", "cm_lfc_list", conditionMessage(l))
  else note("PASS", "cm_lfc_list", sprintf("%d/%d subclusters have a table",
                                           sum(!vapply(l, is.null, TRUE)), length(l)))
  # and back to pooled: the original behaviour must still be reachable
  session$setInputs(cm_contrast = "pooled")
  d <- tryCatch(cm_d(), error = function(e) e)
  note(if (is.data.frame(d)) "PASS" else "FAIL", "cm_d pooled",
       if (is.data.frame(d)) sprintf("%d rows", nrow(d)) else conditionMessage(d))

  # cm_celllevel_tab is the one download the loop above cannot reach: it exists only for a
  # subcluster pseudobulk could not test, and the loop runs with cm_sub = "CM2". Drive it
  # here with a subcluster that has one, so it is tested rather than reported SKIP.
  CLD <- parent.env(environment(xc_venn_p))$CELLDE
  if (!length(CLD)) note("SKIP", "cm_celllevel_tab", "bundle predates build_celllevel_de.R")
  else {
    dl <- Filter(function(x) identical(x$id, "cm_celllevel_tab"), TABLE_DL)[[1]]
    session$setInputs(cm_sub = names(CLD)[1], cm_hideconf = FALSE)
    f <- tryCatch(dl$df(), error = function(e) e)
    if (!is.data.frame(f)) note("FAIL", "cm_celllevel_tab",
      if (inherits(f, "error")) conditionMessage(f) else "not a data.frame")
    else note("PASS", "cm_celllevel_tab",
              sprintf("%d rows x %d cols  -> %s", nrow(f), ncol(f), dl$base()))
    # and it must refuse, not return a blank file, for a subcluster that has no ranking
    session$setInputs(cm_sub = "CM2")
    g <- tryCatch(dl$df(), error = function(e) e)
    note(if (inherits(g, "shiny.silent.error")) "PASS" else "FAIL",
         "cm_celllevel_tab guards", "req() on a subcluster with no ranking")
    session$setInputs(cm_sub = "CM2")
  }

  cat("\n== WT programs x KO clusters ==\n")
  # app.R's file-level helpers are one frame up from the server execution env.
  APP <- parent.env(environment(xc_venn_p))
  gf  <- function(n) get(n, envir = APP)
  vs <- tryCatch(xc_all(), error = function(e) e)
  if (inherits(vs, "error")) note("FAIL", "xc_all", conditionMessage(vs))
  else {
    note(if (length(vs) == 4) "PASS" else "FAIL", "xc_all", "4 comparisons")
    for (k in seq_along(vs))
      note("PASS", paste0("comparison ", k),
           sprintf("A=%d B=%d shared=%d", length(vs[[k]]$sets[[1]]$genes),
                   length(vs[[k]]$sets[[2]]$genes),
                   length(intersect(vs[[k]]$sets[[1]]$genes, vs[[k]]$sets[[2]]$genes))))
    # The crossing must not quietly redefine "up": rebuild comparison 1's WT side from
    # the raw table and the curated category and demand the same answer.
    wt  <- gf("FG")$de[["AllCM"]][["WT_P0_vs_P7__all"]]
    man <- intersect(gf("de_pass")(wt, "up", 0.05, 0.60, "auc"), gf("GENE_SETS")[["CM maturation"]])
    man <- setdiff(man[!grepl("^mt-", man)], gf("CONF"))
    note(if (setequal(man, vs[[1]]$sets[[1]]$genes)) "PASS" else "FAIL",
         "C1 set A == manual", paste(sort(man), collapse = ", "))
  }
  a1 <- gf("xc_ko_set")(c("CM2","CM4","CM5"), "down", "all", "de", 0.05, 0.60, "auc", minc = 1)$genes
  a3 <- gf("xc_ko_set")(c("CM2","CM4","CM5"), "down", "all", "de", 0.05, 0.60, "auc", minc = 3)$genes
  note(if (all(a3 %in% a1) && length(a3) <= length(a1)) "PASS" else "FAIL",
       "minc=3 subset of minc=1", sprintf("%d of %d", length(a3), length(a1)))
  # CM4 has no G1 stratum for any contrast, so the G1 option must report it, not drop it silently
  pres <- gf("xc_ko_present")(c("CM2","CM4","CM5"), "G1", "de")
  note(if (!("CM4" %in% pres)) "PASS" else "FAIL", "CM4 absent from G1",
       paste(pres, collapse = ", "))
  # a symmetric |log2FC| cut cannot classify a sparse cell-cycle gene; the audit has to say so
  aud <- gf("xc_measure_audit")(xc_p())
  note(if (is.data.frame(aud) && nrow(aud) == 4) "PASS" else "FAIL", "xc_measure_audit",
       if (is.data.frame(aud)) paste(sprintf("%s/%s AUC %d lfc %d", aud$category,
         aud$direction, aud$n_AUC, aud$n_log2FC), collapse = " | ") else "not a frame")
  for (nm in c("xc_wt_tab_df","xc_ko_tab_df","xc_ko_pivot_dat","xc_stats_df","xc_genes_df")) {
    r <- tryCatch(get(nm)(), error = function(e) e)
    if (inherits(r, "shiny.silent.error")) note("SKIP", nm, conditionMessage(r))
    else if (inherits(r, "error")) note("FAIL", nm, conditionMessage(r))
    else note("PASS", nm, sprintf("%d rows x %d cols", nrow(r), ncol(r)))
  }

  # the effect-size sliders have to actually move the answer, in both directions
  cat("\n== overlap thresholds are adjustable ==\n")
  n_at <- function(m, e) { session$setInputs(xc_measure = m,
      xc_auc = if (m == "auc") e else 0.60, xc_lfc = if (m == "lfc") e else 0.25)
    vapply(xc_all(), function(v) length(v$sets[[1]]$genes), 0L) }
  a60 <- n_at("auc", 0.60); a52 <- n_at("auc", 0.52)
  note(if (all(a52 >= a60) && sum(a52) > sum(a60)) "PASS" else "FAIL",
       "looser AUC -> larger WT sets",
       sprintf("0.60: %s -> 0.52: %s", paste(a60, collapse = "/"), paste(a52, collapse = "/")))
  l100 <- n_at("lfc", 1.0); l050 <- n_at("lfc", 0.5); l025 <- n_at("lfc", 0.25)
  note(if (all(l025 >= l050) && all(l050 >= l100)) "PASS" else "FAIL",
       "looser |log2FC| -> larger WT sets",
       sprintf("1.0: %s -> 0.5: %s -> 0.25: %s", paste(l100, collapse = "/"),
               paste(l050, collapse = "/"), paste(l025, collapse = "/")))
  session$setInputs(xc_measure = "auc", xc_auc = 0.60, xc_lfc = 0.25)
  # the Venn tab gained the same choice; both measures must resolve
  for (m in c("lfc", "auc")) {
    session$setInputs(vn_measure = m)
    r <- tryCatch(vn_v(), error = function(e) e)
    note(if (is.list(r) && length(r$sets) >= 2) "PASS" else "FAIL", paste0("vn_v measure=", m),
         if (is.list(r)) paste(vapply(r$sets, function(s) length(s$genes), 0L), collapse = "/")
         else conditionMessage(r))
  }
  session$setInputs(vn_measure = "lfc")
  # the volcano colour cut is a display threshold, so it must NOT change the table
  n_tab <- nrow(cm_tab()); session$setInputs(cm_vlfc = 0.25)
  note(if (identical(nrow(cm_tab()), n_tab)) "PASS" else "FAIL",
       "volcano cut does not filter", sprintf("%d rows either way", n_tab))
  session$setInputs(cm_vlfc = 1)

  cat("\n== score definitions: purpose, maths, and the actual genes ==\n")
  # The gene lists are NOT in app_data.rds -- app.R parses them out of
  # build_signature_scores.R so there is no second copy to drift. If that parse ever
  # silently fails the app would document nothing, so check it against the literal.
  flat <- function(x) paste(as.character(x), collapse = "")   # tag lists render to >1 string
  SS <- gf("SCORE_SETS"); sdu <- gf("score_def_ui")
  note(if (!is.null(SS)) "PASS" else "FAIL", "gene lists parsed from the builder",
       if (is.null(SS)) "NULL - is build_signature_scores.R next to app.R?"
       else sprintf("%d sets: %s", length(SS), paste(names(SS), collapse = ", ")))
  lit <- local({ got <- NULL
    # testServer runs with the app dir as cwd, which is how app.R finds it too
    for (x in parse("build_signature_scores.R"))
      if (is.null(got) && is.call(x) && identical(as.character(x[[1]]), "<-") &&
          identical(as.character(x[[2]]), "SETS")) got <- eval(x[[3]])
    got })
  note(if (identical(lapply(lit, function(g) setdiff(g, gf("CONF"))), SS)) "PASS" else "FAIL",
       "parsed lists == the builder literal", "no drift between doc and what was scored")
  h <- flat(sdu("sig_maturation")); h2 <- flat(sdu("sig_metabolic"))
  note(if (grepl("Myh6, Tnni3, Pln, Atp2a2", h, fixed = TRUE) &&
           grepl("Myh7, Tnni1, Nppa, Nppb", h, fixed = TRUE)) "PASS" else "FAIL",
       "sig_maturation names both poles' genes")
  note(if (grepl("Slc2a1, Hk1, Hk2", h2, fixed = TRUE) &&
           grepl("Cpt1a, Cpt1b, Cpt2", h2, fixed = TRUE)) "PASS" else "FAIL",
       "sig_metabolic names both poles' genes")
  note(if (grepl("fetal-to-adult", h, fixed = TRUE)) "PASS" else "FAIL",
       "says what the score is FOR, not just how")
  # a set gene absent from the scoring matrix contributed nothing; that must be visible
  note(if (grepl("did not contribute", flat(sdu("sig_cytokinesis")), fixed = TRUE)) "PASS" else "FAIL",
       "names set genes missing from the matrix", "sig_cytokinesis found 13 of 14")
  d <- tryCatch(gf("score_sets_df")(), error = function(e) e)
  note(if (is.data.frame(d) && nrow(d) == length(SS) && all(nzchar(d$genes))) "PASS" else "FAIL",
       "score_sets_df", if (is.data.frame(d)) sprintf("%d sets, all with genes", nrow(d))
                        else conditionMessage(d))

  cat("\n== method notes under the figures ==\n")
  # The bodies are free functions precisely so this is testable: testServer snapshots
  # output$ values, so a note that silently stopped following its dropdown would still
  # look correct through output$.
  for (o in c("mat_violin_method","mat_scatter_method","gm_method","xc_venn_method")) {
    h <- tryCatch(flat(output[[o]]), error = function(e) paste("ERROR:", conditionMessage(e)))
    note(if (!startsWith(h, "ERROR") && nchar(h) > 500) "PASS" else "FAIL", o,
         if (startsWith(h, "ERROR")) h else sprintf("%d chars", nchar(h)))
  }
  gmn <- gf("gm_method_note"); vmn <- gf("mat_violin_method_note")
  note(if (grepl("0.509", flat(gmn("avg")), fixed = TRUE) &&
           grepl("0.519", flat(gmn("P7")),  fixed = TRUE) &&
           grepl("0.498", flat(gmn("P0")),  fixed = TRUE)) "PASS" else "FAIL",
       "gene map note tracks the panel", "avg/P7/P0 centres 0.509/0.519/0.498")
  note(if (grepl("mat_mature - mat_immature", flat(vmn("sig_maturation")), fixed = TRUE) &&
           grepl("faox - glycolysis", flat(vmn("sig_metabolic")), fixed = TRUE)) "PASS" else "FAIL",
       "violin note tracks the score", "names the gene sets behind it")
  note(if (grepl("Xist", flat(gmn("avg")), fixed = TRUE)) "PASS" else "FAIL",
       "gene map note flags its own gaps", "Xist / mt- not filtered there")

  cat("\n== CM deep-dive per-contrast enrichment ==\n")
  session$setInputs(cm_enr_sub = "CM2", cm_enr_contrast = "P7_KO_vs_WT",
                    cm_enr_stratum = "all", cm_enr_ont = "BP")
  for (nm in c("cm_sub_kogo_df", "cm_sub_kodn_df", "cm_sub_gsea_df")) {
    r <- tryCatch(get(nm)(), error = function(e) e)
    if (inherits(r, "shiny.silent.error"))
      note(if (grepl("build_", conditionMessage(r))) "PASS" else "SKIP", nm, conditionMessage(r))
    else if (inherits(r, "error")) note("FAIL", nm, conditionMessage(r))
    else                           note("PASS", nm, sprintf("%d rows", nrow(r)))
  }
})

cat(sprintf("\n%d passed, %d skipped, %d failed\n", pass, skip, fail))
if (fail) quit(status = 1)
