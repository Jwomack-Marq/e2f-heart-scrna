#!/usr/bin/env Rscript
# test_celllevel_de.R
# ---------------------------------------------------------------------------
# The CM12 case: a subcluster the pseudobulk pipeline cannot test must still be
# reachable, and must say why.
#
# What this guards. Before build_celllevel_de.R, a subcluster with no sub_DE entry was
# not merely empty in the DE panel -- it was absent from the dropdown, so a reader
# clicking through subclusters never learned CM12 existed. The failure mode is silence,
# which is exactly what eyeballing the app does not catch: nothing looks wrong.
#
# It runs in two parts because testServer can see the server's reactives but NOT app.R's
# file-level objects, so the pure dropdown helper is tested against a sourced copy.
#
#   Run from the repo root, in the app image:
#     docker run --rm --network none -v "$PWD:/repo" -u "$(id -u):$(id -g)" -e HOME=/tmp \
#       lab-server-e2f-heart-scrna-dev:latest Rscript /repo/tools/test_celllevel_de.R
# ---------------------------------------------------------------------------
library(shiny)
# app.R reads app_data.rds and sources its helpers by relative path, so the working
# directory has to be shiny_app/. normalizePath because sys.source() will not accept
# the un-normalised ".." form.
.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
setwd(normalizePath(file.path(dirname(.this), "..", "shiny_app")))

FAIL <- 0
ok <- function(lbl, cond) { if (!isTRUE(cond)) FAIL <<- FAIL + 1
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", lbl)) }
# An output silenced by req() throws when read here, where a browser just leaves the slot
# empty. "blank" means either -- NULL, empty, or silenced.
blank <- function(expr) { v <- try(force(expr), silent = TRUE)
  inherits(v, "try-error") || is.null(v) || !nzchar(paste(unlist(v), collapse = "")) }

# --- part 1: the dropdown helper, as a pure function -----------------------
E <- new.env(parent = globalenv())
suppressWarnings(suppressMessages(sys.source("app.R", envir = E)))
cat("== dropdown helper ==\n")
ch <- E$cm_sub_choices(names(E$subDE[["res0.2"]]), names(E$CELLDE))
ok("CM12 is now offered",            "CM12" %in% unname(ch))
ok("CM12 marked not testable",       grepl("not testable", names(ch)[ch == "CM12"]))
ok("tested subclusters unmarked",    !any(grepl("not testable", names(ch)[ch != "CM12"])))
ok("numeric order kept (CM0..CM12)", identical(unname(ch), paste0("CM", 0:12)))
cat(sprintf("    CM12 label: %s\n", names(ch)[ch == "CM12"]))
ok("degrades to the old list on a bundle without the new step",
   identical(unname(E$cm_sub_choices(paste0("CM", 0:11), character(0))), paste0("CM", 0:11)))

cat("\n== the stored ranking ==\n")
d <- E$CELLDE[["CM12"]]
ok("900 rows, 3 strata", nrow(d) == 900 && setequal(unique(d$stratum), c("pooled","P0","P7")))
ok("auc kept alongside p-values", all(c("auc","pval","padj") %in% names(d)))
ok("ranked by |auc-0.5| within stratum",
   { p <- d[d$stratum == "pooled", ]; !is.unsorted(rev(abs(p$auc - 0.5))) })
ok("immune markers are flagged", sum(d$immune_gene) > 0)

# --- part 2: the reactives, driven ----------------------------------------
testServer(shinyAppFile("app.R"), {
  cat("\n== a normal subcluster (CM0) ==\n")
  session$setInputs(cm_sub = "CM0", cm_contrast = "pooled", cm_hideconf = FALSE)
  ok("cell-level mode OFF",         isFALSE(cm_celllevel_on()))
  ok("pseudobulk table has rows",   nrow(cm_d()) > 100)
  ok("pseudobulk table renders",    nchar(output$cm_detab) > 100)
  ok("cell-level note stays blank", blank(output$cm_celllevel_note))
  ok("cell-level table stays blank", blank(output$cm_celllevel_tab))

  cat("\n== the skipped subcluster (CM12) ==\n")
  session$setInputs(cm_sub = "CM12")
  ok("cell-level mode ON",          isTRUE(cm_celllevel_on()))
  ok("cell-level table renders",    nchar(output$cm_celllevel_tab) > 100)
  ok("pseudobulk note steps aside", blank(output$cm_de_note))
  ok("pseudobulk table stays blank", blank(output$cm_detab))
  h <- output$cm_celllevel_note$html
  ok("explanatory note renders",    nchar(h) > 500)
  ok("says RANKING not a test",     grepl("RANKING, not a test", h))
  ok("gives real cell counts",      grepl("36 KO and 63 WT cells", h))
  ok("names the real status",       grepl("skipped_too_few_or_unbalanced", h))
  ok("warns immune, not CM, biology", grepl("must not be read as cardiomyocyte biology", h))
  ok("ranking reactive returns rows", nrow(cm_celllevel_d()) == 900)

  cat("\n== hide-confounders works on the ranking too ==\n")
  session$setInputs(cm_hideconf = TRUE)
  n <- nrow(cm_celllevel_d())
  ok("confounders dropped", n < 900 && n > 0)
  cat(sprintf("    %d of 900 rows kept\n", n))

  cat("\n== back to CM0: no leakage ==\n")
  session$setInputs(cm_sub = "CM0")
  ok("cell-level mode OFF again", isFALSE(cm_celllevel_on()))
  ok("volcano data back",         nrow(cm_d()) > 100)
})
cat(sprintf("\n%s\n", if (FAIL == 0) "ALL PASS" else paste(FAIL, "FAILURES")))
if (FAIL) quit(status = 1)
