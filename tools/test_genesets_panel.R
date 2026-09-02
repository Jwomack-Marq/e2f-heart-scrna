#!/usr/bin/env Rscript
# test_genesets_panel.R
# ---------------------------------------------------------------------------
# The "Gene sets & sources" panel: layout, tab order, and the drift badge.
#
# What this guards. Two things that fail silently:
#
#   1. LAYOUT. The card holding the tabset must wrap its children in
#      card_body(fillable = FALSE). Without it bslib flexes the card body and the
#      tabset -- whose tables ask for 420px of scroll -- gets crushed to a few
#      pixels. Nothing errors; the page just renders unusably, which is how this
#      reached a user twice already.
#
#   2. THE DRIFT BADGE. The duplicate-name check is a real defect report (two panels
#      currently share a name but hold different genes). Moving it off the top of the
#      page into a tab means it is only discoverable via the count badge in the tab
#      title. If the badge silently stops rendering, an active error becomes invisible.
#
# Asserts against the BUILT UI html, not the source text, so a refactor that keeps the
# code looking right but changes what bslib emits still fails.
#
#   docker run --rm --network none -v "$PWD:/repo" -u "$(id -u):$(id -g)" -e HOME=/tmp \
#     lab-server-e2f-heart-scrna-dev:latest Rscript /repo/tools/test_genesets_panel.R
# ---------------------------------------------------------------------------
library(shiny)
.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
setwd(normalizePath(file.path(dirname(.this), "..", "shiny_app")))
E <- new.env(parent = globalenv()); suppressWarnings(suppressMessages(sys.source("app.R", envir = E)))
FAIL <- 0
ok <- function(l, c) { if(!isTRUE(c)) FAIL <<- FAIL+1; cat(sprintf("  [%s] %s\n", if(isTRUE(c))"PASS" else "FAIL", l)) }
h  <- gsub("\\s+", " ", as.character(E$ui))   # collapse ALL whitespace, indentation included
grab <- function(from, n) { i <- regexpr(from, h, fixed = TRUE); substr(h, i, i + n) }

cat("== the badge ==\n")
t1 <- as.character(E$gsp_checks_title())
ok("names the tab",              grepl("Summary &amp; checks", t1))
ok("badge carries the count 2",  grepl('badge bg-danger', t1) && grepl(">2<", t1))
ok("no badge when nothing drifted",
   !grepl("badge", as.character(E$gsp_checks_title(list(drift = data.frame(drifted = c(FALSE, FALSE)))))))
ok("degrades cleanly with no drift table", identical(E$gsp_checks_title(list()), "Summary & checks"))

cat("\n== tab order and default ==\n")
seg  <- grab('id="gsp_tabs"', 1800)
nav  <- substr(seg, 1, regexpr("</ul>", seg, fixed = TRUE))   # the tab strip only, not the panes
vals <- regmatches(nav, gregexpr('data-value="[a-z]+"', nav))[[1]]
vals <- sub('.*"([a-z]+)"', '\\1', vals)
cat(sprintf("    %s\n", paste(vals, collapse = " | ")))
ok("four tabs, checks first",  identical(vals, c("checks","reg","bench","refs")))
ok("Registry is the default",  grepl('<li class="active"> <a [^>]*data-value="reg"', seg))
ok("checks tab is NOT default", !grepl('<li class="active"> <a [^>]*data-value="checks"', seg))

cat("\n== the crush fix: the body holding the tabset must not flex its children ==\n")
# the card-body that directly contains the tabset
cb <- grab('<div class="card-header bslib-gap-spacing">Gene sets and their provenance</div>', 260)
ok("that card-body is a fill ITEM but not a fill CONTAINER",
   grepl('class="card-body html-fill-item"', cb) && !grepl('card-body[^"]*html-fill-container', cb))
ok("matches the known-good PC-dimensions panel",
   { pcd <- grab('id="pcd_label"', 400)
     grepl('class="card-body html-fill-item"', pcd) || grepl('card-body html-fill-item', pcd) })

cat("\n== banners moved, not deleted; tables still there ==\n")
pre <- grab('>Gene sets and their provenance</div>', 200)
ok("banners no longer above the tabset", !grepl('gsp_headline', pre))
ok("banners now inside the checks pane",
   grepl('data-value="checks"[^>]*>.{0,300}gsp_headline', seg, perl = TRUE) ||
   grepl('gsp_headline', grab('data-value="checks" id=', 600)))
ok("registry table present",  grepl('id="gsp_tab"', h))
ok("benchmark table present", grepl('id="gsp_bench"', h))
ok("caveats still under Registry", grepl('gsp_caveats', h))

cat(sprintf("\n%s\n", if (FAIL==0) "ALL PASS" else paste(FAIL,"FAILURES")))
if (FAIL) quit(status=1)
