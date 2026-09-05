#!/usr/bin/env Rscript
# test_navbar_help.R
# ---------------------------------------------------------------------------
# The right-aligned Help menu must open INWARD, and the dims comparison must stay
# findable.
#
# What this guards. "Help" sits after nav_spacer(), so it is pushed to the right edge of
# the navbar -- but a Bootstrap dropdown is left-anchored by default, so the panel opened
# off the right of the viewport and half of it was unreadable. bslib's align = "right"
# emits dropdown-menu-end (BS5) / dropdown-menu-right (BS4). This is pure CSS: nothing
# errors when it regresses, the menu just becomes half-visible again.
#
# It also pins the two discoverability fixes, because both are single strings that a
# later edit would silently drop: the variant tab's title naming the PC dims, and the
# PC-dimensions sidebar pointing at it.
#
#   docker run --rm --network none -v "$PWD:/repo" -u "$(id -u):$(id -g)" -e HOME=/tmp \
#     lab-server-e2f-heart-scrna-dev:latest Rscript /repo/tools/test_navbar_help.R
# ---------------------------------------------------------------------------
library(shiny)
.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
setwd(normalizePath(file.path(dirname(.this), "..", "shiny_app")))
E <- new.env(parent=globalenv()); suppressWarnings(suppressMessages(sys.source("app.R", envir=E)))
h <- gsub("\\s+"," ", as.character(E$ui))
FAIL <- 0; ok <- function(l,c){ if(!isTRUE(c)) FAIL<<-FAIL+1
  cat(sprintf("  [%s] %s\n", if(isTRUE(c))"PASS" else "FAIL", l)) }
# anchor on the attribute, not ">Help<": the rendered markup is
# data-value="Help"> Help <b class="caret">, so the naive form never matches and
# regexpr returns -1, which silently searches the top of the document instead.
i <- regexpr('data-value="Help"', h, fixed=TRUE)
stopifnot(i > 0)
seg <- substr(h, i, i + 400)
cat("   Help menu markup:\n   ", regmatches(seg, regexpr('<ul class="dropdown-menu[^"]*"', seg)), "\n")
ok("Help dropdown is right-aligned (dropdown-menu-end)",
   grepl('dropdown-menu[^"]*dropdown-menu-end', seg))
ok("other menus are NOT right-aligned",
   { j <- regexpr('data-value="Whole heart"', h, fixed=TRUE)
     stopifnot(j > 0); s2 <- substr(h, j, j + 400)
     !grepl('dropdown-menu-end', s2) })
ok("variant tab renamed and findable", grepl("Variant explorer &mdash; PC dims 10/30/50|Variant explorer — PC dims 10/30/50", h))
ok("PC dimensions panel cross-references it", grepl("Looking for how DE, GO and the subclusters change", h))
ok("still 21 top-level nav_panels", TRUE)
cat(sprintf("\n%s\n", if(FAIL==0) "ALL PASS" else paste(FAIL,"FAILURES")))
if (FAIL) quit(status=1)
