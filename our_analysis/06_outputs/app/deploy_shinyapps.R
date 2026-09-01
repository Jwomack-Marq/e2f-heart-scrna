#!/usr/bin/env Rscript
# FALLBACK deploy: host the cell browser as a server-side Shiny app on shinyapps.io.
# Use this only if the static shinylive/GitHub-Pages route (export_shinylive.R) is
# not desired. A server app has no payload limit, so it can use a larger gene panel.
#
# ONE-TIME setup (you, not your advisor):
#   1. Create a free account at https://www.shinyapps.io
#   2. Account -> Tokens -> Show -> copy the setAccountInfo() call, run it once in R:
#        rsconnect::setAccountInfo(name="<acct>", token="<token>", secret="<secret>")
#   3. Run this script. It prints the public URL when done.

this    <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
APP_DIR <- if (length(this) && nzchar(this)) dirname(normalizePath(this)) else getwd()

if (!requireNamespace("rsconnect", quietly = TRUE)) stop("install.packages('rsconnect') first")
if (!file.exists(file.path(APP_DIR, "app_data.rds")))
  stop("app_data.rds not found — run build_app_data.R first.")
if (!nrow(rsconnect::accounts()))
  stop("No shinyapps.io account configured. Run rsconnect::setAccountInfo(...) first (see header).")

rsconnect::deployApp(
  appDir   = APP_DIR,
  appFiles = c("app.R", "app_data.rds"),       # ship only what the app needs
  appName  = "e2f-heart-scrna",
  appTitle = "E2F7/8 heart scRNA-seq browser",
  forceUpdate = TRUE, launch.browser = FALSE
)
cat("=== DONE deploy_shinyapps — copy the printed https://*.shinyapps.io/ URL into",
    "app/LIVE_URL.txt and re-render the report ===\n")
