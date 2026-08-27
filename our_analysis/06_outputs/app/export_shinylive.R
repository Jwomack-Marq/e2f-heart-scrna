#!/usr/bin/env Rscript
# Export the Shiny app to a STATIC shinylive site (HTML/JS/WebAssembly) that runs
# entirely in the browser — no server. The resulting `site/` folder is what gets
# pushed to GitHub Pages (same static-hosting pattern as the contrasena-anki app).
#
#   Rscript export_shinylive.R
#   -> our_analysis/06_outputs/app/site/   (open site/index.html, or push to Pages)
#
# Requires app_data.rds to exist (run build_app_data.R first).

this    <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
APP_DIR <- if (length(this) && nzchar(this)) dirname(normalizePath(this)) else getwd()
SITE    <- file.path(APP_DIR, "site")

if (!requireNamespace("shinylive", quietly = TRUE)) stop("install.packages('shinylive') first")
if (!file.exists(file.path(APP_DIR, "app_data.rds")))
  stop("app_data.rds not found — run build_app_data.R first.")

# Only the files the app needs go into the export (not the build/deploy scripts).
tmp <- file.path(tempdir(), "e2f_app"); unlink(tmp, recursive = TRUE); dir.create(tmp)
file.copy(file.path(APP_DIR, c("app.R", "app_data.rds")), tmp)

unlink(SITE, recursive = TRUE)
shinylive::export(appdir = tmp, destdir = SITE)

cat(sprintf("\nStatic site written to: %s\n", SITE))
cat("Preview locally:  Rscript -e \"httpuv::runStaticServer('", SITE, "', port=8008)\"\n", sep = "")
cat("Then push site/ to a GitHub Pages repo (see DEPLOY_GITHUB_PAGES.md).\n")
cat("=== DONE export_shinylive ===\n")
