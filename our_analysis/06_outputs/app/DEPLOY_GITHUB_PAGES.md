# Deploying the interactive cell browser

Two ways to put the cell browser online. **Option A (recommended)** is what you
already do for `contrasena-anki`: a free static site on GitHub Pages, no server.

---

## What gets built (run these once, in order)

From R (or RStudio), with the working directory at the project:

```r
RS <- "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"   # your R
# 1. Build the slim data object the app reads (a few minutes; loads the big .rds)
system2(RS, "our_analysis/06_outputs/app/build_app_data.R")
# 2. Preview the app locally (optional sanity check)
shiny::runApp("our_analysis/06_outputs/app")
# 3. Export to a static, in-browser site
system2(RS, "our_analysis/06_outputs/app/export_shinylive.R")
```

This produces `our_analysis/06_outputs/app/site/` — a self-contained static website
(HTML + JavaScript + WebAssembly + the data). It runs entirely in the visitor's
browser; there is no server and nothing for your advisor to install.

Preview the exported site locally before deploying:

```r
httpuv::runStaticServer("our_analysis/06_outputs/app/site", port = 8008)
# then open http://localhost:8008
```

---

## Option A — GitHub Pages (free, static, recommended)

Exactly the pattern you used for `contrasena-anki`.

1. Create a new GitHub repo, e.g. `e2f-heart-scrna`.
2. Copy **the contents of** `app/site/` into the repo root (so `index.html` is at
   the top level of the repo).
3. Add the workflow: copy `app/github-pages-static.yml` to the repo at
   `.github/workflows/static.yml`.
4. Commit and push to the `main` branch.
5. The Actions workflow runs automatically and enables Pages. Your live URL appears
   under the repo's **Settings → Pages** (and in the Actions run summary):
   `https://<your-username>.github.io/e2f-heart-scrna/`
6. Put that URL in **`our_analysis/06_outputs/app/LIVE_URL.txt`** (one line), then
   re-render the report so its "Explore the data yourself" link goes live:
   ```r
   Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
   system2(RS, "our_analysis/06_outputs/interactive/render_report.R")
   ```

> Note: the exported `site/` can be tens of MB (it bundles a small R runtime + the
> data). That is fine for GitHub Pages. First load takes a few seconds while the
> WebAssembly runtime initialises; afterwards it is snappy and cached.

---

## Option B — shinyapps.io (server, fallback)

Use this only if you prefer a server app (no gene-panel size limit).

1. Free account at <https://www.shinyapps.io>.
2. In the dashboard: **Account → Tokens → Show**, copy the `setAccountInfo(...)`
   line and run it once in R.
3. Deploy:
   ```r
   system2(RS, "our_analysis/06_outputs/app/deploy_shinyapps.R")
   ```
4. Copy the printed `https://<acct>.shinyapps.io/e2f-heart-scrna/` URL into
   `LIVE_URL.txt` and re-render the report (as in Option A step 6).

Free tier: 5 apps and ~25 active hours/month — plenty for advisor review.

---

## Files in this folder

| File | Purpose |
|------|---------|
| `build_app_data.R` | Extracts the slim `app_data.rds` from the Seurat objects |
| `app.R` | The Shiny cell browser (ggplot + Matrix; webR-compatible) |
| `export_shinylive.R` | Builds the static `site/` for GitHub Pages |
| `github-pages-static.yml` | Copy to `.github/workflows/static.yml` in the Pages repo |
| `deploy_shinyapps.R` | Fallback server deploy to shinyapps.io |
| `LIVE_URL.txt` | Paste the deployed URL here; the report links to it |
| `app_data.rds` | (generated) the data the app reads |
| `site/` | (generated) the static site to publish |
