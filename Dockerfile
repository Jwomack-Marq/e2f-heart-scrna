# syntax=docker/dockerfile:1
FROM rocker/shiny:4.4.1

# System libraries for the R packages the app renders with:
#   ggplot2/ragg -> png/jpeg/tiff/freetype/fontconfig ; svglite/textshaping -> harfbuzz/fribidi ;
#   plotly/DT/curl/openssl -> xml2/curl/ssl. Most exist in rocker/shiny; kept explicit for robustness.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libxml2-dev libcurl4-openssl-dev libssl-dev \
      libpng-dev libjpeg-dev libtiff5-dev \
      libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages the app loads (shiny/htmlwidgets/crosstalk ship in the base image).
# openxlsx backs the multi-sheet .xlsx downloads (shiny_app/download_helpers.R).
RUN install2.r --error --skipinstalled \
      bslib ggplot2 Matrix plotly DT svglite shinycssloaders remotes \
      openxlsx

# presto is GitHub-only; the interactive "Subset & DEGs" tab uses it for Wilcoxon DE.
RUN R -e 'remotes::install_github("immunogenomics/presto", upgrade = "never")'

WORKDIR /srv/shiny-app

# App code + the data-prep build scripts (scripts aren't needed at runtime, kept for reference).
# download_helpers.R IS needed at runtime -- app.R sources it.
COPY shiny_app/app.R shiny_app/download_helpers.R ./
COPY shiny_app/build_communication.R shiny_app/build_refmap.R \
     shiny_app/build_signature_scores.R shiny_app/build_subcluster_enrichment.R \
     shiny_app/build_fourgroup.R ./

# app_data.rds (the ~103 MB enriched bundle) is git-ignored and NOT baked in by default —
# it is supplied at runtime via a volume mount (see docker-compose.yml / DOCKER.md).
# For a fully self-contained image, place the bundle at shiny_app/app_data.rds and
# uncomment the next line (then no runtime mount is needed):
# COPY shiny_app/app_data.rds ./

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/srv/shiny-app', host = '0.0.0.0', port = 3838)"]
