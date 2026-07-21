# syntax=docker/dockerfile:1
FROM rocker/shiny:4.4.1

# System libs commonly needed by ggplot2/plotly/DT (most already in rocker/shiny;
# kept explicit so the build is robust across base-image changes).
RUN apt-get update && apt-get install -y --no-install-recommends \
      libxml2-dev libcurl4-openssl-dev libssl-dev \
      libpng-dev libfontconfig1-dev libfreetype6-dev \
    && rm -rf /var/lib/apt/lists/*

# App R packages (shiny/htmlwidgets/crosstalk are in the base; --skipinstalled avoids rework).
RUN install2.r --error --skipinstalled \
      bslib ggplot2 Matrix plotly DT crosstalk htmlwidgets

WORKDIR /srv/shiny-app
COPY shiny_app/app.R shiny_app/app_data.rds ./

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/srv/shiny-app', host = '0.0.0.0', port = 3838)"]
