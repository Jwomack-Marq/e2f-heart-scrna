# syntax=docker/dockerfile:1
# Renderer for the methods book (docs/).
#
# Deliberately NOT built on e2f-enrich:latest. Every chapter declares
# `engine: markdown`, so no code chunk is ever executed at render time and the
# renderer needs neither R nor app_data.rds -- only the Quarto CLI. Keeping the
# two images apart means a prose fix rebuilds in seconds against a 250 MB image
# instead of a 4.3 GB one, and it makes the no-execution guarantee structural
# rather than a promise: R is not installed here, so a stray executable chunk
# fails the render instead of silently working on one machine.
#
# Figures come from docs/export_assets.R, which DOES need the R image -- see
# docs/export.sh. That step runs once, offline, and its output is committed.
FROM debian:bookworm-slim

ARG QUARTO_VERSION=1.10.18

# ca-certificates for the download; git because Quarto probes for it when
# resolving project metadata and warns loudly otherwise.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl git \
    && curl -fsSL -o /tmp/quarto.deb \
       "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && dpkg -i /tmp/quarto.deb \
    && rm -f /tmp/quarto.deb \
    && apt-get purge -y curl && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Fail the build now rather than at first render.
RUN quarto --version

# Quarto writes its scratch under XDG paths; the render runs with --user so that
# _site/ is not root-owned, which leaves $HOME unwritable (same trap the
# analysis image documents for the msigdbr cache). Point both at /tmp.
ENV XDG_CACHE_HOME=/tmp/.cache \
    XDG_DATA_HOME=/tmp/.local/share

WORKDIR /docs
