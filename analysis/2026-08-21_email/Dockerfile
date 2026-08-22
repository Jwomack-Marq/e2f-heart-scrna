# syntax=docker/dockerfile:1
# Enrichment image for the 2026-08-21 collaborator deliverable.
#
# Built FROM the deployed app image, which is Ubuntu 24.04 / R 4.5.1 with its
# repos already pointed at the p3m.dev noble BINARY mirror, and which already
# carries Matrix + presto (the DE core in 01_de.R is the same presto::wilcoxauc
# call the app's "Subset & DEGs" tab and build_fourgroup.R use).
#
# What this adds is the enrichment stack that build_subcluster_enrichment.R
# declares but that the app image deliberately does not ship, plus openxlsx for
# the workbooks (the repo has no Excel code anywhere -- this is greenfield).
FROM lab-server-e2f-heart-scrna-dev:latest

# clusterProfiler/enrichplot pull in ggraph -> igraph, which needs glpk + gmp;
# ragg needs the freetype/harfbuzz stack (already present, kept explicit).
RUN apt-get update && apt-get install -y --no-install-recommends \
      libglpk-dev libgmp-dev libxml2-dev libcurl4-openssl-dev libssl-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff5-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN side. p3m.dev serves noble binaries, so this is a download, not a build.
RUN install2.r --error --skipinstalled \
      BiocManager openxlsx ggplot2 ragg ggrepel scales viridisLite msigdbr

# The base image pins its CRAN mirror to the p3m noble snapshot of 2025-10-30,
# which is older than current clusterProfiler needs: it ships yulab.utils 0.2.1 and
# ggtangle 0.0.7 against requirements of >= 0.2.3 / >= 0.0.9, so enrichplot and
# clusterProfiler both fail to install with only a *warning*, not an error.
# Pull just those pure-R helpers from CRAN-latest (source, but tiny). p3m stays
# first in the list so everything else still resolves to a noble binary.
RUN R -q -e 'options(repos = c(getOption("repos"), CRAN = "https://cloud.r-project.org")); install.packages(c("yulab.utils","ggtangle","ggfun","tidytree","aplot","gson"), Ncpus = 4)'

# Bioconductor side. org.Mm.eg.db is the mouse annotation DB enrichGO() reads through.
# BiocManager::install() only warns on a failed package, so assert afterwards in the
# same layer -- otherwise a broken image caches and surfaces 40 minutes into a run.
RUN R -q -e 'options(repos = c(getOption("repos"), CRAN = "https://cloud.r-project.org")); \
      BiocManager::install(c("clusterProfiler","org.Mm.eg.db","enrichplot","fgsea","AnnotationDbi"), ask = FALSE, update = FALSE, Ncpus = 4); \
      miss <- setdiff(c("clusterProfiler","org.Mm.eg.db","enrichplot","fgsea","AnnotationDbi"), rownames(installed.packages())); \
      if (length(miss)) stop("Bioconductor install failed for: ", paste(miss, collapse = ", "))'

# Fail the build now rather than 40 minutes into an analysis run.
RUN R -q -e 'for (p in c("Matrix","presto","clusterProfiler","org.Mm.eg.db","enrichplot","fgsea","msigdbr","openxlsx","ggplot2","ggrepel","ragg")) { library(p, character.only = TRUE); cat("ok ", p, " ", as.character(packageVersion(p)), "\n", sep = "") }'

# msigdbr >= 25.x fetches the MSigDB release at first use and caches it under
# R_user_dir(), which is HOME-relative. The analysis containers run with --user so
# that outputs are not root-owned, and that leaves HOME unwritable -- the fetch then
# dies with "Failed to open file //.cache/R/msigdbr/...". Pin the cache to a
# world-readable path and warm it at build time, so the run needs no network at all.
ENV R_USER_CACHE_DIR=/opt/rcache
RUN mkdir -p /opt/rcache \
 && R -q -e 'invisible(msigdbr::msigdbr(species = "Mus musculus", collection = "H")); \
             invisible(msigdbr::msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")); \
             cat("msigdbr cache warmed\n")' \
 && chmod -R a+rX /opt/rcache

WORKDIR /work
