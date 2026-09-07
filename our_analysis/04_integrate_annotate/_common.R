#!/usr/bin/env Rscript
# Forwarder. The real shared library lives at the our_analysis/ root (_common.R).
# Step scripts source "_common.R" from their own folder via dirname(this); this
# stub walks up to the `.projroot` sentinel and sources the real file there, so
# the shared code is never duplicated and works at any depth.
local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- normalizePath(if (length(f) && nzchar(f)) dirname(f) else getwd())
  while (!file.exists(file.path(d, ".projroot"))) {
    up <- dirname(d)
    if (identical(up, d)) stop("our_analysis/.projroot not found above ", d)
    d <- up
  }
  source(file.path(d, "_common.R"))
})
