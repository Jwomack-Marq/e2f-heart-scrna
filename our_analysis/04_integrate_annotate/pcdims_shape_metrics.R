#!/usr/bin/env Rscript
# Scale-free "did the embedding actually change?" metrics, appended to the sweep's
# comparison table.
#
# Why this exists as a second pass rather than being folded into pcdims_sweep.R: the
# sweep's own `umap_procrustes_rmsd` is n-DEPENDENT. It normalises both configurations to
# unit Frobenius norm, so two unrelated embeddings of n cells score sqrt(2/n) -- 0.0069 at
# n = 42,416 -- not 1. Read without that baseline, an rmsd of 0.0029 looks like "almost
# identical" when it is in fact 42% of the way to unrelated. That misreading already
# happened once in this project, so the fix is to publish metrics that cannot be misread:
#
#   procrustes_m2   standard Procrustes m^2, 0 = identical, 1 = unrelated. Scale-free.
#   dist_cor        correlation of pairwise distances on a 4,000-cell sample. Global shape.
#   knn30_kept      share of each cell's 30 nearest neighbours that are still neighbours
#                   in the other embedding. Local structure -- the one that decides whether
#                   "these cells sit together" survives a parameter change.
#
# m2 and knn30 answer different questions and can disagree: UMAP can rearrange clusters
# globally while keeping local neighbourhoods, or the reverse. Report both.
#
#   Rscript pcdims_shape_metrics.R --object=cm|combined
#   Rscript pcdims_shape_metrics.R --embeddings=<file.rds> --out=<comparison.csv>
#
# Reads processing/pcdims_<obj>_embeddings.rds, rewrites
# results/tables/pcdims_<obj>_comparison.csv with the extra rows appended.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))

argval <- function(flag, default = NULL) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
OBJ <- argval("--object", "cm")
EMB <- argval("--embeddings", file.path(PROC, sprintf("pcdims_%s_embeddings.rds", OBJ)))
OUT <- argval("--out",        file.path(OUTTAB, sprintf("pcdims_%s_comparison.csv", OBJ)))
stopifnot(file.exists(EMB), file.exists(OUT))

procrustes_m2 <- function(X, Y) {
  X <- scale(X, TRUE, FALSE); Y <- scale(Y, TRUE, FALSE)
  X <- X / sqrt(sum(X^2)); Y <- Y / sqrt(sum(Y^2))
  max(0, 1 - (sum(svd(crossprod(X, Y))$d))^2)      # clamp: -1e-15 is 0, not a negative
}
dist_cor <- function(X, Y, n = 4000L) {
  set.seed(1); i <- sample(nrow(X), min(n, nrow(X)))
  suppressWarnings(stats::cor(dist(X[i, ]), dist(Y[i, ])))
}
knn_kept <- function(X, Y, k = 30L, n = 3000L) {
  set.seed(1); i <- sample(nrow(X), min(n, nrow(X)))
  nn <- function(M) { D <- as.matrix(dist(M[i, ]))
                      t(apply(D, 1, function(v) order(v)[2:(k + 1)])) }
  a <- nn(X); b <- nn(Y)
  mean(vapply(seq_along(i), function(j) length(intersect(a[j, ], b[j, ])) / k, numeric(1)))
}

r    <- readRDS(EMB)
keys <- grep("^dims", names(r), value = TRUE)
dims <- vapply(keys, function(k) r[[k]]$dims, numeric(1))
n    <- nrow(r[[keys[1]]]$umap)
cmp  <- read.csv(OUT, stringsAsFactors = FALSE, check.names = FALSE)
cmp  <- cmp[!grepl("^(procrustes_m2|dist_cor|knn30_kept|rmsd_unrelated)$", cmp$metric), ]

row <- function(pair, metric, value, note = "") data.frame(
  pair = pair, metric = metric, value = value, note = note, stringsAsFactors = FALSE)
out <- list(row("(all)", "rmsd_unrelated_baseline", sqrt(2 / n),
                sprintf("two UNRELATED embeddings of %d cells score this on umap_procrustes_rmsd; 0 = identical", n)))
cat(sprintf("[%s] n = %d; unrelated-baseline rmsd = %.5f\n", OBJ, n, sqrt(2 / n)))
cat(sprintf("  %-14s %8s %10s %12s\n", "pair", "m2", "dist-cor", "30NN kept"))
for (p in list(c(1, 2), c(2, 3), c(1, 3))) {
  A <- r[[keys[p[1]]]]$umap; B <- r[[keys[p[2]]]]$umap
  tag <- sprintf("dims %d vs %d", dims[p[1]], dims[p[2]])
  m2 <- procrustes_m2(A, B); dc <- dist_cor(A, B); kk <- knn_kept(A, B)
  cat(sprintf("  %-14s %8.4f %10.4f %11.1f%%\n", tag, m2, dc, 100 * kk))
  out <- c(out, list(
    row(tag, "procrustes_m2", m2, "0 = identical, 1 = unrelated; scale-free"),
    row(tag, "dist_cor",      dc, "correlation of pairwise distances, 4,000-cell sample"),
    row(tag, "knn30_kept",    kk, "share of each cell's 30 nearest neighbours also neighbours in the other")))
}
write.csv(rbind(cmp, do.call(rbind, out)), OUT, row.names = FALSE)
cat(sprintf("\nappended %d rows to %s\n", length(out), basename(OUT)))
