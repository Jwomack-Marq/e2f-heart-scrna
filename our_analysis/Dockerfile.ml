# syntax=docker/dockerfile:1
# Image for the KO-signature machine-learning script (05_analyses/ko_signature_ml.R).
#
# WHY A SEPARATE IMAGE. glmnet is in no image in this repo -- checked against every
# env/*.tsv capture -- so the shipped results/models/*_glmnet.rds could not be
# regenerated anywhere on this machine (they came off the Windows box, like the DESeq2
# tables before Dockerfile.seurat grew apeglm). The KO script needs glmnet for the same
# elastic-net logistic regression cell_state_classifier.R uses, plus everything the
# base already carries: Seurat to read the objects, Matrix for sparse features,
# msigdbr (offline cache) for the Hallmark E2F-target set, org.Mm.eg.db for the
# chromosome column in the removed-genes audit, babelgene for cc_lists().
#
# FROM e2f-seurat-full rather than e2f-enrich because the script reads the annotated
# Seurat object and calls NormalizeData / FindVariableFeatures on matrices.
#
# VERSION CAVEAT inherited from the base: this image does not reproduce the package
# versions the shipped objects or the shipped glmnet models were built with. The KO
# script fits every model it reports inside this image, so its results are internally
# consistent; the ONE shipped model it applies (cm_maturation_glmnet.rds) is only
# deserialised and multiplied out, which glmnet's predict() does identically across
# versions for a fixed lambda.
FROM e2f-seurat-full:latest

# CRAN side. p3m.dev serves noble binaries, so this is a download, not a build.
# CRAN-latest is appended (not substituted), exactly as Dockerfile.seurat does.
# install.packages() only WARNS on failure, so assert afterwards.
RUN R -q -e 'options(repos = c(getOption("repos"), CRAN = "https://cloud.r-project.org")); \
      install.packages("glmnet", Ncpus = 4); \
      if (!"glmnet" %in% rownames(installed.packages())) stop("CRAN install failed for: glmnet")'

# Fail the build now rather than 20 minutes into a run. A REAL fit, not library():
# sparse dgCMatrix input, keep = TRUE (the script reads out-of-fold predictions from
# fit.preval), and lambda.1se selection are the three things the script depends on.
RUN R -q -e 'suppressMessages({ library(glmnet); library(Matrix) }); \
      set.seed(1); m <- matrix(rpois(400 * 300, 0.6), nrow = 400, \
                               dimnames = list(paste0("c", 1:400), paste0("g", 1:300))); \
      X <- Matrix::Matrix(m, sparse = TRUE); stopifnot(inherits(X, "dgCMatrix")); \
      y <- factor(rep(c("WT","KO"), each = 200), levels = c("WT","KO")); \
      fit <- cv.glmnet(X, y, family = "binomial", alpha = 0.9, nfolds = 5, \
                       type.measure = "deviance", keep = TRUE); \
      j <- match(fit$lambda.1se, fit$lambda); \
      stopifnot(is.matrix(fit$fit.preval), nrow(fit$fit.preval) == 400, !anyNA(fit$fit.preval[, j])); \
      p <- predict(fit, X, s = "lambda.1se", type = "response"); \
      stopifnot(all(p >= 0 & p <= 1)); \
      cat("glmnet", as.character(packageVersion("glmnet")), "| sparse cv fit ok | nzero at 1se:", fit$nzero[j], "\n")'

# The offline gene-set path: msigdbr must resolve from the cache the base image warmed
# under R_USER_CACHE_DIR (containers run with --user and an unwritable HOME), and
# org.Mm.eg.db must map a symbol to its chromosome.
RUN R -q -e 'h <- msigdbr::msigdbr(species = "Mus musculus", collection = "H"); \
      e2f <- unique(h$gene_symbol[h$gs_name == "HALLMARK_E2F_TARGETS"]); \
      stopifnot(length(e2f) > 150); \
      chr <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = "Xist", \
                                   column = "CHR", keytype = "SYMBOL"); \
      stopifnot(identical(unname(chr), "X")); \
      cat("HALLMARK_E2F_TARGETS:", length(e2f), "mouse symbols | Xist on chr", chr, "\n")'

WORKDIR /work
