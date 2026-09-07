#!/usr/bin/env Rscript
# Cell-state predictor — DESCRIPTIVE, cell-level ML. Genotype is deliberately NOT
# predicted: n=1 per condition, genotype confounded with sex, KO not transcript-
# confirmed -> a genotype classifier would learn sex/animal/batch, not knockout
# biology. Cell identity/state is well-powered at the cell level and confound-free.
#
# Trains two portable, interpretable elastic-net (glmnet) models on the annotated
# atlas and saves them so future/external data can be labelled consistently:
#   1) cell-type classifier (multinomial)            -> results/models/cell_type_glmnet.rds
#   2) CM developmental-stage predictor (P0 vs P7)    -> results/models/cm_maturation_glmnet.rds
# Plus held-out accuracy + confusion, interpretable gene panels, a reusable
# predict_cell_state() applier, and an optional out-of-sample check on the Baniol
# (PRJEB47622) cardiomyocytes.
#
# Note on the CM target: pseudotime isn't stored on this object, and regressing the
# mature-minus-immature module score on its own constituent genes would be circular.
# So the stage model predicts TIMEPOINT (P0 vs P7) — an external label, and NOT sex-
# confounded (both sexes are present at each timepoint) — then we confirm its
# P7-probability tracks the maturation module score.
#
# Reads seurat.combined.annotated.rds; writes results/{models,tables,figures}/.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({ library(glmnet); library(Matrix); library(ggplot2) }))
set.seed(1)
MODELS <- file.path(RESULTS, "models"); if (!dir.exists(MODELS)) dir.create(MODELS, recursive = TRUE)
say <- function(...) { cat(sprintf(...), "\n"); flush(stdout()) }

## ---- 1. load + de-duplicate the lane double-counting -----------------------
comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds"))
# lane1/lane6 are the SAME library sequenced twice (~97-100% barcode overlap). Keep
# ONE lane so a cell's twin can't leak across CV folds; unbiased ~half-depth subsample.
comb <- comb[, comb$lane == "lane1"]
DefaultAssay(comb) <- "RNA"
comb <- NormalizeData(comb, verbose = FALSE)
comb <- FindVariableFeatures(comb, nfeatures = 2000, verbose = FALSE)
hvg  <- VariableFeatures(comb)
lognorm <- GetAssayData(comb, assay = "RNA", layer = "data")
say("after lane-dedup: %d cells x %d genes; %d HVGs used as features",
    ncol(comb), nrow(comb), length(hvg))

## ---- 2. cell-type classifier (multinomial glmnet) --------------------------
y <- factor(comb$celltype)
keep_cls <- names(which(table(y) >= 30))                 # drop ultra-rare classes for stable CV
sel  <- y %in% keep_cls
Xall <- Matrix::t(lognorm[hvg, sel, drop = FALSE])       # cells x genes (sparse)
yall <- droplevels(y[sel])
say("cell-type model: %d cells, %d classes: %s", nrow(Xall), nlevels(yall),
    paste(levels(yall), collapse = ", "))

# honest held-out evaluation via a stratified 80/20 split
te_i <- unlist(lapply(split(seq_along(yall), yall), function(i) sample(i, ceiling(0.2 * length(i)))))
te <- seq_along(yall) %in% te_i; tr <- !te
fit_tr  <- cv.glmnet(Xall[tr, ], yall[tr], family = "multinomial", alpha = 0.9,
                     type.measure = "class", nfolds = 5)
pred_te <- factor(as.vector(predict(fit_tr, Xall[te, ], s = "lambda.1se", type = "class")),
                  levels = levels(yall))
conf <- table(true = yall[te], predicted = pred_te)
acc  <- sum(diag(conf)) / sum(conf)
recall <- diag(conf) / rowSums(conf)
say("cell-type held-out accuracy: %.3f", acc)
write.csv(as.data.frame.matrix(conf), file.path(OUTTAB, "celltype_confusion.csv"))
write.csv(data.frame(celltype = names(recall), recall = round(recall, 3), n_test = rowSums(conf)),
          file.path(OUTTAB, "celltype_recall.csv"), row.names = FALSE)

# portable model on ALL cells + interpretable marker panel (non-zero coefficients)
fit_all <- cv.glmnet(Xall, yall, family = "multinomial", alpha = 0.9,
                     type.measure = "class", nfolds = 5)
co <- coef(fit_all, s = "lambda.1se")
panel <- do.call(rbind, lapply(names(co), function(cl) {
  m <- as.matrix(co[[cl]]); nz <- which(m[, 1] != 0 & rownames(m) != "(Intercept)")
  if (!length(nz)) return(NULL)
  data.frame(celltype = cl, gene = rownames(m)[nz], weight = round(m[nz, 1], 4))
}))
panel <- panel[order(panel$celltype, -abs(panel$weight)), ]
write.csv(panel, file.path(OUTTAB, "celltype_marker_panel.csv"), row.names = FALSE)
ct_bundle <- list(model = fit_all, features = hvg, family = "multinomial", classes = levels(yall))
saveRDS(ct_bundle, file.path(MODELS, "cell_type_glmnet.rds"))

cdf <- as.data.frame(prop.table(conf, 1)); names(cdf) <- c("true", "predicted", "frac")
ggsave(file.path(OUTFIG, "celltype_confusion.png"),
  ggplot(cdf, aes(predicted, true, fill = frac)) + geom_tile() +
    geom_text(aes(label = ifelse(frac > 0.01, sprintf("%.2f", frac), "")), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
    labs(title = sprintf("Cell-type classifier — held-out (acc = %.2f)", acc),
         x = "predicted", y = "true", fill = "row frac"),
  width = 7, height = 6, dpi = 120)

## ---- 3. CM developmental-stage predictor (P0 vs P7) ------------------------
cm <- comb[, comb$celltype == "Cardiomyocyte"]
cm <- AddModuleScore(cm, features = list(intersect(CM_MATURE, rownames(cm))),   name = "matP", seed = 1)
cm <- AddModuleScore(cm, features = list(intersect(CM_IMMATURE, rownames(cm))), name = "immP", seed = 1)
mat_score <- cm$matP1 - cm$immP1
Xcm <- Matrix::t(GetAssayData(cm, assay = "RNA", layer = "data")[hvg, , drop = FALSE])
ycm <- factor(cm$timepoint, levels = c("P0", "P7"))
say("CM stage model: %d CMs (P0=%d, P7=%d)", ncol(cm), sum(ycm == "P0"), sum(ycm == "P7"))

tec_i <- unlist(lapply(split(seq_along(ycm), ycm), function(i) sample(i, ceiling(0.2 * length(i)))))
tec <- seq_along(ycm) %in% tec_i; trc <- !tec
fitc_tr  <- cv.glmnet(Xcm[trc, ], ycm[trc], family = "binomial", alpha = 0.9,
                      type.measure = "class", nfolds = 5)
predc_te <- factor(as.vector(predict(fitc_tr, Xcm[tec, ], s = "lambda.1se", type = "class")),
                   levels = c("P0", "P7"))
confc <- table(true = ycm[tec], predicted = predc_te); accc <- sum(diag(confc)) / sum(confc)
say("CM P0-vs-P7 held-out accuracy: %.3f", accc)

fitc_all <- cv.glmnet(Xcm, ycm, family = "binomial", alpha = 0.9, type.measure = "class", nfolds = 5)
p7prob <- as.vector(predict(fitc_all, Xcm, s = "lambda.1se", type = "response"))
rho <- suppressWarnings(cor(p7prob, mat_score, method = "spearman"))
say("CM: Spearman(P7-probability, maturation module score) = %.3f", rho)

cc <- as.matrix(coef(fitc_all, s = "lambda.1se")); cc <- cc[rownames(cc) != "(Intercept)", , drop = FALSE]
nz <- cc[cc[, 1] != 0, , drop = FALSE]
stage_genes <- data.frame(gene = rownames(nz), weight = round(nz[, 1], 4))
stage_genes <- stage_genes[order(-stage_genes$weight), ]     # + toward P7 (mature), - toward P0
write.csv(stage_genes, file.path(OUTTAB, "cm_stage_genes.csv"), row.names = FALSE)
cm_bundle <- list(model = fitc_all, features = hvg, family = "binomial", positive = "P7")
saveRDS(cm_bundle, file.path(MODELS, "cm_maturation_glmnet.rds"))

sdf <- data.frame(p7prob = p7prob, mat_score = mat_score, timepoint = ycm)
ggsave(file.path(OUTFIG, "cm_stage_vs_maturation.png"),
  ggplot(sdf, aes(mat_score, p7prob, color = timepoint)) + geom_point(alpha = .3, size = .6) +
    scale_color_manual(values = c(P0 = "#1565c0", P7 = "#c62828")) +
    theme_minimal(base_size = 12) +
    labs(title = sprintf("CM stage model (acc = %.2f, rho = %.2f)", accc, rho),
         x = "maturation score (mature - immature)", y = "predicted P(P7)"),
  width = 7, height = 5, dpi = 120)

## ---- 4. portable applier (the reusable "prediction" tool) ------------------
predict_cell_state <- function(bundle, newmat) {
  # newmat: genes (rows) x cells (cols), log-normalized. Missing features -> 0.
  Xn <- matrix(0, nrow = ncol(newmat), ncol = length(bundle$features),
               dimnames = list(colnames(newmat), bundle$features))
  common <- intersect(rownames(newmat), bundle$features)
  Xn[, common] <- t(as.matrix(newmat[common, , drop = FALSE]))
  if (bundle$family == "binomial")
    setNames(as.vector(predict(bundle$model, Xn, s = "lambda.1se", type = "response")), rownames(Xn))
  else
    setNames(as.vector(predict(bundle$model, Xn, s = "lambda.1se", type = "class")), rownames(Xn))
}
saveRDS(predict_cell_state, file.path(MODELS, "predict_cell_state.fun.rds"))

## ---- 5. optional out-of-sample check: Baniol cardiomyocytes ----------------
ban_rds <- "C:/Users/Justi/OneDrive/Documents/GitHub/rna_heart_test/comparison_analysis/processing/seurat.baniol.rds"
if (file.exists(ban_rds)) {
  ban <- readRDS(ban_rds); DefaultAssay(ban) <- "RNA"; ban <- NormalizeData(ban, verbose = FALSE)
  bl <- GetAssayData(ban, assay = "RNA", layer = "data")
  ov <- length(intersect(rownames(bl), hvg))
  say("Baniol external check: %d cells, %d/%d HVGs present", ncol(ban), ov, length(hvg))
  if (ov >= 0.4 * length(hvg)) {
    ct_pred <- predict_cell_state(ct_bundle, bl)
    say("  Baniol cell-type calls: %s",
        paste(names(table(ct_pred)), table(ct_pred), sep = "=", collapse = "  "))
    p7 <- predict_cell_state(cm_bundle, bl)
    tp_col <- intersect(c("timepoint", "stage", "Stage", "age"), colnames(ban@meta.data))
    if (length(tp_col)) {
      agg <- tapply(p7, ban@meta.data[[tp_col[1]]], median)
      say("  Baniol median P(P7) by %s: %s", tp_col[1],
          paste(names(agg), round(agg, 2), sep = "=", collapse = "  "))
    }
  } else say("  gene overlap too low (Smart-seq2 vs 3'); skipping external prediction")
} else say("Baniol object not found; skipping external check")

say("=== DONE cell_state_classifier ===")
