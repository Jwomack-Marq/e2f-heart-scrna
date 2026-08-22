# 05_mito_sensitivity.R
# ---------------------------------------------------------------------------
# The KO-up list of EVERY subcluster contains 3-7 mitochondrially-encoded genes
# (mt-Nd1, mt-Nd2, mt-Nd4, mt-Cytb, mt-Co1, mt-Rnr1/2, mt-Atp8), and NO cluster
# has a single mt- gene on its KO-down side. Seven independent subclusters do not
# agree that perfectly by biology; a one-directional, every-cluster shift in the
# mt- genes is the signature of a difference in mitochondrial read fraction
# between the two libraries. It matters here because those genes ARE the
# oxidative-phosphorylation / electron-transport GO terms: leave them in and the
# headline KO-up result is partly a QC covariate wearing a pathway's name.
#
# So: measure the mito fraction per arm, then re-run GO BP on the primary tables
# with mt- genes removed, and ship both. A term that survives is about nuclear
# genes; a term that vanishes was carried by the mt- block.
# ---------------------------------------------------------------------------
suppressMessages({ library(Matrix); library(clusterProfiler); library(org.Mm.eg.db) })

OUT <- "/out"
app <- readRDS("/in/app_data.rds")
de  <- readRDS(file.path(OUT, "de_tables.rds"))
en  <- readRDS(file.path(OUT, "enrich.rds"))
P   <- de$params
man <- de$manifest

GO_PCUT <- en$params$go_pcut; GO_QCUT <- en$params$go_qcut
MIN_PCT <- en$params$min_pct_universe; MIN_INPUT <- en$params$min_input

# ---- 1. how different is the mitochondrial fraction between the arms? -------
# app$deg_expr is log-normalised, so this is a share of log-normalised signal, not
# a true read fraction -- fine for a relative comparison between arms on the same
# matrix, which is all that is claimed.
M  <- app$deg_expr
MD <- app$deg_meta[match(colnames(M), app$deg_meta$cell), , drop = FALSE]
MD$fourgrp <- paste(MD$genotype, MD$timepoint, sep = "-")
mt <- grep("^mt-", rownames(M), value = TRUE)
share <- Matrix::colSums(M[mt, , drop = FALSE]) / Matrix::colSums(M)
cat(sprintf("%d mitochondrially-encoded genes in the matrix\n\n", length(mt)))

frac <- do.call(rbind, lapply(c("AllCM", P$clusters), function(cl) {
  inCl <- if (cl == "AllCM") !is.na(MD$cm_subcluster) else
          (!is.na(MD$cm_subcluster) & MD$cm_subcluster == cl)
  do.call(rbind, lapply(c("WT-P0", "WT-P7", "KO-P0", "KO-P7"), function(g) {
    i <- which(inCl & MD$fourgrp == g)
    if (!length(i)) return(NULL)
    data.frame(cluster = cl, group = g, n_cells = length(i),
               mt_share_pct = round(100 * mean(share[i]), 2), stringsAsFactors = FALSE)
  }))
}))
wide <- reshape(frac[, c("cluster","group","mt_share_pct")], idvar = "cluster",
                timevar = "group", direction = "wide")
names(wide) <- sub("^mt_share_pct\\.", "", names(wide))
wide$P7_KO_minus_WT <- round(wide$`KO-P7` - wide$`WT-P7`, 2)
wide$P7_KO_over_WT  <- round(wide$`KO-P7` / wide$`WT-P7`, 2)
cat("== mitochondrial share of log-normalised signal, % ==\n")
print(wide, row.names = FALSE)
write.csv(frac, file.path(OUT, "csv", "_mito_fraction.csv"), row.names = FALSE)

# ---- 2. GO BP on the primary tables with mt- genes removed ------------------
pct_cols <- function(d) grep("^pct_", names(d), value = TRUE)
universe_of <- function(d, drop_mt) {
  pc <- pct_cols(d)
  u <- unique(d$gene[do.call(pmax, c(d[pc], list(na.rm = TRUE))) >= MIN_PCT])
  if (drop_mt) setdiff(u, mt) else u
}
run_go <- function(genes, universe) {
  if (length(genes) < MIN_INPUT) return(NULL)
  eg <- tryCatch(suppressWarnings(suppressMessages(
    enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
             universe = unique(universe), pAdjustMethod = "BH",
             pvalueCutoff = GO_PCUT, qvalueCutoff = GO_QCUT,
             minGSSize = 10, maxGSSize = 500, readable = FALSE))), error = function(e) NULL)
  if (is.null(eg)) return(NULL)
  df <- as.data.frame(eg); if (!nrow(df)) NULL else df
}

rows <- list(); cmp <- list()
prim <- man[man$is_primary, ]
for (i in seq_len(nrow(prim))) {
  key <- prim$key[i]; d <- de$tables[[key]]
  uni <- universe_of(d, TRUE)
  for (dirn in c(prim$up_label[i], prim$down_label[i])) {
    sgn <- if (dirn == prim$up_label[i]) 1 else -1
    g <- d$gene[!d$confounder & is.finite(d$padj) & d$padj < P$sig_padj &
                sgn * d$log2FoldChange >= P$sig_lfc]
    n_mt <- sum(g %in% mt); g <- setdiff(g, mt)
    res <- run_go(g, uni)
    nr <- if (is.null(res)) 0L else nrow(res)
    if (!is.null(res)) {
      res$contrast <- prim$contrast[i]; res$cluster <- prim$cluster[i]
      res$stratum <- prim$stratum[i]; res$direction <- dirn
      res$n_input_no_mt <- length(g); res$n_mt_removed <- n_mt
      rows[[paste(key, dirn)]] <- res
    }
    # what had the mt- block been carrying?
    orig <- if (is.null(en$go)) NULL else en$go[en$go$contrast == prim$contrast[i] &
              en$go$cluster == prim$cluster[i] & en$go$stratum == prim$stratum[i] &
              en$go$ontology == "BP" & en$go$direction == dirn, ]
    o <- if (is.null(orig)) character() else orig$Description
    n2 <- if (is.null(res)) character() else res$Description
    cmp[[length(cmp) + 1]] <- data.frame(
      contrast = prim$contrast[i], cluster = prim$cluster[i], direction = dirn,
      n_sig_genes = length(g) + n_mt, n_mt_genes = n_mt,
      terms_with_mt = length(o), terms_without_mt = nr,
      terms_lost = length(setdiff(o, n2)), terms_retained = length(intersect(o, n2)),
      example_lost = paste(head(setdiff(o, n2), 3), collapse = "; "),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-30s %-8s  sig=%3d (mt=%d)  GO terms %3d -> %3d\n",
                key, dirn, length(g) + n_mt, n_mt, length(o), nr))
  }
}

cmp <- do.call(rbind, cmp)
go_nomt <- if (length(rows)) do.call(rbind, rows) else NULL
saveRDS(list(go_no_mt = go_nomt, comparison = cmp, mito_fraction = frac,
             mito_wide = wide, mt_genes = mt),
        file.path(OUT, "mito_sensitivity.rds"))
write.csv(cmp, file.path(OUT, "csv", "_mito_go_comparison.csv"), row.names = FALSE)
if (!is.null(go_nomt)) write.csv(go_nomt, file.path(OUT, "csv", "_go_BP_no_mt.csv"), row.names = FALSE)

cat("\n== GO terms lost when the mt- genes are removed ==\n")
print(cmp[, c("cluster","direction","n_sig_genes","n_mt_genes","terms_with_mt","terms_without_mt","terms_lost")],
      row.names = FALSE)
cat("\nDONE 05_mito_sensitivity.R\n")
