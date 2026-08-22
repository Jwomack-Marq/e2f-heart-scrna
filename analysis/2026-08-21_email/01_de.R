# 01_de.R -- differential expression for the 2026-08-21 collaborator request.
# ---------------------------------------------------------------------------
# Two contrasts, computed UNGATED (every gene in the matrix gets a row, so the
# collaborator gets "complete gene lists") and genome-wide:
#
#   Part 1  P7_KO_vs_WT   KO-P7 vs WT-P7, within CM1/2/3/4/5/7/8
#           Both arms are P7, so the FACS cycling-enrichment is matched and the
#           all-cells stratum is the primary read. A G1-only table ships too.
#
#   Part 2  WT_P7_vs_P0   WT-P7 vs WT-P0, over all CM plus the same 7 subclusters
#           P7 was FACS cycling-enriched 4.5-5.2x and P0 essentially not, so the
#           raw contrast reads out the sort as much as development. G1-matched is
#           the primary read here; the raw one ships labelled.
#
# The DE core is the same presto::wilcoxauc call the app's "Subset & DEGs" tab
# uses (deg_compute(), app.R:566) and that build_fourgroup.R's de_one() uses, so
# these numbers are directly comparable to what the website already shows -- the
# only deliberate difference is that the row gate is removed.
#
# Run:  docker run --rm -v "$PWD/shiny_app/app_data.rds:/in/app_data.rds:ro" \
#                       -v "$PWD/deliverables/2026-08-21:/out" \
#                       -v "$PWD:/repo:ro" e2f-enrich:latest \
#                       Rscript /repo/analysis/2026-08-21_email/01_de.R
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))
stopifnot(requireNamespace("presto", quietly = TRUE))

IN    <- "/in/app_data.rds"
OUT   <- "/out"
REPO  <- "/repo"
dir.create(file.path(OUT, "csv"), recursive = TRUE, showWarnings = FALSE)

RES_COL   <- "SCT_snn_res.0.2"     # the only resolution the app exposes (app.R:23)
MIN_CELLS <- 10                    # same floor as deg_compute() and de_one()
SIG_PADJ  <- 0.05
SIG_LFC   <- 0.25
CLUSTERS  <- paste0("CM", c(1, 2, 3, 4, 5, 7, 8))   # the collaborator's list

cat("Loading", IN, "...\n")
app  <- readRDS(IN)
CONF <- app$confound
stopifnot(!is.null(app$deg_expr), !is.null(app$deg_meta), !is.null(app$cm$meta))

# ---- gene sets ------------------------------------------------------------
# Parsed out of build_signature_scores.R rather than copied, so the two can't
# drift: pull the literal `SETS <- list(...)` and evaluate only that expression.
sets_from_builder <- function(path) {
  for (x in parse(path)) {
    if (is.call(x) && identical(as.character(x[[1]]), "<-") &&
        identical(as.character(x[[2]]), "SETS")) return(eval(x[[3]]))
  }
  stop("could not find the SETS definition in ", path)
}
SETS <- sets_from_builder(file.path(REPO, "shiny_app", "build_signature_scores.R"))
# Seurat's cc.genes lists, as embedded in model/tools/extract_ko_bundle.R (human
# symbols there); title-case for mouse, which is how the bundle's Phase calls
# were made upstream.
S_HUMAN <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
 "CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76",
 "SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1",
 "TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
G2M_HUMAN <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
 "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2",
 "AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20",
 "TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR",
 "AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
titlecase <- function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
SETS$phase_S   <- titlecase(S_HUMAN)
SETS$phase_G2M <- titlecase(G2M_HUMAN)
cat("signature sets:", paste(sprintf("%s(%d)", names(SETS), lengths(SETS)), collapse = " "), "\n")

# sig_sets column: which curated programs a gene belongs to. Built once as a
# lookup over the union of the sets (a few hundred genes) rather than scanning
# every set for each of 24,221 genes in each of ~28 tables.
SET_MAP <- local({
  long <- data.frame(gene = unlist(SETS, use.names = FALSE),
                     set  = rep(names(SETS), lengths(SETS)), stringsAsFactors = FALSE)
  vapply(split(long$set, long$gene), function(x) paste(unique(x), collapse = ";"), "")
})
set_membership <- function(genes) {
  out <- unname(SET_MAP[genes]); out[is.na(out)] <- ""; out
}

# ---- matrices -------------------------------------------------------------
# broad = the genome-wide one GO needs (24,221 genes) but downsampled to 8,026 cells.
# curated = 2,181 genes over all 30,030 cells; used only for a full-depth cross-check.
Mb  <- app$deg_expr
MDb <- app$deg_meta[match(colnames(Mb), app$deg_meta$cell), , drop = FALSE]
Mc  <- app$expr
MDc <- app$meta[match(colnames(Mc), app$meta$cell), , drop = FALSE]

cm_lookup <- setNames(paste0("CM", app$cm$meta[[RES_COL]]), app$cm$meta$cell)
for (nm in c("b", "c")) {
  md <- get(paste0("MD", nm))
  md$cm_subcluster <- unname(cm_lookup[md$cell])
  md$fourgrp <- paste(md$genotype, md$timepoint, sep = "-")
  assign(paste0("MD", nm), md)
}
cat(sprintf("broad   %d genes x %d cells (%d CM)\n", nrow(Mb), ncol(Mb), sum(!is.na(MDb$cm_subcluster))))
cat(sprintf("curated %d genes x %d cells (%d CM)\n", nrow(Mc), ncol(Mc), sum(!is.na(MDc$cm_subcluster))))

# ---- DE core --------------------------------------------------------------
# Identical to de_one() in build_fourgroup.R and deg_compute() in app.R, minus
# the row gate. logFC > 0 means up in A.
de_core <- function(M, idxA, idxB) {
  nA <- length(idxA); nB <- length(idxB)
  if (nA < MIN_CELLS || nB < MIN_CELLS) return(NULL)
  cols <- c(idxA, idxB)
  grp  <- rep(c("A", "B"), c(nA, nB))
  res  <- presto::wilcoxauc(M[, cols, drop = FALSE], grp)
  rA <- res[res$group == "A", ]; rB <- res[res$group == "B", ]
  m  <- match(rA$feature, rB$feature)
  data.frame(
    gene = rA$feature,
    log2FoldChange = round(rA$logFC, 4), auc = round(rA$auc, 4),
    pvalue = rA$pval, padj = rA$padj,
    pct_A = round(rA$pct_in, 1), pct_B = round(rA$pct_out, 1),
    mean_A = round(rA$avgExpr, 4), mean_B = round(rB$avgExpr[m], 4),
    n_A = nA, n_B = nB, stringsAsFactors = FALSE)
}

# cell indices for one arm of one contrast
arm_idx <- function(MD, cluster, group, stratum) {
  inCl <- if (cluster == "AllCM") !is.na(MD$cm_subcluster)
          else !is.na(MD$cm_subcluster) & MD$cm_subcluster == cluster
  ph <- if (identical(stratum, "G1")) as.character(MD$Phase) == "G1" else TRUE
  which(inCl & ph & MD$fourgrp == group)
}

CONTRASTS <- list(
  P7_KO_vs_WT = list(A = "KO-P7", B = "WT-P7", A_lab = "KO_P7", B_lab = "WT_P7",
                     up = "KO_up", down = "KO_down",
                     primary = "all", clusters = CLUSTERS,
                     note = "Both arms are P7, so the FACS cycling-enrichment is matched; the all-cells stratum is the primary read."),
  WT_P7_vs_P0 = list(A = "WT-P7", B = "WT-P0", A_lab = "WT_P7", B_lab = "WT_P0",
                     up = "P7_up", down = "P0_up",
                     primary = "G1", clusters = c("AllCM", CLUSTERS),
                     note = "P7 was FACS cycling-enriched 4.5-5.2x and P0 essentially not; the G1-matched stratum is the primary read and the all-cells one is sort-confounded."))

# ---- run ------------------------------------------------------------------
tables <- list(); skipped <- list(); manifest <- list()

for (ck in names(CONTRASTS)) {
  ct <- CONTRASTS[[ck]]
  for (cl in ct$clusters) for (st in c("all", "G1")) {
    key <- sprintf("%s__%s__%s", ck, cl, st)
    iA <- arm_idx(MDb, cl, ct$A, st); iB <- arm_idx(MDb, cl, ct$B, st)
    d  <- de_core(Mb, iA, iB)
    if (is.null(d)) {
      skipped[[length(skipped) + 1]] <- data.frame(
        contrast = ck, cluster = cl, stratum = st,
        n_A = length(iA), n_B = length(iB),
        reason = sprintf("fewer than %d cells in an arm", MIN_CELLS), stringsAsFactors = FALSE)
      cat(sprintf("  SKIP %-32s A=%-4d B=%-4d\n", key, length(iA), length(iB)))
      next
    }

    d$direction <- ifelse(d$padj < SIG_PADJ & d$log2FoldChange >=  SIG_LFC, ct$up,
                   ifelse(d$padj < SIG_PADJ & d$log2FoldChange <= -SIG_LFC, ct$down, "ns"))
    d$confounder <- d$gene %in% CONF
    d$sig_sets   <- set_membership(d$gene)
    # Genes detected in zero cells on both sides carry no information (presto still
    # returns them at logFC 0, p = 1). They stay in the CSVs but are marked so the
    # workbook can drop them instead of padding every sheet with ~5,000 empty rows.
    d$detected   <- pmax(d$pct_A, d$pct_B) > 0

    # full-depth cross-check on the curated panel (all cells, 2,181 genes)
    jA <- arm_idx(MDc, cl, ct$A, st); jB <- arm_idx(MDc, cl, ct$B, st)
    dc <- de_core(Mc, jA, jB)
    if (!is.null(dc)) {
      k <- match(d$gene, dc$gene)
      d$lfc_fullcells  <- dc$log2FoldChange[k]
      d$padj_fullcells <- dc$padj[k]
      d$n_A_fullcells  <- ifelse(is.na(k), NA_integer_, length(jA))
      d$n_B_fullcells  <- ifelse(is.na(k), NA_integer_, length(jB))
    } else {
      d$lfc_fullcells <- NA_real_; d$padj_fullcells <- NA_real_
      d$n_A_fullcells <- NA_integer_; d$n_B_fullcells <- NA_integer_
    }

    # what the website currently shows for this subcluster (KO vs WT POOLED over
    # P0+P7) -- carried for reference so the collaborator can see the difference.
    sd <- app$tables$sub_DE[["res0.2"]][[cl]]
    if (ck == "P7_KO_vs_WT" && !is.null(sd)) {
      k <- match(d$gene, sd$gene)
      d$lfc_pooled_website  <- sd$log2FoldChange[k]
      d$padj_pooled_website <- sd$padj[k]
    }

    # rename the generic A/B columns to the actual arms
    ren <- c(pct_A = paste0("pct_", ct$A_lab), pct_B = paste0("pct_", ct$B_lab),
             mean_A = paste0("mean_", ct$A_lab), mean_B = paste0("mean_", ct$B_lab),
             n_A = paste0("n_", ct$A_lab), n_B = paste0("n_", ct$B_lab),
             n_A_fullcells = paste0("n_", ct$A_lab, "_fullcells"),
             n_B_fullcells = paste0("n_", ct$B_lab, "_fullcells"))
    names(d)[match(names(ren), names(d))] <- unname(ren)

    d <- d[order(d$padj, -abs(d$log2FoldChange)), , drop = FALSE]
    rownames(d) <- NULL
    tables[[key]] <- d

    nup <- sum(d$direction == ct$up); ndn <- sum(d$direction == ct$down)
    manifest[[length(manifest) + 1]] <- data.frame(
      key = key, contrast = ck, cluster = cl, stratum = st,
      is_primary = st == ct$primary,
      arm_A = ct$A, arm_B = ct$B, n_A = length(iA), n_B = length(iB),
      n_A_fullcells = length(jA), n_B_fullcells = length(jB),
      n_genes_tested = nrow(d), n_up = nup, n_down = ndn,
      up_label = ct$up, down_label = ct$down, stringsAsFactors = FALSE)
    cat(sprintf("  %-32s A=%-4d B=%-4d  genes=%d  up=%-5d down=%-5d%s\n",
                key, length(iA), length(iB), nrow(d), nup, ndn,
                if (st == ct$primary) "  <- primary" else ""))
    write.csv(d, file.path(OUT, "csv", paste0(key, ".csv")), row.names = FALSE, na = "")
  }
}

manifest <- do.call(rbind, manifest)
# If a cluster's primary stratum was skipped (CM4 has no G1 cells at all, so the
# G1-matched read of part 2 does not exist for it), promote the surviving stratum
# so the cluster still gets figures -- flagged, not silently substituted.
manifest$primary_fallback <- FALSE
for (ck in unique(manifest$contrast)) for (cl in unique(manifest$cluster[manifest$contrast == ck])) {
  sel <- manifest$contrast == ck & manifest$cluster == cl
  if (!any(manifest$is_primary[sel]) && any(sel)) {
    j <- which(sel)[1]
    manifest$is_primary[j] <- TRUE; manifest$primary_fallback[j] <- TRUE
    cat(sprintf("  NOTE %s/%s: primary stratum unavailable, promoting '%s' for figures\n",
                ck, cl, manifest$stratum[j]))
  }
}
skipped  <- if (length(skipped)) do.call(rbind, skipped) else
            data.frame(contrast = character(), cluster = character(), stratum = character(),
                       n_A = integer(), n_B = integer(), reason = character())

saveRDS(list(tables = tables, manifest = manifest, skipped = skipped,
             sets = SETS, confound = CONF,
             params = list(sig_padj = SIG_PADJ, sig_lfc = SIG_LFC, min_cells = MIN_CELLS,
                           res_col = RES_COL, clusters = CLUSTERS,
                           contrasts = CONTRASTS,
                           matrix = sprintf("broad app$deg_expr (%d genes x %d cells)", nrow(Mb), ncol(Mb)),
                           crosscheck_matrix = sprintf("curated app$expr (%d genes x %d cells)", nrow(Mc), ncol(Mc)),
                           built = as.character(Sys.time()))),
        file.path(OUT, "de_tables.rds"))
write.csv(manifest, file.path(OUT, "csv", "_manifest.csv"), row.names = FALSE)
write.csv(skipped,  file.path(OUT, "csv", "_skipped.csv"),  row.names = FALSE)

cat("\n== manifest ==\n"); print(manifest[, c("key","n_A","n_B","n_genes_tested","n_up","n_down")], row.names = FALSE)
if (nrow(skipped)) { cat("\n== skipped ==\n"); print(skipped, row.names = FALSE) }
cat("\nDONE 01_de.R\n")
