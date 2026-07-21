# build_communication.R
# ---------------------------------------------------------------------------
# Curated ligand-receptor (cell-cell communication) scoring, focused on the
# E2F7/8 -> Vegfa angiogenic axis, written back into app_data.rds as app$commun.
#
# This is a DESCRIPTIVE, curated-pair score (NATMI/connectome style: mean ligand
# in the sender x mean receptor in the receiver, on log-norm expression) — NOT a
# permutation-tested CellChat/LIANA run (that would need the source Seurat
# objects). No significance; n=1 per condition, so KO-WT deltas are
# hypothesis-generating only.
#
# Run LOCALLY, in a real R session (NOT the bash sandbox):
#     & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" shiny_app/build_communication.R
#
# Output app$commun = list(
#   scores  : data.frame(pathway, ligand, receptor, sender, receiver, timepoint,
#                         WT, KO, delta, matrix)   # delta = KO - WT
#   pairs   : the curated pair table actually used (both genes found)
#   dropped : pairs skipped because a gene was absent from both matrices)
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))

if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

# ---- curated ligand -> receptor pairs (mouse symbols) ---------------------
PAIRS <- rbind(
  data.frame(pathway = "VEGF",        ligand = "Vegfa",  receptor = "Kdr"),
  data.frame(pathway = "VEGF",        ligand = "Vegfa",  receptor = "Flt1"),
  data.frame(pathway = "VEGF",        ligand = "Vegfa",  receptor = "Nrp1"),
  data.frame(pathway = "VEGF",        ligand = "Vegfb",  receptor = "Flt1"),
  data.frame(pathway = "VEGF",        ligand = "Pgf",    receptor = "Flt1"),
  data.frame(pathway = "VEGF",        ligand = "Vegfc",  receptor = "Flt4"),
  data.frame(pathway = "Angiopoietin",ligand = "Angpt1", receptor = "Tek"),
  data.frame(pathway = "Angiopoietin",ligand = "Angpt2", receptor = "Tek"),
  data.frame(pathway = "Notch",       ligand = "Dll4",   receptor = "Notch1"),
  data.frame(pathway = "Notch",       ligand = "Dll4",   receptor = "Notch4"),
  data.frame(pathway = "Notch",       ligand = "Jag1",   receptor = "Notch1"),
  data.frame(pathway = "PDGF",        ligand = "Pdgfb",  receptor = "Pdgfrb"),
  data.frame(pathway = "PDGF",        ligand = "Pdgfa",  receptor = "Pdgfra"),
  data.frame(pathway = "Apelin",      ligand = "Apln",   receptor = "Aplnr"),
  data.frame(pathway = "Ephrin",      ligand = "Efnb2",  receptor = "Ephb4"),
  data.frame(pathway = "BMP",         ligand = "Bmp10",  receptor = "Bmpr2"),
  data.frame(pathway = "Neuregulin",  ligand = "Nrg1",   receptor = "Erbb4"),
  data.frame(pathway = "Neuregulin",  ligand = "Nrg1",   receptor = "Erbb2"),
  data.frame(pathway = "IGF",         ligand = "Igf1",   receptor = "Igf1r"),
  data.frame(pathway = "Natriuretic", ligand = "Nppa",   receptor = "Npr1"),
  stringsAsFactors = FALSE)

cat("Loading app_data.rds ...\n")
app  <- readRDS("app_data.rds")
EXPRc <- app$expr
EXPRb <- app$deg_expr
metaC <- app$meta
metaB <- app$deg_meta        # may be NULL
stopifnot(!is.null(EXPRc), !is.null(metaC),
          all(c("celltype","genotype","timepoint","cell") %in% names(metaC)))

# metadata aligned to a matrix's columns (matched by cell id)
aligned <- function(M, DF) {
  if (is.null(M) || is.null(DF)) return(NULL)
  DF[match(colnames(M), DF$cell), , drop = FALSE]
}
DFc <- aligned(EXPRc, metaC)
# only use the broad matrix if its metadata carries the grouping columns
DFb <- if (!is.null(EXPRb) && !is.null(metaB) &&
           all(c("celltype","genotype","timepoint","cell") %in% names(metaB)))
         aligned(EXPRb, metaB) else NULL

# genes x group mean matrix; group = celltype|genotype|timepoint
group_means <- function(M, DF, genes) {
  genes <- intersect(genes, rownames(M))
  if (!length(genes)) return(NULL)
  grp <- paste(DF$celltype, DF$genotype, DF$timepoint, sep = "|")
  ix_by <- split(seq_len(ncol(M)), grp)
  sub <- M[genes, , drop = FALSE]
  out <- vapply(ix_by, function(ix) Matrix::rowMeans(sub[, ix, drop = FALSE]), numeric(length(genes)))
  if (length(genes) == 1) { out <- matrix(out, nrow = 1); }
  rownames(out) <- genes; out
}

all_genes <- unique(c(PAIRS$ligand, PAIRS$receptor))
GMc <- group_means(EXPRc, DFc, all_genes)
GMb <- if (!is.null(DFb)) group_means(EXPRb, DFb, all_genes) else NULL

has_gene <- function(GM, g) !is.null(GM) && g %in% rownames(GM)
celltypes <- sort(unique(as.character(metaC$celltype)))
timepoints <- sort(unique(as.character(metaC$timepoint)))
geno <- c("WT","KO")

cell_score <- function(GM, gene, ct, gt, tp) {
  key <- paste(ct, gt, tp, sep = "|")
  if (!has_gene(GM, gene) || !key %in% colnames(GM)) return(NA_real_)
  GM[gene, key]
}

cat("\n== scoring curated pairs ==\n")
rows <- list(); used_idx <- integer(0); dropped_idx <- integer(0)
for (i in seq_len(nrow(PAIRS))) {
  lg <- PAIRS$ligand[i]; rc <- PAIRS$receptor[i]
  # prefer the curated panel; fall back to the broad matrix only if BOTH genes live there
  GM <- if (has_gene(GMc, lg) && has_gene(GMc, rc)) GMc
        else if (has_gene(GMb, lg) && has_gene(GMb, rc)) GMb else NULL
  matname <- if (identical(GM, GMc)) "curated" else "broad"
  if (is.null(GM)) { dropped_idx <- c(dropped_idx, i)
                     cat(sprintf("  %-6s -> %-7s  DROP (ligand/receptor absent)\n", lg, rc)); next }
  used_idx <- c(used_idx, i)
  for (s in celltypes) for (r in celltypes) for (tp in timepoints) {
    wt <- cell_score(GM, lg, s, "WT", tp) * cell_score(GM, rc, r, "WT", tp)
    ko <- cell_score(GM, lg, s, "KO", tp) * cell_score(GM, rc, r, "KO", tp)
    if (is.na(wt) && is.na(ko)) next
    rows[[length(rows) + 1]] <- data.frame(
      pathway = PAIRS$pathway[i], ligand = lg, receptor = rc,
      sender = s, receiver = r, timepoint = tp,
      WT = wt, KO = ko, delta = ko - wt, matrix = matname, stringsAsFactors = FALSE)
  }
  cat(sprintf("  %-6s -> %-7s  %s\n", lg, rc, matname))
}

app$commun <- list(
  scores  = if (length(rows)) do.call(rbind, rows) else NULL,
  pairs   = PAIRS[used_idx, , drop = FALSE],
  dropped = if (length(dropped_idx)) PAIRS[dropped_idx, , drop = FALSE] else NULL)

cat(sprintf("\n== summary ==\n  scored rows: %s | pairs used: %d | dropped: %d\n",
            if (is.null(app$commun$scores)) 0 else nrow(app$commun$scores),
            length(used_idx), length(dropped_idx)))

cat("\nBacking up -> app_data.pre_commun.bak.rds\n")
file.copy("app_data.rds", "app_data.pre_commun.bak.rds", overwrite = TRUE)
cat("Saving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
