# build_gene_provenance.R
# ---------------------------------------------------------------------------
# Carry the gene-set provenance audit into app_data.rds as app$genesets, for the
# "Gene sets & sources" tab.
#
# The audit is our_analysis/05_analyses/gene_set_provenance.R, which extracts every gene
# list from the live code (rather than retyping it) and annotates where it came from.
# This script only reads the CSV it writes, so it needs no Seurat.
#
# It also re-runs the duplicate-name drift check here, from the gene lists in the CSV, so
# the app cannot show a clean bill of health that the underlying data does not support.
#
#   Rscript shiny_app/build_gene_provenance.R [--tables=<dir>]
# ---------------------------------------------------------------------------

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

TABLES <- argval("--tables", "../our_analysis/results/tables")
f <- file.path(TABLES, "gene_set_provenance.csv")
if (!file.exists(f))
  stop("provenance audit not found at ", f,
       "\n  run our_analysis/05_analyses/gene_set_provenance.R first")

cat("Reading provenance audit ...\n")
reg <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("set","n_genes","source_type","source","defined_in","note","genes") %in% names(reg)))

EXT  <- "external database"
HAND <- "hand-curated (no source recorded)"
n_ext  <- sum(reg$source_type == EXT)
n_hand <- sum(reg$source_type == HAND)

# Drift: the same biological panel defined in two places under two names. Recomputed from
# the CSV rather than trusted from the audit's console output, so the badge the tab shows
# is always backed by the data actually in the bundle.
norm <- function(g) sort(unique(trimws(strsplit(g, ",\\s*")[[1]])))
PAIRS <- list(c("CM mature program",   "Signature: mat_mature"),
              c("CM immature program", "Signature: mat_immature"),
              c("App dropdown: E2F targets",   "E2F targets"),
              c("App dropdown: Cardiomyocyte", "Cell type: Cardiomyocyte"))
drift <- do.call(rbind, lapply(PAIRS, function(p) {
  a <- reg$genes[reg$set == p[1]]; b <- reg$genes[reg$set == p[2]]
  if (!length(a) || !length(b)) return(NULL)
  ga <- norm(a); gb <- norm(b)
  data.frame(set_a = p[1], n_a = length(ga), set_b = p[2], n_b = length(gb),
             drifted = !identical(ga, gb),
             only_in_a = paste(setdiff(ga, gb), collapse = ", "),
             only_in_b = paste(setdiff(gb, ga), collapse = ", "),
             stringsAsFactors = FALSE) }))
n_drift <- if (is.null(drift)) 0L else sum(drift$drifted)
cat(sprintf("  %d sets | %d external, %d hand-curated | %d of %d checked pairs drifted\n",
            nrow(reg), n_ext, n_hand, n_drift, if (is.null(drift)) 0L else nrow(drift)))

app <- readRDS("app_data.rds")
cat("Backing up -> app_data.pre_genesets.bak.rds\n")
saveRDS(app, "app_data.pre_genesets.bak.rds", compress = "gzip")

# Benchmark: how far our hand panels agree with a published, versioned counterpart.
# Optional -- the tab degrades to registry-only if the benchmark has not been run.
fb <- file.path(TABLES, "gene_set_benchmark.csv")
bench <- if (file.exists(fb)) read.csv(fb, stringsAsFactors = FALSE, check.names = FALSE) else NULL
if (!is.null(bench)) {
  prim <- bench[bench$match == "primary" & !is.na(bench$pct_ours_covered), ]
  cat(sprintf("  benchmark: %d panels, median %.0f%% of our genes covered; %d unanchored\n",
              nrow(prim), median(prim$pct_ours_covered), sum(bench$reference == "(none found)")))
}

# Verified references. Every entry here was checked against the actual record -- these are
# NOT suggested-from-memory citations, which in this domain is how a plausible-looking but
# non-existent paper ends up in a methods section.
REFS <- data.frame(
  topic = c("CM maturation", "CM cell-cycle exit", "Cell-cycle S / G2M",
            "Metabolic panels", "Metabolic panels"),
  reference = c(
    "Uosaki H et al. Transcriptional Landscape of Cardiomyocyte Maturation. Cell Reports 2015;13(8):1705-1716.",
    "Mahmoud AI et al. Meis1 regulates postnatal cardiomyocyte cell cycle arrest. Nature 2013;497(7448):249-253.",
    "Tirosh I et al. 2016, as shipped in Seurat::cc.genes.updated.2019 (2019 update).",
    "MSigDB Hallmark, via msigdbr.",
    "MSigDB KEGG_LEGACY (C2:CP), via msigdbr."),
  link = c("https://doi.org/10.1016/j.celrep.2015.10.032",
           "https://pmc.ncbi.nlm.nih.gov/articles/PMC4159712/",
           "https://satijalab.org/seurat/",
           "https://www.gsea-msigdb.org/gsea/msigdb/",
           "https://www.gsea-msigdb.org/gsea/msigdb/"),
  relevance = c(
    "The canonical maturation transcriptional atlas (>200 arrays, embryonic to adult). We do NOT use its gene list -- our mature/immature panels were written independently and have no external anchor. This is the reference to reconcile them against.",
    "Establishes the postnatal proliferative window closing at P7 and Meis1 driving arrest via p15/p16/p21 -- the same biology, and the same timepoints, as this study. Supports Meis1 and the Cdkn genes in the ccexit panel.",
    "Source of the S and G2M lists. Human symbols mapped to mouse by babelgene, which loses genes.",
    "Reference for prolif, ccexit and E2F targets in the benchmark. Hallmark sets are co-expression-derived, not pathway membership.",
    "Reference for the glycolysis and FAO panels. Membership-curated, and the fairer yardstick for a metabolic panel than Hallmark."),
  stringsAsFactors = FALSE)

app$genesets <- list(
  benchmark = bench, refs = REFS,
  registry = reg, drift = drift, n_ext = n_ext, n_hand = n_hand, n_drift = n_drift,
  headline = paste0(
    "Of the ", nrow(reg), " gene sets this project uses, ", n_ext, " are fetched from an ",
    "external database and ", n_hand, " are hand-curated: typed into a script here with no ",
    "citation, no version, and no record of why a gene is in or out. The genes themselves ",
    "are standard markers, but the selections cannot be traced to a published list or ",
    "version-checked, and that applies to the maturation, metabolic and E2F-target panels ",
    "behind the scores in this app."),
  bench_note = paste(
    "There is no ground truth for a marker panel, so \"is it right\" is not directly",
    "askable. What is askable: for each panel with a published counterpart, how much do",
    "they agree, and which genes are choices nobody else made? Read a low overlap in both",
    "directions -- HALLMARK_GLYCOLYSIS omits Gapdh, Hk1, Pfkm, Pfkl and Slc2a1, so scoring",
    "our panel against it would penalise it for being MORE canonical. Against KEGG, which",
    "curates by pathway membership, the same panel is 10/11."),
  caveats = c(
    "E2F targets is heavily cell-cycle weighted, so a high E2F-target score in a cycling cell is close to tautological.",
    "sig_maturation's immature side contains Ccnd1, Mki67 and Top2a, so it is partly a cell-cycle score; sig_maturation_nocc drops those and is the one to use when the argument involves cycling.",
    "The CM subtype panels are attributed to a prior pipeline's _env.R, but that file is not in this repository, so the pointer cannot be followed.",
    "MSigDB, GO and CellChatDB content follows whatever package version is installed; no version is pinned anywhere in this repo.",
    "The Seurat cell-cycle lists are human symbols mapped to mouse by babelgene, and that mapping loses genes."))

cat("\nSaving app_data.rds (gzip) ...\n")
saveRDS(app, "app_data.rds", compress = "gzip")
cat("Done.\n")
