#!/usr/bin/env Rscript
# Where did every gene set come from? -- an audit, generated FROM the code.
#
# The honest headline, which this script exists to make visible rather than bury: most of
# the gene sets in this project are HAND-CURATED with no recorded source. They are not
# wrong -- Myh6/Myh7, Nppa, the Mcm family and so on are standard cardiac and cell-cycle
# markers that any reviewer would recognise -- but the SELECTION is unauditable. Nothing
# records why Ckm is in the mature panel and Myom1 is not, or what version of anyone's
# list it came from. That is a real limitation and it belongs on the page next to the
# scores those panels produce.
#
# A minority of sets DO have real provenance, because they are fetched from a package at
# runtime rather than typed into a script. Those are marked accordingly.
#
# Gene lists are EXTRACTED from the live definitions (parse + eval of the literal), never
# retyped here -- a registry that is retyped is a registry that goes stale. Only the
# source annotation is maintained by hand. The script also cross-checks sets that are
# defined in more than one place and reports any that have drifted apart.
#
#   Writes results/tables/gene_set_provenance.csv
#
# Runs in e2f-seurat-full:latest (needs Seurat only for cc.genes.updated.2019).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))

REPO <- normalizePath(file.path(OUR_ROOT, ".."))
F_COMMON <- file.path(OUR_ROOT, "_common.R")
F_SIG    <- file.path(REPO, "shiny_app", "build_signature_scores.R")
F_APP    <- file.path(REPO, "shiny_app", "app.R")
F_SUB    <- file.path(OUR_ROOT, "05_analyses", "cm_subtypes.R")
F_CYC    <- file.path(OUR_ROOT, "05_analyses", "cm_cycling_investigate.R")
F_SEX    <- file.path(OUR_ROOT, "05_analyses", "sex_check.R")

# Pull the right-hand side of `name <- <literal>` out of a file and evaluate just that.
# Parsing rather than sourcing keeps this free of the side effects (directory creation,
# library loading) that sourcing _common.R or app.R would drag in.
extract <- function(file, name) {
  if (!file.exists(file)) return(NULL)
  for (e in parse(file, keep.source = FALSE)) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1]]) %in% c("<-", "=") &&
        identical(as.character(e[[2]])[1], name))
      return(tryCatch(eval(e[[3]], envir = baseenv()), error = function(err) NULL))
  }
  NULL
}
line_of <- function(file, pattern) {
  if (!file.exists(file)) return(NA_integer_)
  i <- grep(pattern, readLines(file, warn = FALSE))
  if (length(i)) i[1] else NA_integer_
}
rel <- function(p) sub(paste0("^", REPO, "/"), "", p)

## --- the source annotation: the ONLY hand-maintained part -------------------
# source_type is deliberately blunt. "hand-curated" means: typed into a script in this
# project, with no citation, no database, and no version. Saying "canonical" instead --
# which is what the code comments currently say -- makes an unsourced list sound sourced.
HAND <- "hand-curated (no source recorded)"
EXT  <- "external database"
DERV <- "derived from another set in this table"
NOTE_HAND <- paste("Typed into this project with no citation. The individual genes are",
                   "standard markers, but the selection is not traceable to a published",
                   "list and cannot be version-checked.")

sets <- list()
add <- function(name, genes, file, pat, source_type, source, note = "") {
  if (is.null(genes) || !length(genes)) return(invisible(NULL))
  sets[[length(sets) + 1L]] <<- data.frame(
    set = name, n_genes = length(genes), source_type = source_type, source = source,
    defined_in = sprintf("%s:%s", rel(file), line_of(file, pat)),
    note = note, genes = paste(sort(unique(genes)), collapse = ", "),
    stringsAsFactors = FALSE)
}

## --- sets WITH real provenance ----------------------------------------------
cc <- tryCatch(Seurat::cc.genes.updated.2019, error = function(e) NULL)
if (!is.null(cc)) {
  add("Cell cycle S (Seurat)", cc$s.genes, F_COMMON, "cc_lists <- function", EXT,
      "Seurat::cc.genes.updated.2019$s.genes",
      paste("Tirosh et al. 2016 regev-lab cell-cycle genes, 2019 update, shipped in",
            "Seurat. Human symbols; mapped to mouse by babelgene::orthologs via",
            "to_mouse() -- the mapping is itself a source of loss, see coverage."))
  add("Cell cycle G2M (Seurat)", cc$g2m.genes, F_COMMON, "cc_lists <- function", EXT,
      "Seurat::cc.genes.updated.2019$g2m.genes",
      "Same source and same ortholog caveat as the S list.")
}
# Runtime-fetched collections: recorded as rows even though the genes are not enumerated
# here, because their provenance is the point and their content is version-dependent.
runtime <- data.frame(
  set = c("MSigDB Hallmark", "MSigDB KEGG_LEGACY", "MSigDB GTRD TFT", "GO Biological Process",
          "CellChatDB (mouse)"),
  n_genes = NA_integer_, source_type = EXT,
  source = c("msigdbr(species='Mus musculus', collection='H')",
             "msigdbr(collection='C2', subcollection='CP:KEGG_LEGACY')",
             "msigdbr(collection='C3', subcollection='TFT:GTRD')",
             "clusterProfiler::enrichGO + org.Mm.eg.db",
             "CellChat::CellChatDB.mouse"),
  defined_in = c(sprintf("our_analysis/05_analyses/pathway_msigdb.R:%s", line_of(file.path(OUR_ROOT,"05_analyses","pathway_msigdb.R"), "collection = \"H\"")),
                 sprintf("our_analysis/05_analyses/pathway_msigdb.R:%s", line_of(file.path(OUR_ROOT,"05_analyses","pathway_msigdb.R"), "KEGG_LEGACY")),
                 sprintf("our_analysis/05_analyses/pathway_msigdb.R:%s", line_of(file.path(OUR_ROOT,"05_analyses","pathway_msigdb.R"), "TFT:GTRD")),
                 "shiny_app/build_subcluster_enrichment.R, build_fourgroup_enrichment.R",
                 sprintf("our_analysis/05_analyses/cellchat.R:%s", line_of(file.path(OUR_ROOT,"05_analyses","cellchat.R"), "CellChatDB.mouse"))),
  note = c(rep(paste("Fetched at runtime, so content follows the installed msigdbr /",
                     "MSigDB release. No version is pinned anywhere in this repo --",
                     "see docs/11-confounds-repro.qmd on the missing environment record."), 3),
           "Fetched at runtime from the installed org.Mm.eg.db; no version pinned.",
           "Ships with the CellChat package; no version pinned."),
  genes = "(fetched at runtime -- not enumerated here)", stringsAsFactors = FALSE)

## --- hand-curated sets, extracted from the code -----------------------------
add("CM selection markers", extract(F_COMMON, "CM_MARKERS"), F_COMMON, "^CM_MARKERS",
    HAND, "our_analysis/_common.R", NOTE_HAND)
ct <- extract(F_COMMON, "CELLTYPE_MARKERS")
for (nm in names(ct))
  add(sprintf("Cell type: %s", nm), ct[[nm]], F_COMMON, "^CELLTYPE_MARKERS", HAND,
      "our_analysis/_common.R", paste(NOTE_HAND, "Drives the celltype annotation in",
      "annotate.R by marker-module argmax, so these choices set the cell-type calls."))
add("E2F targets", extract(F_COMMON, "E2F_TARGETS"), F_COMMON, "^E2F_TARGETS", HAND,
    "our_analysis/_common.R",
    paste(NOTE_HAND, "The code comment calls it a 'curated canonical set'; that phrasing",
          "should not be read as a citation. Note it is heavily cell-cycle weighted, so a",
          "high E2F-target score in a cycling cell is close to tautological."))
add("CM mature program", extract(F_COMMON, "CM_MATURE"), F_COMMON, "^CM_MATURE", HAND,
    "our_analysis/_common.R", NOTE_HAND)
add("CM immature program", extract(F_COMMON, "CM_IMMATURE"), F_COMMON, "^CM_IMMATURE", HAND,
    "our_analysis/_common.R", NOTE_HAND)

sig <- extract(F_SIG, "SETS")
sig_notes <- c(
  prolif = "", cytokinesis = "", ccexit = "",
  mat_mature = "Feeds sig_maturation, the maturation score shown in the app.",
  mat_immature = "Feeds sig_maturation. Contains Ccnd1/Mki67/Top2a, so the score is partly a cell-cycle score.",
  mat_immature_nocc = "mat_immature with the three cell-cycle genes removed, to break the circularity above.",
  glycolysis = "Feeds sig_metabolic (FAO/OXPHOS minus glycolysis).",
  faox = "Feeds sig_metabolic.")
for (nm in names(sig))
  add(sprintf("Signature: %s", nm), sig[[nm]], F_SIG, "^SETS <- list",
      if (identical(nm, "mat_immature_nocc")) DERV else HAND,
      "shiny_app/build_signature_scores.R",
      trimws(paste(if (identical(nm, "mat_immature_nocc")) "" else NOTE_HAND,
                   sig_notes[[nm]])))

sub <- extract(F_SUB, "CM_SUBTYPES")
for (nm in names(sub))
  add(sprintf("CM subtype: %s", nm), sub[[nm]], F_SUB, "^CM_SUBTYPES", HAND,
      "our_analysis/05_analyses/cm_subtypes.R",
      paste("The script header attributes these to the prior pipeline's _env.R",
            "CARDIAC_MARKERS, but that file is NOT in this repository -- the pointer",
            "cannot be followed, so this is effectively unsourced."))
add("Cycling markers", extract(F_CYC, "CYCLE_MARKERS"), F_CYC, "^CYCLE_MARKERS", HAND,
    "our_analysis/05_analyses/cm_cycling_investigate.R", NOTE_HAND)
sx <- extract(F_SEX, "SEX")
if (!is.null(sx))
  add("Sex markers", unlist(sx), F_SEX, "^SEX <- list", HAND,
      "our_analysis/05_analyses/sex_check.R",
      paste("Xist plus four Y-linked genes. Conventional and effectively unambiguous for",
            "sexing, unlike the phenotype panels -- listed for completeness."))

app_sets <- extract(F_APP, ".gene_sets_raw")
for (nm in names(app_sets))
  add(sprintf("App dropdown: %s", nm), app_sets[[nm]], F_APP, "^\\.gene_sets_raw", HAND,
      "shiny_app/app.R", paste(NOTE_HAND, "Used only to filter the gene dropdowns; it does",
      "not feed any score."))

reg <- rbind(do.call(rbind, sets), runtime)
reg <- reg[order(reg$source_type != EXT, reg$set), ]

## --- drift check: the same panel defined twice, differently ------------------
# The reason this check is here: the app's "CM maturation" dropdown and the sig_maturation
# score it plots were found to use different gene lists under the same name.
cat("\n=== duplicate-name drift check ===\n")
norm <- function(g) sort(unique(trimws(strsplit(g, ",")[[1]])))
pairs <- list(
  c("CM mature program",              "Signature: mat_mature"),
  c("CM immature program",            "Signature: mat_immature"),
  c("App dropdown: E2F targets",      "E2F targets"),
  c("App dropdown: Cardiomyocyte",    "Cell type: Cardiomyocyte"))
drift <- 0L
for (p in pairs) {
  a <- reg$genes[reg$set == p[1]]; b <- reg$genes[reg$set == p[2]]
  if (!length(a) || !length(b)) next
  ga <- norm(a); gb <- norm(b)
  if (identical(ga, gb)) { cat(sprintf("  ok       %s == %s\n", p[1], p[2])); next }
  drift <- drift + 1L
  cat(sprintf("  DRIFTED  %s (%d) vs %s (%d)\n    only in first : %s\n    only in second: %s\n",
              p[1], length(ga), p[2], length(gb),
              paste(setdiff(ga, gb), collapse = ", "), paste(setdiff(gb, ga), collapse = ", ")))
}
cat(sprintf("  -> %d of %d checked pairs have drifted\n", drift, length(pairs)))

out <- file.path(OUTTAB, "gene_set_provenance.csv")
write.csv(reg, out, row.names = FALSE)
cat(sprintf("\n%d sets: %d from an external database, %d hand-curated, %d derived\n",
            nrow(reg), sum(reg$source_type == EXT), sum(reg$source_type == HAND),
            sum(reg$source_type == DERV)))
cat("wrote: ", out, "\n=== DONE gene_set_provenance ===\n", sep = "")
