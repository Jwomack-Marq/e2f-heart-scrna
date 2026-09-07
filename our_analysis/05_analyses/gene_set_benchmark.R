#!/usr/bin/env Rscript
# Are our hand-curated panels defensible? -- benchmark them against MSigDB Hallmark.
#
# gene_set_provenance.R establishes that 38 of our 46 gene sets were typed into a script
# with no citation. That answers "where did these come from" (nowhere recorded) but not
# "are they right". There is no ground truth for a marker panel, so "right" is not
# directly askable. What IS askable, and what this does: for each hand panel that has a
# published, versioned counterpart, how much do they agree, and which of our genes are
# choices nobody else made?
#
# Read the output this way:
#   - high overlap  -> our panel is a subset of an established set; defensible, and the
#                      honest fix is to cite the reference rather than keep the copy.
#   - our_only      -> genes we added that the reference does not carry. Each one is a
#                      judgement call with no record behind it. Not necessarily wrong --
#                      Hallmark is human-centric and deliberately broad -- but it is
#                      exactly the part a reviewer can question and we cannot defend.
#   - no counterpart-> the panel is ours alone. CM maturation is the important case: there
#                      is no Hallmark set for it, so the maturation score in this app rests
#                      on a list with no external anchor at all.
#
# Hallmark is human symbols; mapped to mouse with the same to_mouse() the rest of the
# pipeline uses, so the mapping loss here is the same loss the cell-cycle scoring carries.
#
#   Writes results/tables/gene_set_benchmark.csv
# Runs in e2f-seurat-full:latest. msigdbr's Hallmark cache is warmed into the image, so
# this needs no network.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages(library(msigdbr)))

REPO  <- normalizePath(file.path(OUR_ROOT, ".."))
F_SIG <- file.path(REPO, "shiny_app", "build_signature_scores.R")

extract <- function(file, name) {
  for (e in parse(file, keep.source = FALSE))
    if (is.call(e) && length(e) >= 3L && as.character(e[[1]]) %in% c("<-", "=") &&
        identical(as.character(e[[2]])[1], name))
      return(tryCatch(eval(e[[3]], envir = baseenv()), error = function(err) NULL))
  NULL
}

h <- msigdbr(species = "Mus musculus", collection = "H")
gs_col <- intersect(c("gs_name"), names(h))[1]
sym_col <- intersect(c("gene_symbol", "mouse_symbol", "symbol"), names(h))[1]
stopifnot(!is.na(gs_col), !is.na(sym_col))
hall <- split(h[[sym_col]], h[[gs_col]])
hall <- lapply(hall, function(g) sort(unique(g)))
cat(sprintf("Hallmark: %d sets, mouse symbols via msigdbr\n\n", length(hall)))

SIG <- extract(F_SIG, "SETS")
ours <- list(
  "glycolysis (sig_metabolic neg)"  = SIG$glycolysis,
  "faox (sig_metabolic pos)"        = SIG$faox,
  "prolif (sig_ploidy pos)"         = SIG$prolif,
  "cytokinesis (sig_ploidy neg)"    = SIG$cytokinesis,
  "ccexit (sig_ccexit)"             = SIG$ccexit,
  "mat_mature (sig_maturation pos)" = SIG$mat_mature,
  "mat_immature (sig_maturation neg)" = SIG$mat_immature,
  "E2F_TARGETS"                     = extract(file.path(OUR_ROOT, "_common.R"), "E2F_TARGETS"))
# Counterparts chosen because the panel claims to measure the same thing, not because they
# maximise overlap. NA = we could find no Hallmark set that claims to measure this.
REF <- list(
  "glycolysis (sig_metabolic neg)"  = "HALLMARK_GLYCOLYSIS",
  "faox (sig_metabolic pos)"        = "HALLMARK_FATTY_ACID_METABOLISM",
  "prolif (sig_ploidy pos)"         = "HALLMARK_G2M_CHECKPOINT",
  "cytokinesis (sig_ploidy neg)"    = "HALLMARK_G2M_CHECKPOINT",
  "ccexit (sig_ccexit)"             = "HALLMARK_P53_PATHWAY",
  "mat_mature (sig_maturation pos)" = NA_character_,
  "mat_immature (sig_maturation neg)" = NA_character_,
  "E2F_TARGETS"                     = "HALLMARK_E2F_TARGETS")
EXTRA <- list("faox (sig_metabolic pos)" = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
# KEGG as a SECOND reference for the metabolic panels, and it is the fairer one. Hallmark
# sets are derived from co-expression, not pathway membership, so HALLMARK_GLYCOLYSIS
# genuinely omits Gapdh, Hk1, Pfkm, Pfkl and Slc2a1 -- core glycolytic enzymes. Scoring our
# panel only against Hallmark would therefore penalise it for being MORE canonical, not
# less. KEGG is curated by pathway membership and is the right yardstick here.
KEGG_REF <- list("glycolysis (sig_metabolic neg)" = "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                 "faox (sig_metabolic pos)"       = "KEGG_FATTY_ACID_METABOLISM")

rows <- list()
one <- function(label, g, refname, tag = "primary") {
  if (is.na(refname)) {
    rows[[length(rows) + 1L]] <<- data.frame(
      our_set = label, n_ours = length(g), reference = "(none found)", n_ref = NA_integer_,
      n_overlap = NA_integer_, pct_ours_covered = NA_real_, jaccard = NA_real_,
      match = tag, our_only = paste(sort(g), collapse = ", "), ref_only = "",
      stringsAsFactors = FALSE)
    cat(sprintf("  %-34s %2d genes   NO Hallmark counterpart -- unanchored\n", label, length(g)))
    return(invisible())
  }
  r <- hall[[refname]]
  ov <- intersect(g, r)
  rows[[length(rows) + 1L]] <<- data.frame(
    our_set = label, n_ours = length(g), reference = refname, n_ref = length(r),
    n_overlap = length(ov), pct_ours_covered = 100 * length(ov) / length(g),
    jaccard = length(ov) / length(union(g, r)), match = tag,
    our_only = paste(sort(setdiff(g, r)), collapse = ", "),
    ref_only = paste(head(sort(setdiff(r, g)), 25), collapse = ", "), stringsAsFactors = FALSE)
  cat(sprintf("  %-34s %2d genes   %-38s %2d/%2d ours in ref (%.0f%%)\n",
              label, length(g), sub("HALLMARK_", "", refname), length(ov), length(g),
              100 * length(ov) / length(g)))
  if (length(setdiff(g, r)))
    cat(sprintf("       ours only: %s\n", paste(sort(setdiff(g, r)), collapse = ", ")))
}
cat("=== our hand panels vs MSigDB Hallmark ===\n")
for (nm in names(ours)) one(nm, ours[[nm]], REF[[nm]])
cat("\n=== secondary counterparts ===\n")
for (nm in names(EXTRA)) one(nm, ours[[nm]], EXTRA[[nm]], tag = "secondary")

kg <- tryCatch(msigdbr(species = "Mus musculus", collection = "C2",
                       subcollection = "CP:KEGG_LEGACY"), error = function(e) NULL)
if (!is.null(kg)) {
  kegg <- lapply(split(kg[[sym_col]], kg[[gs_col]]), function(g) sort(unique(g)))
  hall <- c(hall, kegg)                      # one() looks sets up in `hall`
  cat("\n=== metabolic panels vs KEGG (membership-curated; the fairer reference) ===\n")
  for (nm in names(KEGG_REF)) if (KEGG_REF[[nm]] %in% names(hall))
    one(nm, ours[[nm]], KEGG_REF[[nm]], tag = "kegg")
} else cat("\n(KEGG_LEGACY unavailable -- skipped)\n")

out <- do.call(rbind, rows)
write.csv(out, file.path(OUTTAB, "gene_set_benchmark.csv"), row.names = FALSE)
cov <- out$pct_ours_covered[out$match == "primary" & !is.na(out$pct_ours_covered)]
cat(sprintf("\n%d panels benchmarked; median coverage of our genes by the reference: %.0f%%\n",
            length(cov), median(cov)))
cat(sprintf("%d panel(s) have NO published counterpart at all\n",
            sum(out$reference == "(none found)")))
cat("wrote: ", file.path(OUTTAB, "gene_set_benchmark.csv"), "\n", sep = "")
cat("=== DONE gene_set_benchmark ===\n")
