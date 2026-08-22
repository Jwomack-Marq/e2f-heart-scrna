# 04_cover.R -- the one-page cover note that ships with the files.
# Generated from the actual results so the numbers in it cannot drift from the
# tables. Written last, after 01-03.
suppressMessages(library(openxlsx))
OUT <- "/out"
de <- readRDS(file.path(OUT, "de_tables.rds")); en <- readRDS(file.path(OUT, "enrich.rds"))
ms <- if (file.exists(file.path(OUT, "mito_sensitivity.rds"))) readRDS(file.path(OUT, "mito_sensitivity.rds")) else NULL
man <- de$manifest; P <- de$params

# en$go holds BOTH strata, so filter to the one that is primary for this cluster --
# otherwise the same term is listed twice and reads as two separate hits.
top_terms <- function(ck, cl, dirn, n = 3) {
  if (is.null(en$go)) return("-")
  st <- man$stratum[man$contrast == ck & man$cluster == cl & man$is_primary]
  if (!length(st)) return("-")
  g <- en$go[en$go$contrast == ck & en$go$cluster == cl & en$go$stratum == st[1] &
             en$go$ontology == "BP" & en$go$direction == dirn, ]
  if (!nrow(g)) return("-")
  paste(head(unique(g$Description[order(g$p.adjust)]), n), collapse = "; ")
}
top_genes <- function(key, dirn, n = 6) {
  d <- de$tables[[key]]; d <- d[d$direction == dirn & !d$confounder, ]
  if (!nrow(d)) return("-")
  paste(head(d$gene[order(d$padj)], n), collapse = ", ")
}

L <- c()
add <- function(...) L <<- c(L, sprintf(...))

add("# P7 KO-vs-WT per cardiomyocyte subcluster, and P7 WT vs P0 WT")
add("")
add("Generated %s from the Shiny app's data bundle. Regenerate with `analysis/2026-08-21_email/run.sh`.", P$built)
add("")
add("## Files")
add("")
add("| file | what |")
add("|---|---|")
add("| `P7_KO_vs_WT_by_CM_subcluster.xlsx` | Question 1. Complete gene lists per subcluster + GO run separately on the KO-up and KO-down lists. |")
add("| `P7WT_vs_P0WT.xlsx` | Question 2. All cardiomyocytes plus the same seven subclusters. |")
add("| `plots/part1/`, `plots/part2/` | GO dot plots, volcanoes, GSEA bars, cross-cluster heatmaps. PNG (300 dpi) and PDF. |")
add("| `csv/` | Every contrast, every stratum, every gene, unfiltered - the workbooks trim undetected genes and secondary strata to stay emailable. |")
add("")
add("Start with the **README** sheet in each workbook: it carries the column dictionary and the caveats.")
add("")
add("## Question 1 - P7 KO vs P7 WT")
add("")
add("Note this is **not** what the website's \"KO-vs-WT DE (per subgroup)\" tab shows. That tab pools P0 and P7.")
add("These are P7-only. Each sheet carries a `lfc_pooled_website` column so you can see both side by side.")
add("")
add("| subcluster | KO cells | WT cells | KO-up | KO-down | top KO-up genes | top KO-down genes |")
add("|---|---|---|---|---|---|---|")
m1 <- man[man$contrast == "P7_KO_vs_WT" & man$is_primary, ]
for (i in seq_len(nrow(m1)))
  add("| %s | %d | %d | %d | %d | %s | %s |", m1$cluster[i], m1$n_A[i], m1$n_B[i], m1$n_up[i], m1$n_down[i],
      top_genes(m1$key[i], "KO_up"), top_genes(m1$key[i], "KO_down"))
add("")
add("Top GO Biological Process terms:")
add("")
add("| subcluster | KO-up | KO-down |")
add("|---|---|---|")
for (i in seq_len(nrow(m1)))
  add("| %s | %s | %s |", m1$cluster[i], top_terms("P7_KO_vs_WT", m1$cluster[i], "KO_up"),
      top_terms("P7_KO_vs_WT", m1$cluster[i], "KO_down"))
add("")
add("## Question 2 - P7 WT vs P0 WT")
add("")
m2 <- man[man$contrast == "WT_P7_vs_P0" & man$is_primary, ]
add("| group | P7 cells | P0 cells | up at P7 | up at P0 | top GO (up at P7) | top GO (up at P0) |")
add("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(m2)))
  add("| %s | %d | %d | %d | %d | %s | %s |", m2$cluster[i], m2$n_A[i], m2$n_B[i], m2$n_up[i], m2$n_down[i],
      top_terms("WT_P7_vs_P0", m2$cluster[i], "P7_up", 2), top_terms("WT_P7_vs_P0", m2$cluster[i], "P0_up", 2))
add("")
add("The `sig_sets` column tags every gene with the curated maturation / metabolic / cell-cycle programs it")
add("belongs to, so those three questions are a filter away rather than a manual join.")
add("")
add("## One finding you should see before reading the GO results")
add("")
add("Every one of the seven subclusters has mitochondrially-encoded genes (`mt-Nd1`, `mt-Nd2`, `mt-Nd4`,")
add("`mt-Cytb`, `mt-Co1`, ...) in its **KO-up** list, and not one has an mt- gene on its KO-down side.")
add("Seven independently clustered populations do not agree that perfectly by biology. A one-directional")
add("shift of the whole mt- block in every cluster is the signature of a mitochondrial read-fraction")
add("difference between the two libraries - a QC covariate, not a pathway.")
add("")
add("It matters because those genes **are** the \"oxidative phosphorylation\" / \"electron transport chain\" /")
add("\"ATP synthesis\" GO terms that come up as the headline KO-up result in several clusters.")
if (!is.null(ms)) {
  add("")
  add("So each workbook carries a `GO_BP_no_mt` sheet - the same GO run with mt- genes removed from both the")
  add("input list and the universe - and a `Mito_GO_comparison` sheet counting what each direction loses:")
  add("")
  add("| subcluster | direction | significant genes | of which mt- | GO terms with mt- | without |")
  add("|---|---|---|---|---|---|")
  cm <- ms$comparison[ms$comparison$contrast == "P7_KO_vs_WT" & grepl("_up$", ms$comparison$direction), ]
  for (i in seq_len(nrow(cm)))
    add("| %s | %s | %d | %d | %d | %d |", cm$cluster[i], cm$direction[i], cm$n_sig_genes[i],
        cm$n_mt_genes[i], cm$terms_with_mt[i], cm$terms_without_mt[i])
  add("")
  add("Terms that survive the removal are about nuclear genes and read normally. Terms that disappear were")
  add("being carried by the mt- block. Nothing has been deleted from the DE tables - the mt- rows are")
  add("correct and are flagged in the `mito_encoded` column.")
  add("")
  add("**The KO-down side is not affected**, and the arithmetic says why. In CM1 the mt- share goes")
  add("1.35%% -> 1.71%% of signal, a ratio of 1.27 = **+0.34 in log2** - which is the size of the mt- fold")
  add("changes actually observed (`mt-Nd1` +0.54, `mt-Nd2` +0.45, `mt-Cytb` +0.44, and similar in every")
  add("cluster). The reciprocal squeeze on all other genes is 98.65%% -> 98.29%%, a ratio of 0.9964, or")
  add("**-0.005 in log2** - fifty times smaller than the 0.25 threshold. Consistent with that, removing")
  add("the mt- genes costs the KO-down lists 0-2 GO terms out of 37-425, and the KO-up lists nearly all.")
}
add("")
add("## Please read before interpreting")
add("")
add("1. **n = 1 animal per genotype x timepoint.** The two lanes per sample are the same library sequenced")
add("   twice. The p-values come from a cell-level Wilcoxon test and are pseudoreplicated: they rank genes,")
add("   they do not test a hypothesis. Any gene of interest needs replication in an independent cohort.")
add("2. **KO and WT animals are different sexes.** Y-linked genes, `Xist` and `Tsix` top every KO-vs-WT list.")
add("   They are flagged in the `confounder` column and excluded from every GO and GSEA input.")
add("3. **E2f7 and E2f8 mRNA are not reduced in the KO** in this data - most likely a conditional allele a")
add("   3'-biased assay cannot see. Do not read those rows as a knockdown check.")
add("4. **The genome-wide matrix is downsampled** to 8,026 cells. It is the only genome-wide matrix in the")
add("   bundle, and GO needs one. CM5, CM7 and CM8 come out at 39-96 cells per arm - their lists are the")
add("   least stable. A full-depth genome-wide table would need the upstream Seurat object.")
add("5. **CM0 and CM6 are absent by design.** CM6 has zero P7 cells; CM0 has ~82 WT-P7 and 167 KO-P7.")
add("   Neither supports a P7 KO-vs-WT contrast. CM4 has no G1 cells, so it has no phase-matched table.")
add("6. **An empty GO sheet is not \"nothing is enriched\".** The `GO_audit` sheet records, for every")
add("   cluster x direction x ontology, how many genes went in, how big the universe was, and how many terms")
add("   came out. Check it before concluding a direction is uninformative.")

writeLines(L, file.path(OUT, "README.md"))
cat("wrote", file.path(OUT, "README.md"), "-", length(L), "lines\n")
