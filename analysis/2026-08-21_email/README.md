# 2026-08-21 collaborator request — DE + GO as Excel

Answers two emailed questions that the Shiny app can display but cannot hand over
as files ("I can't download the excel from the website").

1. **P7 KO vs P7 WT** within cardiomyocyte subclusters **CM1, CM2, CM3, CM4, CM5, CM7, CM8** —
   complete gene lists, plus GO enrichment run **separately** on the KO-up and KO-down lists
   of each cluster, with figures.
2. **P7 WT vs P0 WT** — DEGs and enriched pathways, over all cardiomyocytes and the same
   seven subclusters, for a follow-up look at maturation, metabolism and cell-cycle genes.

## Why this isn't just an export of what the website already has

The website's **"KO-vs-WT DE (per subgroup)"** tab reads `app$tables$sub_DE[["res0.2"]]`
(`app.R:1200`), and those tables are KO vs WT **pooled across P0 and P7** — not the
P7-specific contrast that was asked for. The P7-specific numbers do exist, in
`app$fourgroup$de` from `build_fourgroup.R` on the `app/fourgroup-cm-analysis` branch, but
they are row-gated (`max(pct) >= 5% AND (padj < 0.05 OR |log2FC| >= 0.5)`), so they are not
*complete* lists, and GO was never run on them. Hence a recompute.

Every DE table here carries a `lfc_pooled_website` / `padj_pooled_website` column so the
pooled value the collaborator has already seen sits next to the P7-specific one.

## Running it

```bash
analysis/2026-08-21_email/run.sh
```

Outputs land in `deliverables/2026-08-21/` (git-ignored). The first run builds
`e2f-enrich:latest`; the Bioconductor layer takes a while.

| step | does |
|---|---|
| `01_de.R` | both contrasts × {all cells, G1} × clusters, **ungated**, genome-wide → `de_tables.rds` + one CSV each |
| `02_enrich.R` | GO BP/MF/CC per cluster per direction, Hallmark+KEGG GSEA, every figure → `enrich.rds` |
| `03_excel.R` | the two workbooks |

The bundle is mounted **read-only**. Unlike the `shiny_app/build_*.R` scripts, nothing here
rewrites `app_data.rds`.

## Method, and where it comes from

- **DE core** — `presto::wilcoxauc` on the log-normalised matrix, group `"A"` rows, `logFC > 0`
  means up in A. Identical to `deg_compute()` (`app.R:566`) and to `de_one()` in
  `build_fourgroup.R`, **minus the row gate**, so results are directly comparable to the site.
- **Matrix** — `app$deg_expr`, 24,221 genes × 8,026 cells. Genome-wide, which GO needs; the
  cell downsample is the price. The full-depth curated panel (2,181 genes × 30,030 cells) is
  re-run as a cross-check into the `*_fullcells` columns.
- **Signature sets** — parsed out of `shiny_app/build_signature_scores.R` at runtime rather
  than copied, so the two cannot drift.
- **GO** — `enrichGO(OrgDb = org.Mm.eg.db, keyType = "SYMBOL")`, BH, p < 0.05 / q < 0.2, set
  size 10–500, same call as `build_subcluster_enrichment.R`. The **universe** is the genes
  expressed in ≥ 5 % of at least one arm *in that cluster* — not all 24,221 (which would
  inflate every fold enrichment) and not the significant list.
- **GSEA** — `fgsea` over MSigDB mouse Hallmark + `C2 CP:KEGG_LEGACY`, ranked by
  `sign(log2FC) * -log10(p)`; the recipe in `build_subcluster_enrichment.R`.
- **Figure styling** — ported from `go_dotplot_gg()` (`app.R:279`) and `gsea_barplot_gg()`
  (`app.R:265`) so the figures match what the collaborator has been looking at.

## Design decisions worth knowing

- **Part 1 primary stratum is `all` cells; part 2's is `G1`.** Part 1 is unaffected by the
  sort either way — both its arms are P7. For part 2 we measured the confound rather than
  inheriting it: **within cardiomyocytes** the cycling fraction is 16.3 % (WT-P0) vs 25.0 %
  (WT-P7), a ratio of **1.53×**, not the 4.5–5.2× the project notes quote (that figure is not
  the P0-vs-P7 ratio inside the CM compartment of this bundle). Empirically, G1-matched and raw
  log2FCs correlate at **r = 0.99** with a median absolute difference of **0.001**, and no
  cell-cycle gene moves more than 0.07 log2 units between them. G1-matched stays the primary
  read as cheap insurance, but the two are not different biological answers, and the workbook
  says so rather than telling the reader to distrust a table that is in fact fine.
- **CM0 and CM6 are correctly absent from the request.** CM6 has *zero* P7 cells; CM0 has
  82 WT-P7 / 167 KO-P7. Neither supports a P7 KO-vs-WT contrast.
- **CM4 has zero G1 cells** — it is entirely S/G2M. It has no G1 table anywhere. Where that
  removes a cluster's primary stratum (part 2), the surviving stratum is promoted for figures
  and flagged `primary_fallback` in the manifest rather than silently substituted.
- **An empty GO sheet is not the same as "nothing is enriched".** Every cluster × direction ×
  ontology gets a row in `GO_audit` recording the input size, universe size, terms found, and
  which selection rule fired (the fallback ladder relaxes to `padj<0.05` alone, then to top-200
  by `|log2FC|`, when a list is under 10 genes). Placeholder figures say the same thing.
- **Confounders** (`Eif2s3y, Kdm5d, Uty, Ddx3y, Xist, Tsix, Gt(ROSA)26Sor`) stay in the DE
  tables, flagged, and are excluded from every GO and GSEA input.
- **Mitochondrially-encoded genes are a second confounder the project had not flagged.** They
  appear in the KO-up list of *all seven* subclusters and in the KO-down list of *none*. Seven
  independent populations do not agree that cleanly by biology; a one-directional shift of the
  whole `mt-` block is a read-fraction difference between libraries. It matters because those
  genes *are* the OXPHOS / electron-transport / ATP-synthesis GO terms that otherwise headline
  the KO-up result. `05_mito_sensitivity.R` measures the per-arm mitochondrial share and re-runs
  GO BP with `mt-` genes dropped from both the input and the universe; both versions ship, with
  a `Mito_GO_comparison` sheet showing what each direction loses. Nothing is deleted from the DE
  tables — the rows are correct and carry a `mito_encoded` flag.

## Limitations stated in the deliverable itself

n = 1 animal per genotype × timepoint, so cell-level Wilcoxon p-values are pseudoreplicated
and rank rather than test; KO and WT animals are different sexes; E2f7/E2f8 mRNA is not
reduced in the KO; the genome-wide matrix is downsampled, leaving CM5/CM7/CM8 at 39–96 cells
per arm. A full-depth genome-wide table would need the upstream Seurat object, which is not
on this machine.
