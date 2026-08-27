# Replicate structure — resolved: n = 1 animal per condition

**Bottom line:** `lane1` and `lane6` are the **same sequencing library run on two
flow-cell lanes** (technical replicates), **not** two animals. There is therefore
**n = 1 biological replicate per condition**, and KO-vs-WT differential expression
**cannot be given valid p-values** by any method. All DE here is **descriptive /
hypothesis-generating only.**

## Evidence (from the per-sample `run_config.csv`)

Each condition's two "lanes" point to the **identical Novogene library ID**; only the
parent directory differs (`Raw_data/` vs `Raw_data_Lane6/`):

| Condition | lane1 `fastq` | lane6 `fastq` |
|-----------|---------------|---------------|
| P0WT | `…/Raw_data/01.RawData/P0WT/P0WT_CKDL250003755` | `…/Raw_data_Lane6/01.RawData/P0WT/P0WT_CKDL250003755` |
| P7KO | `…/Raw_data/01.RawData/P7KO/P7KO_CKDL250003754` | `…/Raw_data_Lane6/01.RawData/P7KO/P7KO_CKDL250003754` |

A shared `CKDL########` ID = one library prep (one animal / one PIPseq reaction)
demultiplexed from two flow-cell lanes. (Confirm the same pattern for the remaining
conditions by diffing the `fastq` field across each pair of `run_config.csv` files.)

## Why this matters

- Technical replicates share the biological sample, so within-group variance is
  near-zero. Pseudobulk DESeq2 over `lane1`/`lane6` therefore reports thousands of
  genes at `padj ≈ 0` — an artifact of pseudoreplication, not real DE confidence.
  In the rerun: **P0 = 3,777** and **P7 = 3,377** genes at `padj < 0.05`, many at
  `padj = 0`. A male-specific gene (`Ddx3y`) among P0's top hits is a further tell
  that the single KO and single WT animals simply differ (incl. possibly by sex).
- No statistical fix (pseudobulk, mixed models, etc.) recovers valid condition-level
  inference from n = 1. More biological replicates are the only remedy.

## What is still valid from the current data

Treat the current dataset as a **pilot**. The following are legitimate:

- Cell clustering, **cell-type annotation**, and within-sample structure.
- **Descriptive** KO-vs-WT contrasts ranked by **effect size** (shrunken log2FC),
  cross-checked for **consistency across P0 and P7** and against **known E2F biology**.
  (Shared KO-up signal in the pilot: e.g. **Tcf4, Gabbr2, Adamts9** at both timepoints.)
- **Composition description** (cell-type proportions KO vs WT) — *reported*, not tested.
- Hypothesis readouts that don't need p-values: KO verification (`E2f7`/`E2f8`
  reduced in KO), E2F/cell-cycle target **de-repression** scores, cycling-CM fraction,
  and ploidy/maturation **surrogates**.

Everything DE-related must be labeled **"descriptive (n = 1)"** in figures/tables.

## Replicated-cohort spec (for publishable inference)

- **≥ 3 animals per condition** (P0WT, P0KO, P7WT, P7KO), each a **separately prepared
  library** — re-sequencing one library on more lanes does **not** add biological n.
- Balance **sex** across genotypes (or record and model it); the `Ddx3y` signal warns
  that sex may currently be confounded with genotype.
- Power: with ~3–4 biological replicates/group, pseudobulk DESeq2 gives proper
  condition-level p-values; this pilot's effect-size-ranked lists define the priority
  genes/pathways to confirm.
- Keep the same upstream chemistry/PIPseeker settings for comparability.

## Lab-confirmed mechanism + a processing issue it exposes (June 2026)

The lab confirmed the two lanes are the **same library sequenced twice to add depth**
("the same sample/library was sequenced twice on two different lanes... reads from
both lanes for the same sample can be pooled together"). This confirms n=1 and adds a
processing caveat:

- **Barcode overlap lane1↔lane6 is 97–100%** (`rerun/analysis/tables/lane_barcode_overlap.csv`):
  P0WT 99.7%, P0KO 99.7%, P7WT 100.0%, P7KO 99.5%. So essentially *every* cell appears
  in both lanes.
- **The current pipeline counted each lane separately (separate PIPseeker runs) and
  `merge()`d them**, so nearly every cell is represented **twice**, each copy carrying
  only ~half its true depth. Reported cell counts are ~2× inflated (e.g. P0WT "18,924"
  ≈ ~9,500 distinct cells), and `nFeature ≥ 1500` QC was applied to half-depth cells
  (over-aggressive).
- **Standard fix:** combine the two lanes *per library* at the read level — re-run
  PIPseeker giving it BOTH lanes' FASTQs as one sample (it concatenates reads and
  dedupes UMIs across lanes → one matrix/library, each cell once at full depth), then
  re-run the Seurat pipeline on 4 libraries. The FASTQs live on the cluster
  (`/projects/rpci/tliu4/...`), not locally. A local matrix-level pooling (sum shared
  barcodes) is a quick approximation but does NOT dedupe UMIs across lanes (mild
  over-count of abundant genes).
- **Impact on existing results:** the descriptive DE conclusions are robust (pseudobulk
  SUMS over all cells are ~conserved whether cells are split or pooled — n=1, sex
  confound, ROSA26, candidate genes like Gabbr2 all stand). What needs regenerating on
  properly-combined data: cell counts, QC thresholds, clustering/UMAP, doublet calls,
  composition, and per-cell readouts.

**Re-count spec:** `scripts/recount_combine_lanes.sh` (run on the cluster) — one
PIPseeker run per library over both lanes' FASTQs → `<SAMPLE>_combined/`.

**Downstream Seurat changes after re-counting:**
- `scRNA_lane_merge.Rmd` becomes **per-library** (no lane merge): read the single
  `<SAMPLE>_combined/filtered_matrix/sensitivity_5/` matrix, QC, SCTransform, cluster,
  markers, save. Re-derive doublets ONCE per library (scDblFinder on full depth).
- **Re-check the `nFeature ≥ 1500` cutoff** on full-depth cells (it was tuned on
  half-depth, so it likely kept too few cells) — inspect the violins.
- `scRNA_mergeP0/P7.Rmd` are structurally unchanged (merge KO+WT), but the per-cell
  `lane` field disappears, so the pseudobulk "lanes as replicates" path no longer fires
  (`table(condition) >= 2` is false at 1 sample/condition) — it correctly falls back to
  a **descriptive** KO-vs-WT fold-change with no p-values, making n=1 explicit.
- Expect cell counts to ~halve to the true distinct-cell numbers (e.g. P0WT ≈ 9.5k).

## Pilot findings that further constrain interpretation (from rerun/analysis/)

Two results from the descriptive analyses sharpen the caveats above and must shape
any replicated cohort:

1. **The KO and WT animals are different sexes (confound on top of n=1).** The top
   KO-up genes by shrunken log2FC are Y-linked: **Eif2s3y, Kdm5d, Uty, Ddx3y** (at
   *both* P0 and P7). So a large part of the apparent KO-vs-WT signal is simply
   male-vs-female, not genotype. → The replicated cohort **must sex-match** genotypes
   (or balance and model sex).

2. **The KO is not a clean transcript-null, and an engineered ROSA26 allele is
   present.** `E2f7`/`E2f8` are **not reduced** in KO (flat-to-higher: e.g. P7
   E2f8 0.038→0.072), and **`Gt(ROSA)26Sor`** is among the shared KO-up genes —
   consistent with a Cre/reporter knock-in at ROSA26 and an exon-specific
   *conditional* deletion that 3′ scRNA cannot detect (plus E2F7/8 autorepression
   de-repressing the residual locus). → **Confirm the exact KO allele / Cre driver
   with the lab**; a clean 3′ readout of knockout may require targeting the deleted
   exons or a reporter.

After removing the sex (Y-linked) and construct (`Gt(ROSA)26Sor`) genes, the
*biological* candidates consistent across P0 and P7 are: **Gabbr2, Tcf4, Ralyl,
Adamts9, Atp6v0e2, Arhgap36, Stbd1** — the priority list to confirm in a replicated,
sex-matched cohort.

## Status of the original open question

This supersedes the "**Confirm before trusting DE**" caveat in `README.md`: the
question is now answered from the run metadata (n = 1). Lab confirmation is welcome but
the library IDs are dispositive that the two lanes are the same library.
