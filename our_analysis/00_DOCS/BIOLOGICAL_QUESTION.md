# Can we answer the biological question?

**Question (analysis/README.md):** *"Does loss of E2F7/E2F8 disrupt cardiomyocyte cell-cycle exit,
maturation, and binucleation between P0 and P7?"*

**Bottom line:** Of the three axes, **cell-cycle exit** and **maturation** are *addressable
descriptively* with this data (and we have a hypothesis-consistent lead for cell-cycle exit).
**Binucleation is NOT measurable by this assay at all.** And every KO-vs-WT comparison is **descriptive
only** because the design is **n = 1 per condition with sex confounded with genotype**. A definitive
answer needs a new experiment (below).

---

## 1. Cell-cycle exit — ADDRESSABLE (descriptive); strongest signal

What we have (`cell_cycle.R`, `e2f_readouts.R`, `tf_activity.R`, `pathway_msigdb.R`, `cm_subtypes.R`):
- Cell-cycle phase scoring; cycling (S/G2M) fraction per cell type and **per CM subtype**, KO vs WT.
- **Cycling-cardiomyocyte compartment expanded ~1.4× in P7 KO** (25% vs 18%; propeller); CM cycling
  fraction 31.6% (KO) vs 25.6% (WT) at P7 — and **specific to P7**, not P0.
- **E2F-target de-repression in KO** (Hallmark E2F_TARGETS enriched; E2F-target regulon activity up in
  KO cardiomyocytes) — the expected direction for loss of repressor E2Fs.
- Ridge plots of S/G2M scores by cluster (parity with the prior analysis).

**Verdict:** We can characterize cell-cycle exit, and the data *trend* toward the predicted phenotype
(KO cardiomyocytes retain cycling at P7). Descriptive only.

## 2. Maturation — ADDRESSABLE (descriptive); KO effect weak here

What we have (`trajectory_slingshot.R`, `cross_timepoint.R`, CM mature/immature module scores):
- Slingshot **maturation pseudotime** (validated: tracks mature markers up, cycling down).
- P7 cardiomyocytes are more mature than P0 (expected); **KO vs WT pseudotime distributions overlap**
  within each timepoint — no clear maturation *delay* in KO along this axis.
- Cross-timepoint (P0→P7) DE recovers the maturation program (Sln/Bmp10/Myl7 down, etc.).
- **CM subtypes** (ventricular/atrial/trabecular/compact/cycling) now assigned, with composition KO vs WT.

**Verdict:** Maturation is well characterized; the KO-vs-WT maturation difference is ≈ null in this
pilot (descriptive). (This exceeds the prior pipeline, which had no pseudotime.)

## 3. Binucleation — NOT ADDRESSABLE with this data

- Binucleation (two nuclei per cardiomyocyte, from karyokinesis without cytokinesis) and polyploidy are
  **nuclear/DNA-content phenotypes**. 3′ droplet scRNA measures the transcriptome per droplet — it
  **cannot count nuclei per cell or measure DNA ploidy.**
- Cell-cycle/“cycling” markers do **not** equal binucleation: a cycling cell may divide (mononuclear
  proliferation), binucleate (cytokinesis failure), or endoreplicate (polyploidy) — scRNA can't
  distinguish these. Our "ploidy surrogates" (UMI/gene counts, maturation markers) are **indirect** and
  were flagged as such. (The prior MCW pipeline also has no binucleation/ploidy readout.)
- **What's actually required:** an orthogonal assay — confocal imaging of isolated cardiomyocytes with
  nuclei counted (e.g., cTnT/PCM1 + DAPI), or flow-cytometry DNA-content/ploidy, or comparable methods.

**Verdict:** Out of scope for scRNA; needs a dedicated experiment.

---

## Overarching blockers (apply to all three axes)

1. **n = 1 biological replicate per condition** (lane1/lane6 are the same library) → no valid KO-vs-WT
   p-values; everything is descriptive. (See `REPLICATES.md`.)
2. **Sex confounded with genotype** (KO = male, WT = female) → KO-vs-WT differences cannot be separated
   from sex differences.
3. **KO not transcript-confirmed** (E2f7/E2f8 not reduced in KO — likely a conditional/ROSA26 allele
   invisible to 3′ scRNA) → the genotype effect itself needs verification.

## To actually answer the question

- **Replicated, sex-matched cohort:** ≥3 animals/condition, separately prepared libraries (not
  re-sequenced lanes), sex balanced/recorded → valid pseudobulk DE and cycling/composition statistics.
- **Confirm the KO** with an exon-targeted assay (RT-qPCR across floxed exons, long-read, or protein).
- **Add a binucleation/ploidy assay** (imaging or flow) — scRNA cannot supply this.
- **Priority descriptive leads to validate** (from this pilot): the P7 KO cycling-cardiomyocyte expansion
  (cell-cycle exit), E2F-target de-repression in KO cardiomyocytes, and the candidate genes
  Gabbr2 / Tcf4 / Adamts9 (after removing the sex/ROSA26 confounders).

## Coverage vs. the prior MCW pipeline

We now reproduce/extend all of their `01`–`08` analyses: QC, Harmony integration, annotation, cluster
markers, **cardiomyocyte subtyping**, pseudobulk DE (+lfcShrink), pathway/GSEA, abundance (propeller),
cell-cycle scoring (+ ridge plots), and the E2F atlas (FeaturePlot + family DotPlot). Binucleation is the
one axis neither pipeline can address from scRNA.
