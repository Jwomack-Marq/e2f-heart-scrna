# Scripts — what each one does

Run order: **pipseeker_scripts.txt → scrublets.py → scRNA_lane_merge.Rmd (×4) →
scRNA_mergeP0.Rmd / scRNA_mergeP7.Rmd**. The `.Rmd` files set their working
directory to the project root automatically, so knit them from `scripts/`.

---

## 1. `pipseeker_scripts.txt` — preprocessing commands (reference only)
The PIPseeker v3.3.0 command lines used to turn raw FASTQ into count matrices,
one per sample/lane (8 total). Each aligns reads to **GRCm38** (STAR index),
chemistry "V", and writes `filtered_matrix/sensitivity_5/`, `metrics/`, and
`star/` into the sample folder. This is documentation of an already-completed
step — the raw FASTQ are not stored here, so it is not re-run as-is.

## 2. `scrublets.py` — doublet detection (Python)
Runs **Scrublet** on one sample/lane at a time to flag barcodes that are likely
two cells captured together.
- **Input:** `<sample>/filtered_matrix/sensitivity_5/` (matrix + features).
- **Output:** `<sample>/predicted_doublets.csv` (a `doublet_scores` and a
  True/False `predicted_doublets` column) and a score histogram PNG.
- **Run:** `Rscript`-free; `python scripts/scrublets.py -s P0WT_lane1` (repeat for
  all 8). Project root is auto-detected; override with `-p`.
- Key params: `expected_doublet_rate=0.08`, 30 PCs.

## 3. `scRNA_lane_merge.Rmd` — per-group QC + merge the two lanes (R)
The reusable template; set `samplename` to one of P0WT / P0KO / P7WT / P7KO and
knit, once per group (4 runs).
- **Input:** both lanes' `filtered_matrix/` + their `predicted_doublets.csv`.
- **Steps:** load each lane → compute mouse mito % (`^mt-`) → **QC filter**
  (genes ≥ 1500, a per-lane upper gene cap, mito ≤ `cutoff.mt`, and **remove
  predicted doublets**) → tag `$lane` → `merge()` the two lanes →
  `SCTransform` → PCA → cluster → UMAP → cluster markers (`FindAllMarkers`) →
  cell-type annotation (MAESTRO Heart_and_Aorta signatures).
- **Output:** `processing/merge.lanes.<group>.rds`,
  `results/figures/<group>cluster.umap.png`, and two marker tables in
  `results/markers/`.

## 4. `scRNA_mergeP0.Rmd` / `scRNA_mergeP7.Rmd` — merge conditions + compare KO vs WT (R)
One per timepoint. Combines the WT and KO objects for that timepoint and runs
the differential-expression comparison.
- **Input:** the two `processing/merge.lanes.<tp>{KO,WT}.rds` for that timepoint.
- **Steps:** `merge()` KO + WT → `SCTransform` → PCA → cluster → UMAP (with a
  lane-effect check) → cluster markers → cell-type annotation → cell-cycle
  scoring → **select cardiomyocytes by annotation/markers** →
  **KO-vs-WT differential expression**: an *exploratory* single-cell volcano
  **plus** a **pseudobulk DESeq2** test (the one to trust, valid only if the
  lanes are separate animals). `scRNA_mergeP0.Rmd` also runs **GO enrichment**
  on KO-up genes.
- **Output:** `processing/seurat.<tp>.merge.rds`,
  `results/markers/<tp>.cardiac.pseudobulk.DESeq2.csv`, volcano/GO figures.

---

### Running R on this machine
R 4.5.2 is installed at `C:\Program Files\R\R-4.5.2`. Either add its `bin\` to
PATH, or call it directly, e.g.
`& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e 'rmarkdown::render("scripts/scRNA_lane_merge.Rmd")'`.
Required packages incl. Seurat, glmGamPoi, dittoSeq, MAESTRO, clusterProfiler,
EnhancedVolcano, **DESeq2**, **babelgene** (see the main README for the full list).
