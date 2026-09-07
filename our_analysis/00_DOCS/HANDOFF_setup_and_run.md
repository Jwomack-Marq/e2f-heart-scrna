# Handoff: set up & run the scRNA re-analysis on the dedicated 32 GB laptop

Goal: regenerate the Seurat clustering + marker + DE outputs into `rerun/` so they
can be compared against the original run — without crashing (the previous attempt
on a busy 32 GB laptop ran out of RAM). This machine is dedicated to the job, so
the main requirements are: **install `presto`** (prevents the FindAllMarkers
memory blow-up) and **run one sample at a time**.

> Context for whoever runs this: the project lives in OneDrive, so this folder
> (`Han_scRNA_2025/`) with all scripts and data should already be synced here.
> Confirm the 8 sample folders (e.g. `P0WT_lane1/filtered_matrix/sensitivity_5/`)
> are fully downloaded (not just cloud placeholders) before starting.

---

## 0. Prerequisites to confirm/install

1. **R 4.5.x** (the other machine used 4.5.2). Check:
   `& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" --version`
   (adjust the path/version if different.)
2. **OneDrive fully synced** — right-click the `Han_scRNA_2025` folder →
   "Always keep on this device" so the matrices are real files, not placeholders.
3. These CRAN/Bioconductor packages (install once — see step 1 below):
   Seurat, SeuratObject, glmGamPoi, dplyr, tibble, patchwork, dittoSeq,
   clusterProfiler, org.Mm.eg.db, enrichplot, EnhancedVolcano, DT, DESeq2,
   babelgene, rmarkdown, knitr, future, **presto**.
   - MAESTRO is intentionally **skipped** (deprecated, hard to build). The
     scripts already guard it and fall back to marker-based cell selection.

---

## 1. Install the packages (one time)

Run from a terminal (adjust the Rscript path if needed):

```powershell
$rs = "C:\Program Files\R\R-4.5.2\bin\Rscript.exe"

# Bioconductor + CRAN stack
& $rs -e "if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('glmGamPoi','dittoSeq','clusterProfiler','enrichplot','EnhancedVolcano','DESeq2','babelgene'), update=FALSE, ask=FALSE)"

# presto — install the prebuilt binary from R-universe (no compiler needed).
& $rs -e "install.packages('presto', repos=c('https://immunogenomics.r-universe.dev','https://cloud.r-project.org'))"
```

Verify everything is present (should print OK for each):

```powershell
& $rs -e "p<-c('Seurat','glmGamPoi','dittoSeq','clusterProfiler','enrichplot','EnhancedVolcano','DESeq2','babelgene','presto','rmarkdown','knitr'); i<-rownames(installed.packages()); for(x in p) cat(sprintf('%-16s %s\n',x,ifelse(x%in%i,'OK','MISSING')))"
```

If `presto` shows MISSING (R-universe binary unavailable for this R version) and
Rtools is installed, build it from source — **note this compiles external code, so
it requires explicit approval**:
`& $rs -e "if(!requireNamespace('remotes',quietly=TRUE)) install.packages('remotes'); remotes::install_github('immunogenomics/presto', upgrade='never')"`

---

## 2. Free up memory before running

- Close OneDrive sync (pause it), browsers, Teams, Outlook, etc.
- Confirm free RAM is high:
  `Get-CimInstance Win32_OperatingSystem | %% { "{0:N1} GB free" -f ($_.FreePhysicalMemory/1MB) }`
  Aim for **>25 GB free** before starting.

---

## 3. Run the pipeline — ONE sample at a time

The runner `rerun/run_rerun.R` already: sets the working dir to the project root,
sends output to `rerun/`, forces a **sequential** future plan, and keeps a 3 GiB
safety brake. Run each command and let it finish before starting the next
(each is one isolated R process, so memory is released between them).

```powershell
$rs   = "C:\Program Files\R\R-4.5.2\bin\Rscript.exe"
$run  = "$PWD\rerun\run_rerun.R"     # run from inside the Han_scRNA_2025 folder

# Step A — lane merges (4 groups). ~5-15 min each.
& $rs $run lane P0WT
& $rs $run lane P0KO
& $rs $run lane P7WT
& $rs $run lane P7KO

# Step B — condition merges + DE (needs the 4 .rds from Step A). ~15-30 min each.
& $rs $run merge P0
& $rs $run merge P7
```

Watch the console: each ends with `=== DONE ... in N min ===`. If any step is
killed for memory, close more apps and re-run just that one command.

---

## 4. What success looks like (compare to original)

New files appear under `rerun/` (originals stay untouched at the project root /
`processing/` / `results/`):

- `rerun/processing/merge.lanes.{P0WT,P0KO,P7WT,P7KO}.rds`
- `rerun/processing/seurat.{P0,P7}.merge.rds`
- `rerun/results/markers/*.allmarkers.csv`, `*.cluster.top30.markers.csv`
- `rerun/results/markers/{P0,P7}.cardiac.pseudobulk.DESeq2.csv`  ← new pseudobulk DE
- `rerun/results/figures/*cluster.umap.png`
- `rerun/qc_comparison.csv`  ← already generated (doublet+mito filtering impact)

Then I (Claude) can diff the rerun marker/DE tables and cell counts against the
original `results/` to summarize what changed.

---

## 5. Things to keep in mind

- **Seurat version:** this machine likely has Seurat v5; the original run was v4.
  So numbers won't be bit-identical even aside from our QC fixes — some difference
  is just the version (clustering/SCTransform defaults). The scripts include a v5
  `JoinLayers` shim so they run correctly under v5.
- **MAESTRO skipped:** cell-type annotation by MAESTRO is bypassed; cardiomyocytes
  are selected by canonical markers (Tnnt2/Myh6/Actc1/Nppa) instead — reproducible
  and doesn't need MAESTRO.
- **Pseudobulk DE caveat unchanged:** still pending lab confirmation of whether
  lane1/lane6 are separate animals. Based on the FASTQ paths they appear to be the
  *same* library on two flow-cell lanes (n=1 per condition), which would make the
  KO-vs-WT statistics descriptive only.
- **Never set `future.globals.maxSize` to `Inf`** — that's what caused the crash.
  The 3 GiB brake is intentional.
