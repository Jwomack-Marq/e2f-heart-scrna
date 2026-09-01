#!/usr/bin/env Rscript
# Shared helpers/gene-sets for the our_analysis step scripts (04_/05_/06_).
# Source it via the tiny `_common.R` forwarder that sits in each step folder
# (it walks up to find this file). Sequential future plan + a bounded future-globals brake.
#
# Paths are anchored on the `our_analysis/` root, located by walking up from the
# running script until the `.projroot` sentinel file is found -- so scripts work
# at any folder depth and regardless of the caller's working directory.

suppressWarnings(suppressMessages({
  library(Seurat)
  library(dplyr)
  library(future)
}))
future::plan("sequential")
# 16 GiB, and the reasoning matters because the number was 3 GiB and that broke things.
#
# The limit exists to stop a large object being shipped to a parallel worker by accident.
# But the line above sets a SEQUENTIAL plan: there are no workers and nothing is ever
# transferred, so on this pipeline the check can only ever refuse to run something that
# would have been fine. And it did -- SCTransform on the 42,416-cell cardiomyocyte
# compartment needs 3.35 GiB for its chunking closure, so cm_subcluster_build.R died with
# "the total size of the 19 globals ... exceeds the maximum allowed size" on any machine
# with a clean environment. The shipped object predates that and was built elsewhere.
#
# Still bounded rather than Inf, per the original intent: if a future version of this
# pipeline does go parallel, a runaway transfer should still hit a ceiling.
options(future.globals.maxSize = 16 * 1024^3)

# ---- paths -----------------------------------------------------------------
.find_root <- function() {
  a0 <- commandArgs(FALSE); sp <- sub("^--file=", "", a0[grep("^--file=", a0)])
  d <- normalizePath(if (length(sp) && nzchar(sp)) dirname(sp) else getwd())
  while (!file.exists(file.path(d, ".projroot"))) {
    up <- dirname(d)
    if (identical(up, d)) stop("our_analysis/.projroot not found above ", d)
    d <- up
  }
  d
}
OUR_ROOT <- .find_root()
PROJ     <- OUR_ROOT                                  # back-compat alias
INPUT    <- file.path(OUR_ROOT, "01_input")           # copied raw PIPseeker matrices
PROC     <- file.path(OUR_ROOT, "processing")         # our intermediate Seurat .rds
RESULTS  <- file.path(OUR_ROOT, "results")
OUT      <- file.path(OUR_ROOT, "06_outputs")         # report/deck destination
OUTFIG   <- file.path(RESULTS, "figures")
OUTTAB   <- file.path(RESULTS, "tables")
for (d in c(PROC, OUT, OUTFIG, OUTTAB, file.path(RESULTS, "markers")))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

TIMEPOINTS <- c("P0", "P7")
merged_path <- function(tp) file.path(PROC, paste0("seurat.", tp, ".merge.rds"))

# ---- human -> mouse ortholog helper (same approach as the merge scripts) ----
to_mouse <- function(genes) {
  if (requireNamespace("babelgene", quietly = TRUE)) {
    o <- babelgene::orthologs(genes = genes, species = "mouse", human = TRUE)
    unique(o$symbol)
  } else paste0(substr(genes, 1, 1), tolower(substr(genes, 2, nchar(genes))))
}

# ---- cardiomyocyte selection (mirrors scRNA_mergeP0.Rmd marker-score path) --
# Returns a logical vector over cells: TRUE = cardiomyocyte cluster (top-quartile
# CM-marker score). Reproducible; no MAESTRO needed.
CM_MARKERS <- c("Tnnt2", "Myh6", "Actc1", "Nppa", "Myh7", "Ttn")
select_cardiac <- function(obj) {
  mk <- intersect(CM_MARKERS, rownames(obj))
  obj <- AddModuleScore(obj, features = list(mk), name = "CMsel")
  cl.score <- tapply(obj$CMsel1, obj$seurat_clusters, median)
  cm.clusters <- names(cl.score)[cl.score > quantile(cl.score, 0.75, na.rm = TRUE)]
  obj$seurat_clusters %in% cm.clusters
}

# ---- canonical cardiac cell-type markers (mouse) for annotation/dotplots ----
CELLTYPE_MARKERS <- list(
  Cardiomyocyte   = c("Tnnt2","Myh6","Actc1","Ttn","Tnni3","Nppa"),
  Fibroblast      = c("Col1a1","Col1a2","Dcn","Pdgfra","Gsn"),
  Endothelial     = c("Pecam1","Cdh5","Kdr","Fabp4","Egfl7"),
  Endocardial     = c("Npr3","Plvap"),
  Immune_Myeloid  = c("Ptprc","Cd68","Lyz2","C1qa","Csf1r"),
  Mural_Pericyte  = c("Rgs5","Pdgfrb","Myh11","Acta2","Tagln"),
  Epicardial      = c("Wt1","Msln","Upk3b"),
  Neuronal_Glial  = c("Plp1","Kcna1","Sox10"),
  RBC             = c("Hba-a1","Hbb-bs")
)

# ---- E2F / cell-cycle target genes (mouse symbols) --------------------------
# E2F7/E2F8 are repressors of classic activating-E2F (cell-cycle) targets, so KO
# should DE-REPRESS these. Curated canonical set + Seurat's S/G2M lists.
E2F_TARGETS <- c("Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7","Pcna","Cdc6","Cdt1",
                 "Ccne1","Ccne2","Ccna2","Ccnb1","Ccnb2","Cdk1","Cdc20","E2f1",
                 "Mki67","Top2a","Birc5","Bub1","Rrm2","Tk1","Cdc25a","Foxm1",
                 "Aurkb","Plk1","Cenpa","Cenpe")
cc_lists <- function() {
  list(S   = to_mouse(Seurat::cc.genes.updated.2019$s.genes),
       G2M = to_mouse(Seurat::cc.genes.updated.2019$g2m.genes))
}

# ---- CM maturation markers (ploidy surrogate context) -----------------------
# Polyploidization tracks CM maturation; report these as INDIRECT surrogates.
CM_MATURE   <- c("Myh6","Tnni3","Pln","Atp2a2","Ckm","Myl2","Cox6a2")
CM_IMMATURE <- c("Myh7","Tnni1","Nppa","Nppb","Ccnd1","Mki67","Top2a")

genotype_of <- function(orig) ifelse(grepl("KO$", orig), "KO", "WT")

cat(sprintf("[_common] OUR_ROOT=%s\n", OUR_ROOT))
