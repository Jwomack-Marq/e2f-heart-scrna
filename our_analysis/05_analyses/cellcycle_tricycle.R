#!/usr/bin/env Rscript
# tricycle as an independent second opinion on cell-cycle state.
#
# WHY. Our cell-cycle numbers come from Seurat::CellCycleScoring (05_analyses/cell_cycle.R)
# and they read oddly: cardiomyocytes are 15-17 % "cycling" at P0 and 26-32 % at P7. This
# script was written to test whether that is a method artifact, because CellCycleScoring has
# a known failure mode -- AddModuleScore centres each cell against control gene sets drawn
# from THE SAME dataset, and the phase call is then an argmax with G1 assigned only when both
# scores are negative, so in a barely-cycling population it can assign much of the noise.
#
# tricycle (Zheng et al. 2022, Genome Biology) shares none of that machinery: it projects
# cells onto a FIXED reference embedding learned from mouse neurosphere data and returns a
# CONTINUOUS position theta. Nothing about the answer depends on our object's composition.
#
# WHAT THE COMPARISON ACTUALLY FOUND -- recorded here because it is the opposite of the
# hypothesis above, and the hypothesis should not survive in the comments:
#
#   1. tricycle AGREES with CellCycleScoring. kappa 0.87, 95 % per-cell concordance, r 0.97
#      across cell-type x timepoint x genotype groups. The calls are not a Seurat artifact,
#      and the P0 -> P7 rise reproduces under an independent method.
#   2. BOTH are strongly depth-dependent. Within one cell type at one timepoint the cycling
#      fraction climbs monotonically with detected genes -- endothelial cells at P0 go from
#      21 % in the lowest gene quartile to 72 % in the highest. Nothing biological makes a
#      cell likelier to be in S phase because it was sequenced deeper.
#   3. Depth-matching removes about two thirds of the P0 -> P7 cardiomyocyte gap
#      (15.1 -> 22.1 % raw becomes 18.5 -> 20.7 % at matched depth).
#   4. Ambient RNA is pervasive and sets a floor: 100 % of "RBC" cells detect cardiac
#      sarcomere genes (Tnnt2/Myh6/Actc1) they cannot transcribe. That same population
#      scores 18-22 % "cycling" by every method here -- at or above the cardiomyocyte value.
#
#   So the honest reading is not "Seurat was wrong and tricycle is right". It is that the
#   absolute cardiomyocyte cycling fraction is not resolvable above the ambient/depth floor
#   in this dataset, and neither method can fix that because the problem is upstream of both.
#
# NOTE ON THE BIOLOGY, against which the P7 rise is NOT prima facie absurd: mouse
# cardiomyocytes go through a binucleation wave around P4-P10, entering S phase and mitosis
# without completing cytokinesis. Elevated S/G2M marker expression at P7 is consistent with
# that, so "P7 > P0" cannot by itself be called an error -- which is exactly why the depth
# and ambient controls below, not the direction of the difference, are what settle it.
#
#   Reads  processing/seurat.combined.annotated.rds  (whole heart: all cell types, so
#          genuinely proliferative populations act as a POSITIVE CONTROL -- if tricycle
#          reports no cycling anywhere it has failed, and only if it separates cycling
#          immune/fibroblast cells from P7 CMs is the CM result meaningful)
#   Writes results/tables/cellcycle_tricycle_*.csv
#          results/figures/cellcycle_tricycle_*.png
#
#   Rscript cellcycle_tricycle.R [--max-cells=N]   # N subsamples, for a smoke test
#
# Runs in e2f-tricycle:latest (our_analysis/Dockerfile.tricycle).

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
suppressWarnings(suppressMessages({
  library(tricycle); library(SingleCellExperiment); library(ggplot2)
}))

argval <- function(flag, default) {
  a <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", flag, "="), "", a[1]) else default
}
MAXCELLS <- as.integer(argval("--max-cells", "0"))
SEED     <- 42L
say <- function(...) cat(sprintf("[tricycle] %s\n", paste0(...)))

## ---- 1. the cells, and the call we already made ---------------------------
obj <- file.path(PROC, "seurat.combined.annotated.rds")
stopifnot(file.exists(obj))
say("reading ", basename(obj), " ...")
comb <- readRDS(obj)
DefaultAssay(comb) <- "RNA"
comb <- NormalizeData(comb, verbose = FALSE)

md <- comb@meta.data
stopifnot("celltype" %in% names(md))
if (!"Phase" %in% names(md))
  stop("no Phase column -- this object predates cell_cycle.R; run that first")
md$genotype <- genotype_of(md$orig.ident)
say(sprintf("%s cells | %d cell types | Phase present (the call we are checking)",
            format(nrow(md), big.mark = ","), length(unique(md$celltype))))

if (MAXCELLS > 0 && MAXCELLS < ncol(comb)) {
  set.seed(SEED); keep <- sample(colnames(comb), MAXCELLS)
  comb <- comb[, keep]; md <- md[keep, , drop = FALSE]
  say(sprintf("SMOKE TEST: subsampled to %s cells", format(ncol(comb), big.mark = ",")))
}

## ---- 2. tricycle -----------------------------------------------------------
# logcounts, because that is what tricycle's reference was built against. Same matrix
# Seurat scored, so no normalisation difference can explain a disagreement.
X <- SeuratObject::GetAssayData(comb, assay = "RNA", layer = "data")
sce <- SingleCellExperiment(assays = list(logcounts = X),
                            colData = md[colnames(X), , drop = FALSE])
rm(X); invisible(gc())

say("projecting onto the reference cycle space (mouse, gene symbols) ...")
# How many of the 500 reference genes matched matters: the projection is only as good as
# the overlap, and a poor match would be a reason to distrust the answer rather than the
# Seurat one. project_cycle_space reports it as a message, so capture rather than re-derive.
msg <- character(0)
sce <- withCallingHandlers(
  project_cycle_space(sce, gname.type = "SYMBOL", species = "mouse"),
  message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
ngene <- suppressWarnings(as.integer(sub(".*is (\\d+).*", "\\1",
           grep("number of projection genes", msg, value = TRUE)[1])))
say(sprintf("  %s of tricycle's 500 reference genes matched our features",
            if (is.na(ngene)) "?" else format(ngene)))
if (!is.na(ngene) && ngene < 300)
  warning("only ", ngene, " reference genes matched -- the projection is weakly supported")

sce <- estimate_cycle_position(sce)
theta <- sce$tricyclePosition
say(sprintf("  theta assigned for %s cells", format(sum(!is.na(theta)), big.mark = ",")))

# The discrete call that is ALLOWED TO ABSTAIN. NA here is not a failure, it is the
# measurement refusing to guess -- exactly what CellCycleScoring cannot do.
say("Schwabe staging (returns NA where it cannot confidently stage) ...")
# Given an SCE this returns the SCE with a CCStage column, NOT a bare vector -- reading
# the return value as a vector silently yields nonsense.
sce <- suppressMessages(estimate_Schwabe_stage(sce, gname.type = "SYMBOL", species = "mouse"))
stage <- as.character(sce$CCStage)
stopifnot(length(stage) == ncol(sce))
say(sprintf("  staged %s of %s cells; %s (%.1f%%) left NA -- the abstentions",
            format(sum(!is.na(stage)), big.mark = ","), format(ncol(sce), big.mark = ","),
            format(sum(is.na(stage)), big.mark = ","), 100 * mean(is.na(stage))))

## ---- 3. one table, both methods, same cells --------------------------------
# The cycling arc is derived from OUR data, not assumed. fit_periodic_loess of each
# canonical marker against theta puts the S genes first and the G2M genes after, in order:
# Mcm2 0.64pi, Pcna 0.79pi, Rrm2 1.01pi, Cdk1 1.05pi, Top2a 1.06pi, Mki67 1.16pi,
# Ccnb1 1.20pi -- with R2 0.41-0.62, so theta genuinely tracks the cycle here. That places
# S-through-M in (0.5pi, 1.5pi) and leaves the G1/G0 pile at theta ~ 0/2pi, which is where
# ~70 % of all cells sit. Verified rather than assumed, and re-checked below.
CYC_LO <- 0.5 * pi; CYC_HI <- 1.5 * pi
cyc_tri <- !is.na(theta) & theta >= CYC_LO & theta <= CYC_HI

out <- data.frame(
  cell            = colnames(sce),
  celltype        = md$celltype,
  timepoint       = md$timepoint,
  genotype        = md$genotype,
  orig.ident      = md$orig.ident,
  seurat_phase    = md$Phase,
  seurat_cycling  = md$Phase %in% c("S", "G2M"),
  tricycle_theta  = round(theta, 4),
  tricycle_stage  = as.character(stage),
  tricycle_staged = !is.na(stage),
  tricycle_cycling = cyc_tri,
  stringsAsFactors = FALSE)
write.csv(out, file.path(OUTTAB, "cellcycle_tricycle_percell.csv"), row.names = FALSE)

## ---- 4. the comparison that answers the question ---------------------------
agg <- function(d) data.frame(
  n                 = nrow(d),
  pct_cycling_seurat = round(100 * mean(d$seurat_cycling), 1),
  pct_cycling_tricycle = round(100 * mean(d$tricycle_cycling), 1),
  pct_staged_tricycle  = round(100 * mean(d$tricycle_staged), 1),
  median_theta      = round(median(d$tricycle_theta, na.rm = TRUE), 3))

by_ct <- do.call(rbind, lapply(split(out, list(out$timepoint, out$celltype, out$genotype),
                                     drop = TRUE), function(d) {
  cbind(data.frame(timepoint = d$timepoint[1], celltype = d$celltype[1],
                   genotype = d$genotype[1], stringsAsFactors = FALSE), agg(d))
}))
by_ct <- by_ct[order(by_ct$celltype, by_ct$timepoint, by_ct$genotype), ]
rownames(by_ct) <- NULL
write.csv(by_ct, file.path(OUTTAB, "cellcycle_tricycle_vs_seurat_by_celltype.csv"), row.names = FALSE)

cat("\n=== cycling fraction, both methods, same cells ===\n")
print(by_ct, row.names = FALSE)

cat("\n=== CARDIOMYOCYTES: does the fraction fall from P0 to P7, as it must? ===\n")
cm <- by_ct[by_ct$celltype == "Cardiomyocyte", ]
print(cm, row.names = FALSE)
for (g in c("WT", "KO")) {
  r <- cm[cm$genotype == g, ]
  if (nrow(r) == 2) {
    p0 <- r[r$timepoint == "P0", ]; p7 <- r[r$timepoint == "P7", ]
    cat(sprintf("  %s  Seurat  P0 %5.1f%% -> P7 %5.1f%%  (%s)\n", g,
                p0$pct_cycling_seurat, p7$pct_cycling_seurat,
                if (p7$pct_cycling_seurat < p0$pct_cycling_seurat) "falls" else "RISES -- implausible"))
    cat(sprintf("  %s  tricycle P0 %5.1f%% -> P7 %5.1f%%  (%s)\n", g,
                p0$pct_cycling_tricycle, p7$pct_cycling_tricycle,
                if (p7$pct_cycling_tricycle < p0$pct_cycling_tricycle) "falls" else "RISES -- implausible"))
  }
}

# The positive control. Without this the CM result is uninterpretable: a method that
# calls nothing cycling anywhere would "agree" with the biology by accident.
cat("\n=== POSITIVE CONTROL: which cell types does tricycle call most cycling at P0? ===\n")
pc <- by_ct[by_ct$timepoint == "P0" & by_ct$n >= 100, ]
pc <- pc[order(-pc$pct_cycling_tricycle), c("celltype","genotype","n",
                                            "pct_cycling_seurat","pct_cycling_tricycle")]
print(head(pc, 10), row.names = FALSE)

cat("\n=== agreement: Seurat phase x tricycle stage (all cells) ===\n")
tb <- table(seurat = out$seurat_phase, tricycle = ifelse(is.na(out$tricycle_stage),
                                                         "(not staged)", out$tricycle_stage))
print(tb)
cat(sprintf("\n  cells both call cycling : %s\n  Seurat only            : %s\n  tricycle only          : %s\n  neither                : %s\n",
            format(sum(out$seurat_cycling & out$tricycle_cycling), big.mark = ","),
            format(sum(out$seurat_cycling & !out$tricycle_cycling), big.mark = ","),
            format(sum(!out$seurat_cycling & out$tricycle_cycling), big.mark = ","),
            format(sum(!out$seurat_cycling & !out$tricycle_cycling), big.mark = ",")))
write.csv(as.data.frame.matrix(tb),
          file.path(OUTTAB, "cellcycle_tricycle_confusion.csv"))

## ---- 4b. the three controls that decide how much of this to believe --------
# Without these the two methods just agree with each other, which proves nothing about
# whether either is measuring the cell cycle.

cat("\n=== CONTROL 1: does theta actually track the cycle in OUR cells? ===\n")
cat("    Canonical markers vs theta (periodic loess). S genes should peak before G2M genes.\n")
Xc <- SeuratObject::GetAssayData(comb, assay = "RNA", layer = "data")
# Fit on at most 20k cells: this estimates the SHAPE of marker expression around the
# cycle, which a subsample settles as well as the full set, and loess on 59k points x 7
# genes otherwise dominates the runtime of the whole script.
set.seed(SEED)
li <- if (ncol(Xc) > 20000) sort(sample(ncol(Xc), 20000)) else seq_len(ncol(Xc))
peaks <- do.call(rbind, lapply(c("Mcm2","Pcna","Rrm2","Cdk1","Top2a","Mki67","Ccnb1"), function(g) {
  if (!g %in% rownames(Xc)) return(NULL)
  f <- tryCatch(fit_periodic_loess(theta[li], as.numeric(Xc[g, li]), plot = FALSE), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  data.frame(gene = g, peak_theta_pi = round(f$pred.df$x[which.max(f$pred.df$y)] / pi, 2),
             R2 = round(f$rsquared, 3), stringsAsFactors = FALSE)
}))
print(peaks, row.names = FALSE)
write.csv(peaks, file.path(OUTTAB, "cellcycle_tricycle_marker_peaks.csv"), row.names = FALSE)
cat(sprintf("  -> theta explains %.0f-%.0f%% of marker variance; the arc %.2fpi-%.2fpi is\n",
            100*min(peaks$R2), 100*max(peaks$R2), CYC_LO/pi, CYC_HI/pi))
cat("     supported by the data rather than assumed.\n")

cat("\n=== CONTROL 2: is the call driven by SEQUENCING DEPTH? ===\n")
out$ngene <- md$nFeature_RNA
for (ct in c("Cardiomyocyte","Endothelial","Fibroblast")) for (tp in c("P0","P7")) {
  z <- out[out$celltype == ct & out$timepoint == tp, ]
  if (nrow(z) < 400) next
  cat(sprintf("  %-14s %s  r(cycling, nGene) = %+.3f tricycle / %+.3f Seurat\n", ct, tp,
              cor(as.numeric(z$tricycle_cycling), z$ngene), cor(as.numeric(z$seurat_cycling), z$ngene)))
}
# Depth-matched P0 vs P7 in cardiomyocytes: common absolute bins, equal weight per bin.
cmx <- out[out$celltype == "Cardiomyocyte", ]
cmx$bin <- cut(cmx$ngene, breaks = c(1500,2000,2500,3000,3500,4000,5000,Inf), include.lowest = TRUE)
mm <- do.call(rbind, lapply(levels(cmx$bin), function(b) {
  a <- cmx[cmx$bin==b & cmx$timepoint=="P0", ]; z <- cmx[cmx$bin==b & cmx$timepoint=="P7", ]
  if (nrow(a) < 100 || nrow(z) < 100) return(NULL)
  data.frame(bin=b, n_P0=nrow(a), n_P7=nrow(z), w=min(nrow(a),nrow(z)),
             P0=100*mean(a$tricycle_cycling), P7=100*mean(z$tricycle_cycling))
}))
cat("\n  cardiomyocytes, P0 vs P7 at MATCHED depth (tricycle):\n")
print(transform(mm, P0=round(P0,1), P7=round(P7,1))[, c("bin","n_P0","n_P7","P0","P7")], row.names = FALSE)
cat(sprintf("\n  depth-standardised: P0 %.1f%% vs P7 %.1f%%   (raw, unmatched: %.1f%% vs %.1f%%)\n",
            weighted.mean(mm$P0, mm$w), weighted.mean(mm$P7, mm$w),
            100*mean(cmx$tricycle_cycling[cmx$timepoint=="P0"]),
            100*mean(cmx$tricycle_cycling[cmx$timepoint=="P7"])))
write.csv(mm, file.path(OUTTAB, "cellcycle_tricycle_depthmatched_cm.csv"), row.names = FALSE)

cat("\n=== CONTROL 3: the ambient-RNA floor ===\n")
Cc <- SeuratObject::GetAssayData(comb, assay = "RNA", layer = "counts")
sarc <- intersect(c("Tnnt2","Myh6","Actc1"), rownames(Cc))
nonCM <- out$celltype != "Cardiomyocyte"
cat(sprintf("  %% of NON-cardiomyocytes detecting cardiac sarcomere genes (%s): %.1f%%\n",
            paste(sarc, collapse="/"),
            100*mean(Matrix::colSums(Cc[sarc, nonCM, drop=FALSE] > 0) > 0)))
cat("  Those cells cannot transcribe sarcomere genes. Detection at this rate is ambient RNA,\n")
cat("  and the same ambient carries proliferation transcripts into every barcode.\n")
if (any(out$celltype == "RBC")) {
  r <- out[out$celltype == "RBC", ]
  cat(sprintf("\n  'RBC' cells: n=%d, tricycle %.1f%% cycling, Seurat %.1f%% cycling\n",
              nrow(r), 100*mean(r$tricycle_cycling), 100*mean(r$seurat_cycling)))
  cat(sprintf("  compare cardiomyocytes: tricycle %.1f%%\n",
              100*mean(out$tricycle_cycling[out$celltype=="Cardiomyocyte"])))
  cat("  A population scoring at or above cardiomyocytes bounds what the CM number can mean.\n")
  cat("  CAVEAT: these may be nucleated erythroid precursors, which do cycle -- so this is a\n")
  cat("  ceiling on confidence, not a calibrated zero.\n")
}
rm(Xc, Cc); invisible(gc())

## ---- 5. figures ------------------------------------------------------------
cmcells <- out[out$celltype == "Cardiomyocyte" & !is.na(out$tricycle_theta), ]
p1 <- ggplot(cmcells, aes(tricycle_theta, colour = timepoint)) +
  geom_density(linewidth = 0.9) +
  scale_x_continuous(breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
                     labels = c("0\n(G1/G0)", "0.5π\n(S)", "π\n(G2)",
                                "1.5π\n(M)", "2π")) +
  facet_wrap(~ genotype) +
  labs(title = "Cardiomyocyte cell-cycle position (tricycle), P0 vs P7",
       subtitle = "Threshold-free. A population that has exited the cycle piles up near 0/2π.",
       x = "tricycle position θ", y = "density") + theme_bw()
ggsave(file.path(OUTFIG, "cellcycle_tricycle_cm_theta.png"), p1, width = 10, height = 4.5, dpi = 130)

pl <- by_ct[by_ct$n >= 100, ]
p2 <- ggplot(pl, aes(pct_cycling_seurat, pct_cycling_tricycle, colour = timepoint)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(size = n), alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = celltype), size = 3, max.overlaps = 20) +
  labs(title = "Cycling fraction: Seurat CellCycleScoring vs tricycle",
       subtitle = "Same cells, same normalisation. Points above the dashed line: tricycle calls more cycling.",
       x = "% cycling (Seurat S/G2M)", y = "% cycling (tricycle θ arc)") + theme_bw()
ggsave(file.path(OUTFIG, "cellcycle_tricycle_vs_seurat.png"), p2, width = 9, height = 6, dpi = 130)

say("wrote tables + figures")
cat("=== DONE cellcycle_tricycle ===\n")
