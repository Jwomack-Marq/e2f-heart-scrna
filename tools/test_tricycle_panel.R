#!/usr/bin/env Rscript
# test_tricycle_panel.R
# ---------------------------------------------------------------------------
# The Cell cycle (tricycle) tab must draw the cycle it claims to draw.
#
# What this guards, and why none of it is obvious from reading the panel:
#
#   1. THE ANGLE IS THE MEASUREMENT. The wheel places each cell at
#      (r·cos θ, r·sin θ) and then annotates fixed G1/S/G2/M anchors at 0, π/2, π, 3π/2.
#      That is only correct because tricycle's θ IS atan2(PC2, PC1) of the same embedding.
#      If a future run ever stored θ against a different centring, the plot would still
#      render, still look like a cycle, and put every phase label in the wrong place.
#   2. PERCENTAGES MUST NOT DEPEND ON THE DISPLAY. The wheel thins to `Cells drawn` for
#      speed while the subtitle and notes quote fractions. Those must come from the full
#      slice, so this checks the reported number against the unthinned data.
#   3. WEDGE LENGTH IS A SHARE, NOT A COUNT. The four groups differ in size by ~2x
#      (P7 KO 10,597 vs P7 WT 6,537). On raw counts the knockout's wedges would be longer
#      before any biology entered -- in the one figure built to compare them.
#   4. ABSTENTIONS SURVIVE. 54% of cells are ones tricycle declined to stage. NA would be
#      silently dropped by most of the plotting path; they must reach the app as a level.
#   5. THE PROSE READS THE RUN. κ, concordance and the ambient floor are quoted in the
#      panel and must come from the controls table, not from a remembered number.
#
#   docker run --rm --network none -v "$PWD:/repo" -u "$(id -u):$(id -g)" -e HOME=/tmp \
#     lab-server-e2f-heart-scrna-dev:latest Rscript /repo/tools/test_tricycle_panel.R
# ---------------------------------------------------------------------------
library(shiny)
.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
setwd(normalizePath(file.path(dirname(.this), "..", "shiny_app")))
FAIL <- 0; ok <- function(l,c){ if(!isTRUE(c)) FAIL<<-FAIL+1
  cat(sprintf("  [%s] %s\n", if(isTRUE(c))"PASS" else "FAIL", l)) }

E <- new.env(parent = globalenv()); suppressWarnings(suppressMessages(sys.source("app.R", envir = E)))
TRI <- E$TRI
if (is.null(TRI)) { cat("  [SKIP] no app$tricycle in this bundle -- run build_tricycle.R\n"); quit(status = 0) }
pc <- TRI$percell

cat("\n== the bundle carries what the tab needs ==\n")
ok(sprintf("per-cell table present (%s cells)", format(nrow(pc), big.mark = ",")), nrow(pc) > 50000)
ok("cycle-space coordinates kept, not just the angle",
   all(c("tricycle_pc1","tricycle_pc2") %in% names(pc)) && !all(is.na(pc$tricycle_pc1)))
ok("depth carried for the confound panel", "ngene" %in% names(pc) && !all(is.na(pc$ngene)))
ok("controls table present", !is.null(TRI$controls) && nrow(TRI$controls) > 10)

cat("\n== 1. the phase anchors sit where the angle says ==\n")
# The whole wheel hangs on this identity. Tolerance is for the 4-dp rounding in the CSV.
th <- atan2(pc$tricycle_pc2, pc$tricycle_pc1); th[th < 0] <- th[th < 0] + 2*pi
dif <- abs(th - pc$tricycle_theta); dif <- pmin(dif, 2*pi - dif)
cat(sprintf("   max circular |atan2(PC2,PC1) - theta| = %.5f rad\n", max(dif, na.rm = TRUE)))
ok("theta IS the angle of the stored embedding", max(dif, na.rm = TRUE) < 0.01)
# and the markers must run S-then-G2/M around it, which is what licenses the cycling arc
pk <- TRI$marker_peaks
sp <- mean(pk$peak_theta_pi[pk$gene %in% c("Mcm2","Pcna","Rrm2")])
gp <- mean(pk$peak_theta_pi[pk$gene %in% c("Cdk1","Top2a","Mki67","Ccnb1")])
cat(sprintf("   mean peak: S genes %.2fpi, G2/M genes %.2fpi\n", sp, gp))
ok("S markers peak BEFORE G2/M markers along theta", sp < gp)
ok("both fall inside the 0.5pi-1.5pi arc called cycling", sp > 0.5 && gp < 1.5)

cat("\n== 4. abstentions reach the app as a level, not as NA ==\n")
ok("'(not staged)' is a stage level", TRI$not_staged %in% levels(pc$tricycle_stage))
ok("no NA stages left", !any(is.na(pc$tricycle_stage)))
cat(sprintf("   abstained: %.1f%%\n", 100*mean(pc$tricycle_stage == TRI$not_staged)))

cat("\n== 3. the rose reports shares, so unequal groups compare ==\n")
d <- E$tri_slice("Cardiomyocyte")
# fill_by is deliberately ALSO a facet variable here: that is the KO-vs-WT view, and it is
# the case where a duplicated key silently doubled every share.
g  <- E$tri_rose_gg(d, facet = c("timepoint","genotype"), fill_by = "genotype")
tb <- attr(g, "rose_data")
tot <- tapply(tb$share, interaction(tb$timepoint, tb$genotype), sum)
cat(sprintf("   per-panel share totals: %s\n", paste(sprintf("%.3f", tot), collapse = ", ")))
ok("every panel sums to 1 (shares, not counts)", all(abs(tot - 1) < 1e-9))
# and again with a fill that is NOT a facet variable, so both paths are covered
tb2 <- attr(E$tri_rose_gg(d, facet = c("timepoint","genotype")), "rose_data")
ok("same when the fill is not a facet variable",
   all(abs(tapply(tb2$share, interaction(tb2$timepoint, tb2$genotype), sum) - 1) < 1e-9))
# and the sizes really are unequal, so the normalisation is doing work rather than being moot
n <- table(d$timepoint, d$genotype)
ok(sprintf("group sizes genuinely differ (max/min = %.2fx)", max(n)/min(n)), max(n)/min(n) > 1.3)

cat("\n== 2. quoted fractions come from all cells, not the thinned draw ==\n")
full <- E$tri_slice("Cardiomyocyte", "P7", "KO")
thin <- E$tri_slice("Cardiomyocyte", "P7", "KO", maxn = 2000L)
ok("thinning actually thins", nrow(thin) == 2000L && nrow(full) > 2000L)
ref <- TRI$by_celltype
r <- ref[ref$celltype == "Cardiomyocyte" & ref$timepoint == "P7" & ref$genotype == "KO", ]
cat(sprintf("   full slice %.1f%% cycling vs summary table %.1f%%\n",
            100*mean(full$tricycle_cycling), r$pct_cycling_tricycle[1]))
ok("full slice reproduces the analysis' own summary",
   abs(100*mean(full$tricycle_cycling) - r$pct_cycling_tricycle[1]) < 0.15)

cat("\n== 5. the prose reads the run ==\n")
k <- E$tri_ctrl("kappa"); cc <- E$tri_ctrl("pct_concordance")
amb <- E$tri_ctrl("pct_nonCM_detecting_sarcomere", "Tnnt2/Myh6/Actc1")
cat(sprintf("   kappa %.3f | concordance %.1f%% | ambient %.1f%%\n", k, cc, amb))
ok("kappa read from the controls table", !is.na(k) && k > 0.5)
ok("concordance read from the controls table", !is.na(cc) && cc > 50)
ok("ambient floor read from the controls table", !is.na(amb) && amb > 50)

cat("\n== every figure in the panel builds ==\n")
mk <- function(l, expr) ok(l, { p <- try(expr, silent = TRUE)
  !inherits(p, "try-error") && inherits(p, "ggplot") &&
    !inherits(try(ggplot2::ggplot_build(p), silent = TRUE), "try-error") })
mk("wheel — by stage",      E$tri_wheel_gg(d, "tricycle_stage", show_peaks = TRUE))
mk("wheel — by theta (cyclic scale)", E$tri_wheel_gg(d, "tricycle_theta"))
mk("wheel — raw radius",    E$tri_wheel_gg(d, "tricycle_stage", rscale = "raw"))
mk("wheel — faceted WT/KO", E$tri_wheel_gg(d, "tricycle_stage", facet = c("timepoint","genotype")))
mk("rose — by stage",       E$tri_rose_gg(d))
mk("theta density",         E$tri_theta_density_gg(d))
mk("stage composition",     E$tri_comp_gg(d))
mk("tricycle vs Seurat",    E$tri_vs_gg(TRI$by_celltype))
mk("confusion matrix",      E$tri_confusion_gg(TRI$confusion))
mk("marker peaks",          E$tri_peaks_gg(TRI$marker_peaks))
mk("depth quartiles",       E$tri_depth_gg(E$tri_depth_df(c("Cardiomyocyte","Endothelial"))))

cat("\n== the depth-matched KO-vs-WT comparison ==\n")
M <- TRI$matched
if (is.null(M)) {
  cat("  [SKIP] no depth-matched table -- run cellcycle_tricycle_depthmatched.R\n")
} else {
  # The matching's own success criterion. Without this the raw-vs-matched arrows are just
  # two numbers; 0.5 is what says genotype can no longer be told from depth at all.
  cat(sprintf("   depth AUC after matching: %s\n",
              paste(sprintf("%.3f", M$auc_depth_matched), collapse = ", ")))
  ok("matching drove depth AUC to 0.5 in every stratum",
     all(abs(M$auc_depth_matched - 0.5) < 0.01))
  ok("and it was NOT already 0.5 before (so the match did work)",
     any(abs(M$auc_depth_raw - 0.5) > 0.05))
  cm7 <- M[M$celltype == "Cardiomyocyte" & M$timepoint == "P7", ]
  cat(sprintf("   P7 cardiomyocytes: KO-WT %+.1f raw -> %+.1f matched (WT was the deeper library: %.0f vs %.0f UMIs)\n",
              cm7$gap_raw, cm7$gap_matched, cm7$median_numi_WT_raw, cm7$median_numi_KO_raw))
  # The finding worth pinning: at P7 the WT is deeper, depth inflates the cycling call, so
  # the raw gap understates rather than exaggerates. If a rebuild ever flips this, the
  # panel's prose ("raw was understating it") becomes a lie.
  ok("at P7 the WT cardiomyocyte library is the deeper one",
     cm7$median_numi_WT_raw > cm7$median_numi_KO_raw)
  ok("so the P7 KO-WT gap does not shrink at matched depth",
     cm7$gap_matched >= cm7$gap_raw - 0.05)
  ok("gap stays positive (KO cycles more at P7)", cm7$gap_matched > 0)
  mk("depth-matched gap figure", E$tri_matched_gg(M))

  # The panel's prose and chapter 12 both assert the cross-cell-type SHAPE, not just the
  # cardiomyocyte number: no genotype difference at P0, and a graded rather than uniform
  # spread at P7. That spread is the argument that a genotype-wide sort or library artifact
  # is not what produces this, so if a rebuild flattened it the prose would be a lie.
  z  <- M[M$celltype != "RBC", ]                 # ambient-floor population, excluded
  p0 <- z[z$timepoint == "P0", ]; p7 <- z[z$timepoint == "P7", ]
  cat(sprintf("   P0 gaps span %+.1f..%+.1f; P7 gaps span %+.1f..%+.1f (%.1f points)\n",
              min(p0$gap_matched), max(p0$gap_matched),
              min(p7$gap_matched), max(p7$gap_matched),
              max(p7$gap_matched) - min(p7$gap_matched)))
  ok("no genotype cycling difference at P0 (all within +/-4 points)",
     all(abs(p0$gap_matched) < 4))
  ok("the P7 effect is graded, not a uniform shift (spread > 10 points)",
     max(p7$gap_matched) - min(p7$gap_matched) > 10)
  ok("fibroblasts carry the largest P7 gap, as the prose says",
     p7$celltype[which.max(p7$gap_matched)] == "Fibroblast")
  ok("and at least one cell type reverses sign at P7", any(p7$gap_matched < 0))
}

cat("\n== the depth confound is visible, which is why the panel exists ==\n")
dd <- E$tri_depth_df("Endothelial")
p0 <- dd[dd$timepoint == "P0", ]
p0 <- p0[order(p0$quartile), ]
cat(sprintf("   endothelial P0 across depth quartiles: %s\n",
            paste0(p0$pct_cycling_tricycle, "%", collapse = " -> ")))
ok("cycling fraction rises steeply with depth (>2x across quartiles)",
   max(p0$pct_cycling_tricycle) > 2 * min(p0$pct_cycling_tricycle))

cat(sprintf("\n%s\n", if (FAIL == 0) "ALL PASS" else paste(FAIL, "FAILURES")))
if (FAIL) quit(status = 1)
