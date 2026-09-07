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
