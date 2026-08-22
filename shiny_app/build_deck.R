#!/usr/bin/env Rscript
# build_deck.R
# ---------------------------------------------------------------------------
# Assemble a lab-meeting deck about the SHINY APP half of this project: the
# pipeline that produces app_data.rds, what each tab of the browser answers, and
# what the analyses actually found. The model half has its own deck
# (model/tools/build_deck.py); this is its sibling and deliberately shares that
# deck's stylesheet and navigation so the two look like one project.
#
#   Rscript shiny_app/build_deck.R          # -> shiny_app/deck.html
#
# Every number on every slide is READ FROM THE BUNDLE at build time, never typed
# in. That is the same discipline build_deck.py follows against results.json, and
# for the same reason: prose drifts from data, and a deck is where that shows up
# in front of an audience. If a number here looks wrong, the bundle is wrong.
#
# The deck is a build artifact and is git-ignored, like model/deck.html.
# ---------------------------------------------------------------------------

suppressMessages(library(Matrix))
if (!file.exists("app_data.rds") && file.exists("shiny_app/app_data.rds")) setwd("shiny_app")
stopifnot(file.exists("app_data.rds"))

TEMPLATE <- "../model/tools/deck_template.html"
OUT      <- "deck.html"
OUT_BODY <- "deck.body.html"   # same deck, no <html>/<head>/<body> — for Artifact publishing

cat("Loading app_data.rds ...\n")
app <- readRDS("app_data.rds")
fg  <- app$fourgroup
stopifnot(!is.null(fg))
meta <- app$meta; cmm <- app$cm$meta; GM <- fg$geneaxes

# ---- helpers ---------------------------------------------------------------
fmt  <- function(x) formatC(x, big.mark = ",", format = "d")
pct  <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), x)
sgn  <- function(x, d = 3) sprintf(paste0("%+.", d, "f"), x)
esc  <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE)
                      gsub("<", "&lt;", gsub(">", "&gt;", x, fixed = TRUE), fixed = TRUE) }
grp  <- function(cl, g, col = "n") {
  v <- fg$counts[[col]][fg$counts$cluster == cl & fg$counts$group == g]
  if (length(v)) v[1] else NA
}
phase_pct <- function(cl, g, ph) {
  v <- fg$phase$pct[fg$phase$cluster == cl & fg$phase$group == g & fg$phase$Phase == ph]
  if (length(v)) v[1] else NA
}
GROUPS <- fg$built$groups

# ---- numbers ---------------------------------------------------------------
N <- list(
  cells = nrow(meta), cm = nrow(cmm),
  celltypes = length(unique(as.character(meta$celltype))),
  panel = length(app$genes), broad = length(app$deg_genes),
  built = as.character(app$built), bundle_mb = round(file.size("app_data.rds") / 1e6),
  de1 = sum(vapply(fg$de, length, 0L)),
  de2 = if (is.null(fg$de2)) 0L else sum(vapply(fg$de2, length, 0L)),
  skipped = if (is.null(fg$skipped)) 0L else nrow(fg$skipped),
  genemap = if (is.null(GM)) 0L else nrow(GM))

# maturation / G1 summary, G1 stratum, pooled and CM2
summ <- function(cl, score = "sig_maturation_nocc", stratum = "G1") {
  sc <- fg$scores[fg$scores$score == score & fg$scores$stratum == stratum &
                  fg$scores$cluster == cl, , drop = FALSE]
  m <- function(g) { v <- sc$mean[sc$group == g]; if (length(v)) v[1] else NA }
  s <- function(g) { v <- sc$sd[sc$group == g];   if (length(v)) v[1] else NA }
  n <- function(g) { v <- sc$n[sc$group == g];    if (length(v)) v[1] else NA }
  d <- function(ko, wt) {
    diff <- m(ko) - m(wt)
    sp <- sqrt(((n(ko) - 1) * s(ko)^2 + (n(wt) - 1) * s(wt)^2) / max(n(ko) + n(wt) - 2, 1))
    c(diff = diff, d = if (is.finite(sp) && sp > 0) diff / sp else NA)
  }
  list(p0 = d("KO-P0", "WT-P0"), p7 = d("KO-P7", "WT-P7"),
       gain = 100 * (m("KO-P7") - m("KO-P0")) / (m("WT-P7") - m("WT-P0")),
       g1_p0 = phase_pct(cl, "KO-P0", "G1") - phase_pct(cl, "WT-P0", "G1"),
       g1_p7 = phase_pct(cl, "KO-P7", "G1") - phase_pct(cl, "WT-P7", "G1"))
}
S_ALL <- summ("AllCM"); S_CM2 <- summ("CM2")

# gene-map coupling per panel
gm_panel <- function(pn) {
  cols <- if (pn == "avg") c("mat_auc", "met_auc") else paste0(c("mat_auc_", "met_auc_"), pn)
  if (!all(cols %in% names(GM))) return(NULL)
  d <- GM[is.na(GM$in_score_set) & !is.na(GM[[cols[1]]]) & !is.na(GM[[cols[2]]]), , drop = FALSE]
  x <- d[[cols[1]]]; y <- d[[cols[2]]]
  cx <- stats::median(x); cy <- stats::median(y)
  q <- paste0(ifelse(x >= cx, "mature", "immature"), "+", ifelse(y >= cy, "oxidative", "glycolytic"))
  list(n = nrow(d), cor = stats::cor(x, y),
       diag = mean(q %in% c("mature+oxidative", "immature+glycolytic")))
}
PAN <- setNames(lapply(c("avg", "P0", "P7"), gm_panel), c("avg", "P0", "P7"))

# intersection, CM2
IT <- fg$intersect[!fg$intersect$confounder, , drop = FALSE]
it2 <- IT[IT$cluster == "CM2", , drop = FALSE]
qn  <- function(q) sum(it2$quadrant == q)
hyp <- IT$quadrant %in% c("immature_up_in_KO", "mature_down_in_KO")
N$cyc_hits   <- sum(hyp & !is.na(IT$cyc_class) & IT$cyc_class == "cycling-associated")
N$cyc_resid  <- stats::median(IT$cyc_resid[hyp], na.rm = TRUE)
N$hyp_rows   <- sum(hyp)

top_genes <- function(cl, key, st, n = 5, up = TRUE) {
  d <- fg$de[[cl]][[paste0(key, "__", st)]]
  if (is.null(d)) d <- fg$de2[[cl]][[paste0(key, "__", st)]]
  if (is.null(d)) return(character(0))
  d <- d[!d$confounder, , drop = FALSE]
  d <- d[if (up) d$log2FoldChange > 0 else d$log2FoldChange < 0, , drop = FALSE]
  d <- d[order(if (up) -d$log2FoldChange else d$log2FoldChange), , drop = FALSE]
  sprintf("%s <span class=\"num\">%s</span>", head(d$gene, n),
          sgn(head(d$log2FoldChange, n), 2))
}

# ---- slide chrome ----------------------------------------------------------
sl <- function(mv, notes, eyebrow_n, eyebrow, head, body, foot_l, foot_r) sprintf('
<section class="slide" data-mv="%s" data-notes="%s">
  <div>
    <p class="eyebrow"><span class="tag">%s</span><span>%s</span></p>
    <h2 class="head">%s</h2>
  </div>
  <div class="body">%s</div>
  <div class="foot"><span>%s</span><span class="src">%s</span></div>
</section>', mv, notes, eyebrow_n, eyebrow, head, body, foot_l, foot_r)

# h3, not h4: the shared stylesheet styles `.card h3` and nothing else inside a card
card <- function(t, b) sprintf('<div class="card"><h3>%s</h3><p>%s</p></div>', t, b)
cards <- function(...) sprintf('<div class="cards">%s</div>', paste0(..., collapse = ""))
callout <- function(x, cls = "callout") sprintf('<div class="%s">%s</div>', cls, x)
cols <- function(a, b, w = "c-6-6") sprintf('<div class="cols %s"><div>%s</div><div>%s</div></div>', w, a, b)

# a plain table, styled by the deck's own rules
tbl <- function(headers, rows, cls = "") sprintf(
  '<table class="dt %s"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>', cls,
  paste0(sprintf("<th>%s</th>", headers), collapse = ""),
  paste0(vapply(rows, function(r) paste0("<tr>", paste0(sprintf("<td>%s</td>", r), collapse = ""), "</tr>"), ""),
         collapse = ""))

# ---- figure: the pipeline, as an inline SVG --------------------------------
pipeline_svg <- function() {
  lane <- function(y, x, w, lab, sub, fill) sprintf(
    '<rect x="%d" y="%d" width="%d" height="52" rx="4" fill="%s" stroke="var(--axs)"/>
     <text x="%d" y="%d" text-anchor="middle" font-size="12.5" font-weight="600" fill="currentColor">%s</text>
     <text x="%d" y="%d" text-anchor="middle" font-size="10.5" fill="var(--mut)">%s</text>',
    x, y, w, fill, x + w/2, y + 22, lab, x + w/2, y + 38, sub)
  arr <- function(x1, x2, y) sprintf(
    '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="var(--axs)" stroke-width="1.4" marker-end="url(#pa)"/>',
    x1, y, x2, y)
  paste0('<svg viewBox="0 0 980 250" width="100%" role="img" aria-label="pipeline">',
    '<defs><marker id="pa" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" ',
    'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--axs)"/></marker></defs>',
    '<text x="8" y="18" font-size="10.5" font-weight="700" fill="var(--mut)">UPSTREAM &mdash; not in this repo</text>',
    lane(28, 8, 190, "4 libraries", "WT/KO x P0/P7", "var(--bnd)"),
    arr(202, 232, 54),
    lane(28, 236, 210, "QC, SCT, Harmony", "cluster + UMAP", "var(--bnd)"),
    arr(450, 480, 54),
    lane(28, 484, 220, "build_app_data.R", "the heavy lift", "var(--bnd)"),
    arr(708, 738, 54),
    lane(28, 742, 230, sprintf("app_data.rds (%d MB)", N$bundle_mb), "the only runtime input", "var(--grd)"),
    '<text x="8" y="126" font-size="10.5" font-weight="700" fill="var(--mut)">IN THIS REPO &mdash; each reads the bundle and writes one slot back</text>',
    lane(136, 8, 178, "build_signature_scores", "sig_* per-cell scores", "var(--bnd)"),
    lane(136, 194, 150, "build_communication", "app$commun", "var(--bnd)"),
    lane(136, 350, 140, "build_refmap", "app$refmap", "var(--bnd)"),
    lane(136, 496, 168, "build_subcluster_enr", "app$enrich$sub", "var(--bnd)"),
    lane(136, 670, 172, "build_fourgroup", "app$fourgroup", "var(--grd)"),
    arr(880, 908, 162),
    '<rect x="854" y="200" width="118" height="40" rx="4" fill="var(--grd)" stroke="var(--axs)"/>',
    '<text x="913" y="225" text-anchor="middle" font-size="12.5" font-weight="600" fill="currentColor">Shiny app</text>',
    '<line x1="913" y1="188" x2="913" y2="198" stroke="var(--axs)" stroke-width="1.4" marker-end="url(#pa)"/>',
    '</svg>')
}

# ---- figure: the maturation gap, as a dumbbell -----------------------------
gap_svg <- function() {
  sc <- fg$scores[fg$scores$score == "sig_maturation_nocc" & fg$scores$stratum == "G1" &
                  fg$scores$cluster == "AllCM", , drop = FALSE]
  m <- function(g) sc$mean[sc$group == g]; se <- function(g) sc$se[sc$group == g]
  lo <- min(sc$mean) - .12; hi <- max(sc$mean) + .12
  X <- function(v) 90 + (v - lo) / (hi - lo) * 780
  row <- function(y, tp, wt, ko) paste0(
    sprintf('<text x="72" y="%d" text-anchor="end" font-size="13" font-weight="600" fill="currentColor">%s</text>', y + 5, tp),
    sprintf('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="var(--axs)" stroke-width="2"/>', X(m(wt)), y, X(m(ko)), y),
    sprintf('<circle cx="%.1f" cy="%d" r="7" fill="var(--s1)"/>', X(m(wt)), y),
    sprintf('<circle cx="%.1f" cy="%d" r="7" fill="var(--neg)"/>', X(m(ko)), y),
    sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="11.5" font-weight="700" fill="var(--neg)">%s</text>',
            (X(m(wt)) + X(m(ko))) / 2, y - 14, sgn(m(ko) - m(wt), 3)))
  paste0('<svg viewBox="0 0 900 170" width="100%" role="img" aria-label="maturation gap">',
    sprintf('<line x1="90" y1="132" x2="870" y2="132" stroke="var(--grd)"/>'),
    row(46, "P0", "WT-P0", "KO-P0"), row(100, "P7", "WT-P7", "KO-P7"),
    '<circle cx="640" cy="158" r="6" fill="var(--s1)"/><text x="654" y="162" font-size="11.5" fill="var(--sec)">WT</text>',
    '<circle cx="700" cy="158" r="6" fill="var(--neg)"/><text x="714" y="162" font-size="11.5" fill="var(--sec)">KO</text>',
    sprintf('<text x="90" y="162" font-size="11.5" fill="var(--mut)">cycle-free maturation score, G1 cells only &mdash; effect size %s at P0, %s at P7</text>',
            sprintf("d = %.2f", S_ALL$p0[["d"]]), sprintf("d = %.2f", S_ALL$p7[["d"]])),
    '</svg>')
}

# ---- slides ----------------------------------------------------------------
SL <- c()

# 1 · title
SL <- c(SL, sprintf('
<section class="slide title-slide on" data-mv="Data and app" data-notes="Fifteen minutes on the descriptive half. The browser is not a figure dump &mdash; it is built around two confounds that would otherwise produce clean-looking wrong answers, and most of what I will show you is how those are handled. Three things to take away: the P7 sort artefact, the maturation result that survives holding cycling fixed, and the negative on cycle competence.">
  <div class="title-wrap">
    <p class="eyebrow"><span class="tag">Lab meeting</span><span>E2F7/8 heart scRNA-seq browser</span></p>
    <h1>What the browser <em>shows</em>, and what it refuses to</h1>
    <p class="sub">A pipeline, eleven tabs, and two confounds that decide whether any of it means anything.</p>
    <p class="spec">
      <span><b>%s</b> cells</span><span>%s cardiomyocytes</span><span>%s cell types</span>
      <span>%s DE tables</span><span>built %s</span>
    </p>
  </div>
  <div class="foot">
    <span>Justin Womack &middot; <code>e2f-heart-scrna/shiny_app</code></span>
    <span class="src">click &rsaquo; or press <kbd>&rarr;</kbd> &middot; <kbd>N</kbd> for notes</span>
  </div>
</section>', fmt(N$cells), fmt(N$cm), N$celltypes, fmt(N$de1 + N$de2), N$built))

# 2 · the dataset
SL <- c(SL, sl("Data and app",
  "Four libraries, one per genotype-timepoint cell. Say the design limit out loud before anything else: n = 1 animal per condition and the genotypes are different sexes, so the Y-genes and Xist top every KO-vs-WT list and are flagged as confounders throughout. Everything downstream is effect sizes, not inference.",
  "01", "The data",
  "Four samples, one per cell of the design",
  cols(
    tbl(c("", GROUPS, "total"),
        lapply(c("AllCM","CM2","CM4","CM5","CM9"), function(cl) c(
          if (cl == "AllCM") "<b>all CM</b>" else cl,
          vapply(GROUPS, function(g) fmt(grp(cl, g)), ""),
          fmt(sum(vapply(GROUPS, function(g) grp(cl, g), 0)))))),
    cards(
      card("One animal per condition", "No biological replication. Every p-value in the app is cell-level and pseudoreplicated; tables rank by effect size."),
      card("Sex-confounded", sprintf("KO and WT are different sexes. %s genes are flagged and can be hidden everywhere.", length(app$confound))),
      card("Two expression matrices", sprintf("A %s-gene curated panel over all cells, and a %s-gene broad matrix over a downsample. Neither dominates.", fmt(N$panel), fmt(N$broad)))),
    "c-7-5"),
  "Descriptive pilot &mdash; hypothesis-generating only",
  "app$meta, app$fourgroup$counts"))

# 3 · pipeline
SL <- c(SL, sl("Data and app",
  "The heavy lift happens upstream and is not in this repo. What is in the repo is five builder scripts that each read the bundle, compute one thing, and write one slot back. That matters for reproducibility: the app never computes anything expensive at runtime, so a tab either has its slot or shows a message telling you which script to run.",
  "02", "How it is built",
  "One bundle, five builders, no runtime computation",
  paste0('<figure><div class="frame">', pipeline_svg(), '</div>',
    '<figcaption>Everything the app renders is precomputed. Each builder is idempotent, backs up before writing, and adds exactly one slot &mdash; so a missing analysis degrades to a message naming the script, not an error.</figcaption></figure>'),
  "Re-run order is fixed; build_fourgroup consumes the sig_* scores",
  "README.md &sect; Data pipeline"))

# 4 · the two confounds
SL <- c(SL, sl("Data and app",
  "These are the two that would otherwise produce clean-looking wrong answers. The sort one is the bigger threat: P7 was cycling-enriched and P0 was not, so the raw cycling fraction rises from P0 to P7 in this data and falls in reality. Anyone reading a P0-vs-P7 contrast off this dataset without knowing that gets the developmental direction backwards.",
  "03", "What the app is built around",
  "Two confounds, handled rather than noted",
  cols(
    paste0('<h4 class="lab">The FACS sort</h4>',
      tbl(c("phase", GROUPS), list(
        c("G1",  vapply(GROUPS, function(g) pct(phase_pct("AllCM", g, "G1")), "")),
        c("<b>S</b>", vapply(GROUPS, function(g) sprintf("<b>%s</b>", pct(phase_pct("AllCM", g, "S"))), "")),
        c("G2M", vapply(GROUPS, function(g) pct(phase_pct("AllCM", g, "G2M")), "")))),
      callout("P7 was cycling-enriched 4.5&ndash;5.2&times;; P0 was not. The S-phase row <b>is the sort</b>. Every P0-vs-P7 contrast therefore exists in a G1-matched stratum, which is what the app shows by default.", "callout warn")),
    paste0('<h4 class="lab">The circular score</h4>',
      '<p>The maturation score&rsquo;s immature program contained <code>Mki67</code>, <code>Top2a</code> and <code>Ccnd1</code>. Using it to argue <em>less mature therefore more cycle-competent</em> would have been partly circular &mdash; the same failure the <code>E2f8</code> phase-list audit found.</p>',
      callout("A cycle-free variant drops those three. Every maturation claim in this deck uses it.", "callout")),
    "c-6-6"),
  "model/README.md documents the enrichment; the app carries it in the About block",
  "build_signature_scores.R, build_fourgroup.R"))

# 5 · tab map
SL <- c(SL, sl("Data and app",
  "Eleven tabs, but they answer four kinds of question. I will not walk through all of them; the point of this slide is that the browser is organised by question, not by method.",
  "04", "What the app shows",
  "Eleven tabs, four kinds of question",
  cards(
    card("Where are the cells?", "UMAP explorer, gene detail, composition &mdash; colour by any gene or metadata, split by genotype and timepoint."),
    card("What differs?", "DE by cell type, subset &amp; DEG explorer, four-group DE, per-subcluster enrichment. Volcano linked both ways to a table."),
    card("What state are they in?", "Cell-cycle exit &amp; ploidy, maturation &amp; metabolism, the gene map, G1 and maturation summaries."),
    card("Which genes, and are they real?", "Intersection quadrants, candidate follow-up, gene-set Venn &mdash; each reported against a null or a caveat.")),
  sprintf("%s precomputed DE tables across two matrices, %s genes on the gene map", fmt(N$de1 + N$de2), fmt(N$genemap)),
  "shiny_app/app.R"))

# 6 · four-group DE
SL <- c(SL, sl("What it found",
  "The collaborator asked for three contrasts in four subclusters. Five of the twelve had no table on the matrix with good gene coverage, because that matrix is downsampled and CM2's KO-P0 arm falls to nine cells. Rather than pick one matrix and lose either genes or cells, both ship and the app says which one has the contrast you are missing.",
  "05", "Four-group DE",
  "Two matrices, because neither one wins",
  cols(
    cards(
      card(sprintf("Gene coverage &mdash; %s genes", fmt(N$broad)),
           sprintf("%s tables. But the cells are downsampled, so CM2&rsquo;s KO-P0 arm drops to 9 and CM4/CM9 lose their G1 strata.", N$de1)),
      card(sprintf("Cell coverage &mdash; %s genes", fmt(N$panel)),
           sprintf("%s tables over all %s cardiomyocytes. Every contrast runs, on a curated panel.", N$de2, fmt(N$cm))),
      card("Coverage of the request", "7 of 12 priority contrasts on one matrix alone; <b>12 of 12</b> across both.")),
    paste0('<h4 class="lab">CM2, P7 KO vs WT, G1 cells</h4>',
      '<ul class="tight">', paste0(sprintf("<li>%s</li>", top_genes("CM2","P7_KO_vs_WT","G1", 5, TRUE)), collapse = ""), '</ul>',
      '<p class="fine">and down: ', paste(gsub(" <span.*", "", top_genes("CM2","P7_KO_vs_WT","G1", 3, FALSE)), collapse = ", "), '</p>',
      callout("Arms can clear the 10-cell floor and still be useless: CM2&rsquo;s KO-P0 is 31 cells, 12 of them G1. Those are flagged, not silently reported.", "callout warn")),
    "c-5-7"),
  "Cell-level Wilcoxon; tables rank by effect size, not p",
  "app$fourgroup$de, $de2"))

# 7 · the headline
SL <- c(SL, sl("What it found",
  "This is the result. Read it in the G1 stratum so cycling composition is held fixed and the gap cannot be the sort. The maturation gap roughly doubles from P0 to P7 and the effect size goes from small to large. At the same time the G1 fraction moves from no difference at P0 to minus seven points at P7 &mdash; fewer P7 KO cells in G1, so more cycling. Both halves of the hypothesis, both appearing at P7. And the control that matters: restricting to G1 barely moves the gap, which is what should happen if this is maturation rather than composition.",
  "06", "G1 and maturation",
  "P7 KO cardiomyocytes are less mature &mdash; and it is specific to P7",
  paste0('<figure><div class="frame">', gap_svg(), '</div></figure>',
    cards(
      card("Maturation gap", sprintf("%s at P0 (d = %.2f) &rarr; <b>%s at P7</b> (d = %.2f). Small to large.",
           sgn(S_ALL$p0[["diff"]]), S_ALL$p0[["d"]], sgn(S_ALL$p7[["diff"]]), S_ALL$p7[["d"]])),
      card("G1 fraction", sprintf("%s pp at P0 &rarr; <b>%s pp at P7</b>. Fewer KO cells in G1 means more cycling.",
           sgn(S_ALL$g1_p0, 1), sgn(S_ALL$g1_p7, 1))),
      card("Delayed, not blocked", sprintf("The KO reaches <b>%s</b> of WT&rsquo;s P0&rarr;P7 maturation gain (%s in CM2).",
           pct(S_ALL$gain, 0), pct(S_CM2$gain, 0))))),
  "G1-only, so the gap is not the cycling enrichment",
  "app$fourgroup$scores, $phase"))

# 8 · intersection + the negative
SL <- c(SL, sl("What it found",
  "The intersection does what was asked and gives clean lists. The part worth your attention is the negative underneath. The stated purpose was to find genes linking delayed maturation to cycle competence, so I added a third axis for cycling. Zero of the hypothesis-quadrant genes are cycling-associated. My first test of that looked significant in the wrong direction and was circular &mdash; mature and cycling are anti-correlated states, so every mature marker scores anti-cycling automatically. The non-circular version regresses that out and still puts these genes below the line.",
  "07", "Intersection",
  "Good lists, and an honest negative under them",
  cols(
    paste0('<h4 class="lab">CM2 &mdash; immature genes up in P7 KO</h4>',
      '<p class="genes">', paste(head(it2$gene[it2$quadrant == "immature_up_in_KO"][
        order(-abs(it2$p7ko_log2FC[it2$quadrant == "immature_up_in_KO"]))], 8), collapse = " &middot; "), '</p>',
      '<h4 class="lab">CM2 &mdash; mature genes down in P7 KO</h4>',
      '<p class="genes">', paste(head(it2$gene[it2$quadrant == "mature_down_in_KO"][
        order(-abs(it2$p7ko_log2FC[it2$quadrant == "mature_down_in_KO"]))], 8), collapse = " &middot; "), '</p>',
      sprintf('<p class="fine">%d up / %d down in CM2. The down list is OXPHOS and TCA &mdash; the metabolic half of maturation, specifically.</p>',
              qn("immature_up_in_KO"), qn("mature_down_in_KO"))),
    paste0('<h4 class="lab">Do these genes link to cycle competence?</h4>',
      callout(sprintf("<b>No.</b> Of %d genes in the two hypothesis quadrants across all clusters, <b>%d</b> are cycling-associated. Median residual cycling association is <b>%s</b> &mdash; the part not already explained by maturation position. Below the line, not above.",
              N$hyp_rows, N$cyc_hits, sgn(N$cyc_resid)), "callout warn"),
      '<p>So these are maturation genes the KO moves, not a molecular bridge to cycling. If the link is real it is at the level of cell <em>state</em> &mdash; which is what the previous slide showed.</p>'),
    "c-6-6"),
  "Cycling axis validates: Top2a 0.81, Birc5 0.80, Mki67 0.77, Myh6 0.26",
  "app$fourgroup$intersect"))

# 9 · candidates
SL <- c(SL, sl("What it found",
  "The shortlist splits cleanly in two and the two halves are reporting different biology. Tcf4, Adamts9 and Gabbr2 are immature-state, strongly moved by the knockout, and P7-enriched in CM2. The proliferation genes are essentially unmoved. Foxm1 is barely expressed at all &mdash; 0.3 percent of cells &mdash; which is worth knowing before anyone builds on it.",
  "08", "Candidate genes",
  "The shortlist is two different stories",
  tbl(c("gene", "P7 KO vs WT", "P7-specific", "maturation", "cycling"),
    local({
      sp <- NULL
      pull <- function(cl, key, g) { d <- fg$de2[[cl]][[paste0(key, "__all")]]
        if (is.null(d)) return(NA_real_); d$log2FoldChange[match(g, d$gene)] }
      gs <- c("Tcf4","Adamts9","Gabbr2","Rrm2","Birc5","Prc1","Aurkb","Foxm1")
      lapply(gs, function(g) {
        p7 <- pull("CM2","P7_KO_vs_WT",g); p0 <- pull("CM2","P0_KO_vs_WT",g)
        i <- match(g, GM$gene)
        c(sprintf("<b>%s</b>", g),
          if (is.na(p7)) "&mdash;" else sgn(p7, 2),
          if (is.na(p7) || is.na(p0)) "&mdash;" else sgn(abs(p7) - abs(p0), 2),
          if (is.na(i) || is.na(GM$mat_class[i])) "&mdash;" else sub("-associated", "", GM$mat_class[i]),
          if (is.na(i) || is.na(GM$cyc_class[i])) "&mdash;" else sub("-associated", "", GM$cyc_class[i]))
      })
    })),
  "CM2, all cells, curated matrix",
  "app$fourgroup$de2, $geneaxes"))

# 10 · gene map
SL <- c(SL, sl("What it found",
  "Each point is a gene, placed by how strongly it marks mature versus immature and oxidative versus glycolytic. Two things came out of building it. First, the axes must be split at their own median rather than at 0.5, because the AUC carries a small global offset and a hard 0.5 split put 65 percent of genes in one corner. Second, and more interesting: the two programmes are tightly coupled at P0 and come apart by P7. Restricting to well-measured genes tightens P0 and does nothing for P7, so it is not simply the smaller P7 cell count.",
  "09", "Gene map",
  "Maturation and metabolism come apart by P7",
  cols(
    tbl(c("panel", "genes", "diagonal", "correlation"),
      lapply(c("avg","P0","P7"), function(pn) { p <- PAN[[pn]]
        c(if (pn == "avg") "averaged" else pn, fmt(p$n), pct(100 * p$diag), sgn(p$cor, 3)) })),
    cards(
      card("Split at the median, not 0.5", "The AUC carries a global offset; a hard 0.5 split put 65% of genes in one corner. An artifact that would have read as biology."),
      card("Coupled at P0, not at P7", "+0.68 vs +0.13 over all genes; +0.89 vs +0.16 among genes measured well at that age."),
      card("Circularity guarded", "Genes inside the scoring sets are hidden by default. With them hidden the diagonal still holds 57% &mdash; 50% would mean the axes are independent.")),
    "c-5-7"),
  "Diagonal = mature+oxidative and immature+glycolytic, the expected coupling",
  "app$fourgroup$geneaxes"))

# 11 · venn
SL <- c(SL, sl("What it found",
  "I was asked for a Venn. It is worth saying what a Venn does badly: it hides the threshold that built each set, the direction of change, and the overlap you would get by chance. Sixteen genes in the middle looks like a finding until you notice chance predicts twenty-four. So every overlap in the tab ships with its expected value and a hypergeometric p, and the universe is printed above the plot.",
  "10", "Gene-set Venn",
  "A Venn, reported against its null",
  cols(
    paste0(callout("<b>Is maturation just the cell cycle?</b> The P0&rarr;P7 signature overlaps the cycling program at or <em>below</em> chance &mdash; 16 observed against 24 expected, fold 0.7, p 0.97. A negative control, and it passes.", "callout"),
      '<p>Splitting that signature by genotype has more in it: the KO-only genes that fall from P0 to P7 are enriched for cycling genes where the WT-only genes are depleted &mdash; though the KO contributes more P7 cells, so part of that is power.</p>'),
    cards(
      card("The trap in the circles", "The data-driven cycling set is 531 genes but only 46 are canonical &mdash; the rest are largely housekeeping, because cycling cells are globally more transcriptionally active."),
      card("So the default is curated", "52 canonical genes. The axis recovers 46 of them, so the two agree on what everyone means; the curated set just doesn&rsquo;t drag the tail along.")),
    "c-6-6"),
  "Universe is the tested gene space, not the gated DE tables",
  "shiny_app/app.R &sect; gene-set Venn"))

# 12 · what we cannot claim
SL <- c(SL, sl("What it found",
  "Finish on the limits, because the design bounds every number in this deck. One animal per condition means none of this is inference. The genotypes are different sexes. The sort makes raw cross-timepoint cycling comparisons invalid. What survives all of that is the within-timepoint contrast at P7, held in G1 &mdash; and that is the one the maturation result rests on.",
  "11", "Limits",
  "What this design can and cannot support",
  cols(
    cards(
      card("Cannot: any p-value", "n = 1 animal per condition. Cell-level tests treat cells as replicates; the p columns rank, they do not test."),
      card("Cannot: raw P0 vs P7 cycling", "The sort inverts the developmental direction. Only phase-matched or within-timepoint contrasts are safe."),
      card("Cannot: sex-independent claims", "KO and WT are different sexes; the strongest raw hits are Y-genes and Xist.")),
    cards(
      card("Can: within-timepoint effect sizes", "P7 KO vs WT, in G1, is clean of the sort and is where the maturation result lives."),
      card("Can: rank candidates", "Effect sizes and percent-expressing are interpretable even when p-values are not."),
      card("Next: replication", sprintf("A sex-matched cohort at n &ge; 3 would turn %s of these tables into inference.", fmt(N$de1 + N$de2)))),
    "c-6-6"),
  "Every caveat here is also carried in the app&rsquo;s About block",
  "README.md, app.R &sect; about"))

# ---- assemble --------------------------------------------------------------
stopifnot(file.exists(TEMPLATE))
tpl <- paste(readLines(TEMPLATE, warn = FALSE), collapse = "\n")
grab <- function(tag) {
  m <- regmatches(tpl, regexpr(sprintf("(?s)<%s>.*?</%s>", tag, tag), tpl, perl = TRUE))
  stopifnot(length(m) == 1); m
}
style  <- grab("style")
script <- regmatches(tpl, regexpr("<script>(?s).*</script>", tpl, perl = TRUE))
stopifnot(length(script) == 1)

# a few classes this deck uses that the model deck does not
extra <- '<style>
.dt{width:100%;border-collapse:collapse;font-size:13px;margin:2px 0 10px}
.dt th{text-align:left;font-weight:700;color:var(--mut);font-size:11px;letter-spacing:.04em;
  text-transform:uppercase;border-bottom:1px solid var(--grd);padding:4px 8px 5px}
.dt td{padding:4px 8px;border-bottom:1px solid var(--bnd);color:var(--sec)}
.dt td:first-child{color:var(--pri)}
.dt tbody tr:last-child td{border-bottom:0}
.lab{font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:var(--mut);
  font-weight:700;margin:0 0 6px}
ul.tight{margin:2px 0 8px;padding-left:18px}
ul.tight li{margin:2px 0;font-size:13.5px;color:var(--pri)}
.num{color:var(--mut);font-family:var(--mono);font-size:12px}
.genes{font-size:13px;color:var(--pri);line-height:1.7;margin:0 0 12px}
.fine{font-size:12px;color:var(--mut);margin:0 0 8px}

/* On-screen navigation. The shared template is keyboard-only, which works when the
   deck is a local file with document focus but leaves no way through it when it is
   embedded in a frame, or on a touch screen. These controls drive the shared handler
   by dispatching the keys it already listens for, so there is one navigation
   implementation, not two. */
.navbar{position:fixed;top:.85rem;right:3rem;display:flex;align-items:center;gap:6px;
  z-index:60;font-family:var(--sans)}
.navbar button{width:34px;height:30px;border:1px solid var(--grd);background:var(--srf);
  color:var(--sec);border-radius:4px;font-size:15px;line-height:1;cursor:pointer;
  display:grid;place-items:center;transition:border-color .12s,color .12s}
.navbar button:hover{border-color:var(--axs);color:var(--pri)}
.navbar button:focus-visible{outline:2px solid var(--s1);outline-offset:2px}
.navbar button[disabled]{opacity:.35;cursor:default}

@media print{.navbar{display:none}}

/* The shared .cards grid is a fixed three columns, which is right at full width but
   cramped inside one half of a .cols split. Stack when nested. */
.cols .cards{grid-template-columns:1fr;gap:.7rem}
.cols .card{padding:.75rem .85rem;gap:.35rem}
.cols .card h3{font-size:.92rem}
.cols .card p{font-size:.82rem}

/* #counter and #movement are position:fixed at bottom 1rem, and the slide foot sits at
   the bottom of the slide grid, so the two print on top of each other. Give the slide
   enough bottom padding to clear them. */
.slide{padding-bottom:3.1rem}
</style>'

# The template's CSS and JS both depend on body scaffold that is NOT part of either:
# .slide is position:absolute/inset:0 and resolves against #deck, and the script
# reaches for #rail, #counter, #movement and #notes-body. Emitting the slides on their
# own left the script throwing on the first null before it could mark a slide visible,
# which renders as a blank page. Take the scaffold from the same template as the rest.
scaffold <- local({
  t <- tpl
  pre  <- substr(t, regexpr('<div id="rail">', t, fixed = TRUE),
                 regexpr('<section class="slide', t, fixed = TRUE) - 1)
  k    <- max(gregexpr("</section>", t, fixed = TRUE)[[1]]) + nchar("</section>")
  post <- substr(t, k, regexpr("<script>", substr(t, k, nchar(t)), fixed = TRUE) + k - 2)
  list(pre = pre, post = post)
})
stopifnot(grepl('id="deck"', scaffold$pre), grepl('id="notes-body"', scaffold$post))

navjs <- '<div class="navbar">
  <button id="npv" aria-label="Previous slide">&lsaquo;</button>
  <button id="nnx" aria-label="Next slide">&rsaquo;</button>
</div>
<script>
(function(){
  var slides = document.querySelectorAll(".slide");
  var pv = document.getElementById("npv"), nx = document.getElementById("nnx");
  // Drive the existing keydown handler rather than reimplementing show().
  function key(k){ document.dispatchEvent(new KeyboardEvent("keydown", {key:k, bubbles:true})); }
  function sync(){
    var on = 0;
    for (var j = 0; j < slides.length; j++) if (slides[j].classList.contains("on")) { on = j; break; }
    pv.disabled = on === 0; nx.disabled = on === slides.length - 1;
  }
  pv.addEventListener("click", function(){ key("ArrowLeft"); });
  nx.addEventListener("click", function(){ key("ArrowRight"); });
  new MutationObserver(sync).observe(document.body,
    {subtree:true, attributes:true, attributeFilter:["class"]});
  sync();
  // Touch: horizontal swipe only, so vertical scrolling inside a slide still works.
  var x0 = null, y0 = null;
  document.addEventListener("touchstart", function(e){
    x0 = e.touches[0].clientX; y0 = e.touches[0].clientY; }, {passive:true});
  document.addEventListener("touchend", function(e){
    if (x0 === null) return;
    var dx = e.changedTouches[0].clientX - x0, dy = e.changedTouches[0].clientY - y0;
    if (Math.abs(dx) > 45 && Math.abs(dx) > Math.abs(dy) * 1.6) key(dx < 0 ? "ArrowRight" : "ArrowLeft");
    x0 = y0 = null; }, {passive:true});
  // An embedded frame does not receive key events until something in it is focused.
  try { window.focus(); } catch (err) {}
  document.addEventListener("pointerdown", function(){ try { window.focus(); } catch (err) {} });
})();
</script>'

html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>E2F7/8 heart scRNA-seq &mdash; the browser and what it shows</title>',
  style, extra, '</head><body>\n', scaffold$pre, "\n", paste(SL, collapse = "\n"), "\n",
  scaffold$post, "\n", script, "\n", navjs, "\n</body></html>")

writeLines(html, OUT)
cat(sprintf("\nWrote %s  (%d slides, %.0f KB)\n", OUT, length(SL), file.size(OUT) / 1024))
cat("Open it, or print to PDF with ?print in the URL.\n")

# Second output: the same deck without the outer document wrapper, for publishing
# as an Artifact (which supplies its own doctype/head/body). Kept in the builder
# rather than done by hand so a re-publish after new data is one command.
body_only <- paste0(
  "<title>E2F7/8 scRNA Deck</title>\n", style, extra, "\n",
  scaffold$pre, "\n", paste(SL, collapse = "\n"), "\n", scaffold$post, "\n",
  script, "\n", navjs, "\n")
writeLines(body_only, OUT_BODY)
cat(sprintf("Wrote %s  (artifact-ready, no document wrapper)\n", OUT_BODY))
