# 99_verify.R -- checks the deliverable rather than trusting it.
# Run after 01-04. Any FAIL line means do not send the files.
suppressMessages({ library(openxlsx); library(Matrix) })
OUT <- "/out"
app <- readRDS("/in/app_data.rds")
de  <- readRDS(file.path(OUT, "de_tables.rds"))
en  <- readRDS(file.path(OUT, "enrich.rds"))
fails <- 0L
chk <- function(ok, what, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", what,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
  if (!ok) fails <<- fails + 1L
}

cat("\n== 1. regression against the bundle's own four-group grid ==\n")
ref <- app$fourgroup$de$CM1$P7_KO_vs_WT__all
new <- de$tables$P7_KO_vs_WT__CM1__all
k <- match(ref$gene, new$gene)
chk(!anyNA(k), "every gated reference gene is present in the ungated recompute")
chk(max(abs(ref$log2FoldChange - new$log2FoldChange[k])) == 0, "log2FC identical to the bundle")
chk(max(abs(ref$auc - new$auc[k])) == 0, "AUC identical to the bundle")
chk(identical(signif(ref$pvalue, 3), signif(new$pvalue[k], 3)), "p-values identical to the bundle")

cat("\n== 2. arm sizes match app$fourgroup$counts$n_in_de_matrix ==\n")
cc <- app$fourgroup$counts
bad <- character()
for (i in seq_len(nrow(de$manifest))) {
  m <- de$manifest[i, ]
  if (m$stratum != "all") next
  for (side in c("A", "B")) {
    grp <- m[[paste0("arm_", side)]]
    exp <- cc$n_in_de_matrix[cc$cluster == m$cluster & cc$group == grp]
    got <- m[[paste0("n_", side)]]
    if (length(exp) == 1 && exp != got) bad <- c(bad, sprintf("%s/%s %d!=%d", m$key, grp, got, exp))
  }
}
chk(!length(bad), "all all-cells arm sizes match the bundle", paste(bad, collapse = "; "))

cat("\n== 3. completeness ==\n")
ng <- nrow(app$deg_expr)
chk(all(vapply(de$tables, nrow, 0L) == ng), sprintf("every table has all %d genes", ng))
req <- c("E2f7","E2f8","Mki67","Top2a","Myh6","Myh7","Nppa","Ect2","Nppb","Pln")
chk(all(vapply(de$tables, function(d) all(req %in% d$gene), TRUE)), "all watch-list genes present in every table")

cat("\n== 4. confounders ==\n")
chk(all(vapply(de$tables, function(d) all(d$confounder == (d$gene %in% app$confound)), TRUE)),
    "confounder flag matches app$confound exactly")
gochk <- TRUE
if (!is.null(en$go)) {
  gochk <- !any(vapply(strsplit(en$go$geneID, "/"), function(g) any(g %in% app$confound), TRUE))
}
chk(gochk, "no confounder gene appears in any enriched GO term's gene list")

cat("\n== 5. GO sanity: term gene lists are subsets of their input list ==\n")
if (is.null(en$go)) chk(FALSE, "GO results exist") else {
  bad <- 0L; n <- 0L
  for (i in sample(seq_len(nrow(en$go)), min(60, nrow(en$go)))) {
    r <- en$go[i, ]
    key <- sprintf("%s__%s__%s", r$contrast, r$cluster, r$stratum)
    d <- de$tables[[key]]
    sgn <- if (grepl("_up$", r$direction) && r$direction != "P0_up") 1 else -1
    if (r$direction == "P0_up") sgn <- -1
    pool <- d$gene[!d$confounder & sgn * d$log2FoldChange > 0]
    n <- n + 1L
    if (!all(strsplit(r$geneID, "/")[[1]] %in% pool)) bad <- bad + 1L
  }
  chk(bad == 0, sprintf("%d sampled GO terms: every gene is in that direction's candidate pool", n),
      if (bad) sprintf("%d bad", bad) else "")
  chk(all(en$go$n_universe > en$go$n_input), "universe is always larger than the input list")
  chk(all(en$go$n_universe < ng), "universe is the expressed-gene space, not all genes")
}

cat("\n== 6. workbooks ==\n")
for (f in c("P7_KO_vs_WT_by_CM_subcluster.xlsx", "P7WT_vs_P0WT.xlsx")) {
  p <- file.path(OUT, f)
  chk(file.exists(p), paste("exists:", f))
  if (!file.exists(p)) next
  sn <- getSheetNames(p)
  chk(all(nchar(sn) <= 31), paste0(f, ": all sheet names <= 31 chars"))
  chk(!any(duplicated(sn)), paste0(f, ": no duplicate sheet names"))
  chk("README" %in% sn && "Summary" %in% sn && "GO_audit" %in% sn,
      paste0(f, ": README / Summary / GO_audit present"))
  chk(sum(grepl("^DE_", sn)) >= 7, paste0(f, ": at least 7 DE sheets"), paste(sum(grepl("^DE_", sn)), "found"))
  sz <- file.size(p) / 1e6
  chk(sz < 25, sprintf("%s: %.1f MB (emailable)", f, sz))
  first <- read.xlsx(p, sheet = grep("^DE_", sn, value = TRUE)[1])
  chk(nrow(first) > 10000, sprintf("%s: first DE sheet has %d rows", f, nrow(first)))
  chk("direction" %in% names(first) && "sig_sets" %in% names(first),
      paste0(f, ": DE sheets carry direction and sig_sets"))
}

cat("\n== 7. figures ==\n")
for (part in c("part1", "part2")) {
  fs <- list.files(file.path(OUT, "plots", part), pattern = "\\.png$")
  chk(length(fs) > 0, sprintf("%s has PNGs", part), sprintf("%d files", length(fs)))
  pdfs <- list.files(file.path(OUT, "plots", part), pattern = "\\.pdf$")
  chk(length(fs) == length(pdfs), sprintf("%s: PNG and PDF counts match", part))
  chk(all(file.size(file.path(OUT, "plots", part, fs)) > 5000), sprintf("%s: no truncated PNGs", part))
}
prim <- de$manifest[de$manifest$is_primary, ]
missing <- character()
for (i in seq_len(nrow(prim))) {
  part <- if (prim$contrast[i] == "P7_KO_vs_WT") "part1" else "part2"
  for (dirn in c(prim$up_label[i], prim$down_label[i])) {
    f <- file.path(OUT, "plots", part, sprintf("go_bp_%s_%s_%s.png", prim$contrast[i], prim$cluster[i], dirn))
    if (!file.exists(f)) missing <- c(missing, basename(f))
  }
}
chk(!length(missing), "every primary cluster x direction has a GO panel (real or placeholder)",
    paste(head(missing, 5), collapse = ", "))

cat("\n== 8. part 2 biology smoke test (all CM, G1-matched) ==\n")
g1 <- de$tables$WT_P7_vs_P0__AllCM__G1
lf <- function(g) g1$log2FoldChange[match(g, g1$gene)]
chk(all(lf(c("Myh7","Nppa","Nppb","Myl7","Tnni1")) < 0), "immature markers are down at P7")
chk(all(lf(c("Acadm","Acadvl","Hadha","Ppargc1a")) > 0), "fatty-acid oxidation genes are up at P7")
chk(lf("Pln") > 0 && lf("Atp2a2") > 0, "Pln and Atp2a2 (mature) are up at P7")
chk(all(lf(c("Pgk1","Aldoa")) < 0), "glycolytic genes are down at P7")

cat(sprintf("\n===== %d failure(s) =====\n", fails))
if (fails) quit(status = 1)
