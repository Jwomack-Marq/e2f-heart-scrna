#!/usr/bin/env Rscript
# One-time enrichment: add a per-gene info table (app$geneinfo) to app_data.rds so
# the browser can show name / aliases / location / summary + external links for the
# gene picked on a volcano plot. Build-time only (needs internet); the deployed app
# makes NO runtime network calls.
#
#   Rscript enrich_geneinfo.R [path/to/app_data.rds]
#
# Mouse fields come from MyGene.info (mouse); the functional summary is borrowed
# from the human ortholog (mouse genes rarely carry an NCBI summary).

suppressPackageStartupMessages({ library(httr); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
RDS  <- if (length(args)) args[1] else
  "C:/Users/Justi/OneDrive/Documents/GitHub/e2f-heart-scrna/shiny_app/app_data.rds"
stopifnot(file.exists(RDS))

app <- readRDS(RDS)

# ---- gene universe = every gene the volcanoes can show -----------------------
g <- character(0)
for (d in app$tables$ct_DE) if (!is.null(d$gene)) g <- c(g, d$gene)
for (res in app$tables$sub_DE) for (d in res) if (!is.null(d$gene)) g <- c(g, d$gene)
genes <- sort(unique(g))
cat(sprintf("gene universe: %d symbols\n", length(genes)))

# ---- helper: MyGene batch POST -> flattened data.frame -----------------------
mygene_post <- function(terms, species, fields, scopes = "symbol") {
  out <- list()
  chunks <- split(terms, ceiling(seq_along(terms) / 1000))
  for (i in seq_along(chunks)) {
    resp <- tryCatch(POST("https://mygene.info/v3/query",
                          body = list(q = paste(chunks[[i]], collapse = ","),
                                      scopes = scopes, fields = fields, species = species),
                          encode = "form", timeout(60)),
                     error = function(e) NULL)
    if (is.null(resp) || http_error(resp)) { cat("  chunk", i, "FAILED\n"); next }
    df <- tryCatch(fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE),
                   error = function(e) NULL)
    if (is.data.frame(df)) out[[length(out) + 1]] <- df
    cat(sprintf("  %s chunk %d/%d (%d rows)\n", species, i, length(chunks),
                if (is.data.frame(df)) nrow(df) else 0L))
    Sys.sleep(0.3)
  }
  if (!length(out)) return(NULL)
  # union columns across chunks
  cols <- unique(unlist(lapply(out, names)))
  do.call(rbind, lapply(out, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}

# scalar-ify possibly list/multi-valued columns; take first non-NA
first1 <- function(x) { if (is.list(x)) x <- x[[1]]; x <- x[!is.na(x)]; if (length(x)) x[1] else NA }
getcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else rep(NA, nrow(df))

# ---- 1) mouse annotation -----------------------------------------------------
cat("Querying MyGene (mouse) ...\n")
m <- mygene_post(genes, "mouse",
       "name,alias,type_of_gene,genomic_pos.chr,genomic_pos.start,genomic_pos.end,entrezgene,ensembl.gene,MGI")
stopifnot(!is.null(m), "query" %in% names(m))
m <- m[!duplicated(m$query), ]                      # first (best-scoring) hit per symbol
rownames(m) <- m$query

# ---- 2) human-ortholog summary (uppercased symbol) ---------------------------
cat("Querying MyGene (human summary) ...\n")
hsym <- toupper(genes)
h <- mygene_post(unique(hsym), "human", "summary")
hsummary <- setNames(rep("", length(genes)), genes)
if (!is.null(h) && "query" %in% names(h) && "summary" %in% names(h)) {
  h <- h[!duplicated(h$query), ]
  hs <- setNames(vapply(h$summary, function(s) { s <- first1(s); if (is.na(s)) "" else s }, ""), h$query)
  hit <- hsym %in% names(hs)
  hsummary[hit] <- hs[hsym[hit]]
}

# ---- 3) assemble app$geneinfo data.frame keyed by mouse symbol ---------------
pull <- function(sym, col) if (sym %in% rownames(m)) first1(m[sym, col]) else NA
mk <- function(sym, col) { v <- pull(sym, col); if (is.na(v)) "" else as.character(v) }
gi <- data.frame(
  name    = vapply(genes, mk, "", col = "name"),
  alias   = vapply(genes, function(s) { a <- if (s %in% rownames(m)) m[s, "alias"] else NA
                     if (is.list(a)) a <- a[[1]]; a <- a[!is.na(a)]
                     if (length(a)) paste(a, collapse = ", ") else "" }, ""),
  type    = vapply(genes, mk, "", col = "type_of_gene"),
  chr     = vapply(genes, mk, "", col = "genomic_pos.chr"),
  start   = vapply(genes, mk, "", col = "genomic_pos.start"),
  end     = vapply(genes, mk, "", col = "genomic_pos.end"),
  entrez  = vapply(genes, mk, "", col = "entrezgene"),
  ensembl = vapply(genes, mk, "", col = "ensembl.gene"),
  mgi     = vapply(genes, mk, "", col = "MGI"),
  summary = unname(hsummary[genes]),
  stringsAsFactors = FALSE)
rownames(gi) <- genes

cat(sprintf("resolved names: %d/%d | summaries: %d/%d\n",
            sum(nzchar(gi$name)), nrow(gi), sum(nzchar(gi$summary)), nrow(gi)))
cat("spot-check Myh7:\n"); print(gi["Myh7", c("name","type","mgi","summary")])

# ---- 4) merge + save (back up first) -----------------------------------------
bak <- sub("\\.rds$", paste0(".bak.rds"), RDS)
if (!file.exists(bak)) file.copy(RDS, bak)
app$geneinfo <- gi
saveRDS(app, RDS, compress = "xz")
cat("saved:", RDS, "\n=== DONE enrich_geneinfo ===\n")
