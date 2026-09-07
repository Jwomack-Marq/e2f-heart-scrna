#!/usr/bin/env Rscript
# Determine the sex of each sample directly from canonical sex markers:
#   Xist  -> high in FEMALES (X-inactivation), ~silent in males
#   Ddx3y/Eif2s3y/Uty/Kdm5d -> Y-linked, expressed ONLY in males
# Reported per condition AND per lane (to show it's not a lane artifact).
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))

SEX <- list(female_Xist = "Xist",
            male_Y = c("Ddx3y","Eif2s3y","Uty","Kdm5d"))
markers <- unlist(SEX)

rows <- list()
for (tp in TIMEPOINTS) {
  obj <- readRDS(merged_path(tp)); DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
  obj$grp  <- obj$orig.ident
  obj$grpl <- paste(obj$orig.ident, obj$lane, sep = "_")
  present <- intersect(markers, rownames(obj))
  dat <- FetchData(obj, vars = present, layer = "data")
  cnt <- FetchData(obj, vars = present, layer = "counts")
  for (g in unique(obj$grpl)) {
    m <- obj$grpl == g
    for (gene in present) {
      rows[[paste(g, gene)]] <- data.frame(
        sample = sub("_lane[0-9]+$", "", g), lane = sub("^.*_(lane[0-9]+)$", "\\1", g),
        gene = gene, mean_logexpr = round(mean(dat[m, gene]), 3),
        pct_expressing = round(100 * mean(cnt[m, gene] > 0), 1))
    }
  }
  rm(obj); gc(verbose = FALSE)
}
res <- do.call(rbind, rows)
# call sex per sample from Xist vs Y-gene signal
library(dplyr)
calls <- res %>% group_by(sample) %>%
  summarise(Xist_pct = pct_expressing[gene == "Xist"][1],
            Ymax_pct = max(pct_expressing[gene %in% SEX$male_Y]),
            .groups = "drop") %>%
  mutate(sex_call = ifelse(Ymax_pct > Xist_pct, "MALE (Y+, Xist-low)", "FEMALE (Xist+, Y-low)"))

write.csv(res,   file.path(OUTTAB, "sex_markers_by_sample_lane.csv"), row.names = FALSE)
write.csv(calls, file.path(OUTTAB, "sex_calls.csv"), row.names = FALSE)
cat("\n--- sex markers (% of cells expressing), per sample x lane ---\n")
print(res, row.names = FALSE)
cat("\n--- sex call per sample ---\n"); print(as.data.frame(calls), row.names = FALSE)
cat("\n=== DONE sex_check ===\n")
