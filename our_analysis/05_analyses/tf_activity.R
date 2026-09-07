#!/usr/bin/env Rscript
# TF / E2F regulon activity (DESCRIPTIVE, n=1, sex-confounded).
# Robust substitute for full SCENIC: decoupleR ULM over CollecTRI (mouse) regulons,
# on per-(celltype x sample x lane) pseudobulk profiles -> TF activity per group.
# Tests whether ACTIVATING-E2F regulon activity is higher in KO (de-repression, the
# E2F7/8-repressor hypothesis). Full SCENIC (GENIE3 + cisTarget DBs) is noted as a
# heavier pySCENIC/cluster follow-up.
this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))
note_scenic <- "Full SCENIC (GENIE3 co-expression + RcisTarget motif pruning + AUCell) needs ~1GB+ cisTarget feather DBs and long runtime; run via pySCENIC on the cluster if a motif-pruned regulon is required. Here: decoupleR ULM over CollecTRI regulons (literature/curated), which gives per-group TF activity including the E2F family."
if (!requireNamespace("decoupleR", quietly = TRUE)) {
  writeLines(c("decoupleR not installed; TF-activity skipped.", note_scenic), file.path(OUTTAB, "tf_activity_blocked.txt"))
  cat("decoupleR missing.\n=== DONE tf_activity ===\n"); quit(save = "no")
}
suppressWarnings(suppressMessages({ library(decoupleR); library(ggplot2); library(dplyr); library(tidyr); library(Matrix) }))

net <- tryCatch(decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE), error = function(e) NULL)
net_src <- "CollecTRI"
if (is.null(net) || nrow(net) < 1000) {
  # CollecTRI/Omnipath unavailable -> robust offline fallback: MSigDB GTRD TF->target
  # regulons + Hallmark E2F targets (all mor=1). Still yields E2F regulon activity.
  net_src <- "MSigDB_TFT_GTRD+Hallmark_E2F"
  suppressWarnings(suppressMessages(library(msigdbr)))
  m_tft <- msigdbr(species = "Mus musculus", collection = "C3", subcollection = "TFT:GTRD")
  m_h   <- msigdbr(species = "Mus musculus", collection = "H")
  m_h   <- m_h[m_h$gs_name == "HALLMARK_E2F_TARGETS", ]
  net <- unique(rbind(
    data.frame(source = m_tft$gs_name, target = m_tft$gene_symbol, mor = 1),
    data.frame(source = m_h$gs_name,   target = m_h$gene_symbol,   mor = 1)))
}
cat(sprintf("Regulon source: %s -- %d edges, %d regulons\n", net_src, nrow(net), length(unique(net$source))))

comb <- readRDS(file.path(PROC, "seurat.combined.annotated.rds")); DefaultAssay(comb) <- "RNA"
if (inherits(comb[["RNA"]], "Assay5")) comb[["RNA"]] <- SeuratObject::JoinLayers(comb[["RNA"]])
# Build pseudobulk by (celltype x sample x lane) with an explicit sparse model matrix
# so the group labels survive exactly (AggregateExpression sanitizes underscores).
counts <- GetAssayData(comb, assay = "RNA", layer = "counts")
grp <- factor(paste(comb$celltype, comb$orig.ident, comb$lane, sep = "@@"))
mm <- Matrix::sparse.model.matrix(~ 0 + grp); colnames(mm) <- sub("^grp", "", colnames(mm))
keep <- Matrix::colSums(mm) >= 10                    # drop tiny groups (noisy)
mm <- mm[, keep, drop = FALSE]
pb <- as.matrix(counts %*% mm)                       # genes x groups
cn <- colnames(pb)
sp <- do.call(rbind, strsplit(cn, "@@", fixed = TRUE))
meta <- data.frame(celltype = sp[,1], orig.ident = sp[,2], lane = sp[,3], stringsAsFactors = FALSE)
meta$genotype  <- genotype_of(meta$orig.ident)
meta$timepoint <- substr(meta$orig.ident, 1, 2)

# logCPM for ULM
cpm <- sweep(as.matrix(pb), 2, pmax(colSums(pb), 1), "/") * 1e6
logcpm <- log1p(cpm)

acts <- decoupleR::run_ulm(mat = logcpm, net = net, .source = "source", .target = "target",
                           .mor = "mor", minsize = 5)
acts <- acts[acts$statistic == "ulm", ]
acts <- merge(acts, cbind(condition = cn, meta), by = "condition")

## E2F family focus -----------------------------------------------------------
# E2F regulons: CollecTRI uses TF symbols (E2f1..); MSigDB uses E2F#_TARGET_GENES / HALLMARK_E2F_TARGETS
e2f <- c("E2f1","E2f2","E2f3","E2f4","E2f5","E2f6","E2f7","E2f8","Tfdp1","Tfdp2")
e2f_act <- acts[acts$source %in% e2f | grepl("E2F", acts$source, ignore.case = TRUE), ]
e2f_summ <- e2f_act |> group_by(source, celltype, timepoint, genotype) |>
  summarise(mean_activity = round(mean(score), 3), .groups = "drop")
e2f_wide <- e2f_summ |> pivot_wider(names_from = genotype, values_from = mean_activity) |>
  mutate(KO_minus_WT = round(KO - WT, 3))
write.csv(e2f_wide, file.path(OUTTAB, "e2f_regulon_activity.csv"), row.names = FALSE)

## all-TF summary by celltype x genotype (mean over timepoints/lanes) ---------
tf_summ <- acts |> group_by(source, celltype, genotype) |>
  summarise(mean_activity = round(mean(score), 3), .groups = "drop")
write.csv(tf_summ, file.path(OUTTAB, "tf_activity_by_celltype_genotype.csv"), row.names = FALSE)
writeLines(note_scenic, file.path(OUTTAB, "tf_activity_SCENIC_note.txt"))

## figures --------------------------------------------------------------------
png(file.path(OUTFIG, "tf_activity_E2F_heatmap.png"), 1000, 600)
print(ggplot(e2f_summ, aes(paste(celltype, timepoint, sep="@"), source, fill = mean_activity)) +
        geom_tile() + facet_wrap(~ genotype) +
        scale_fill_gradient2(low="steelblue", mid="white", high="firebrick") +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "E2F-family regulon activity (decoupleR/CollecTRI; descriptive)", x=NULL, y=NULL))
dev.off()
# top TFs by |KO-WT| difference (cardiomyocyte)
cm_diff <- tf_summ |> filter(celltype == "Cardiomyocyte") |>
  pivot_wider(names_from = genotype, values_from = mean_activity) |>
  mutate(diff = KO - WT) |> arrange(desc(abs(diff))) |> slice_head(n = 20)
png(file.path(OUTFIG, "tf_activity_top_KOvsWT.png"), 800, 600)
print(ggplot(cm_diff, aes(reorder(source, diff), diff, fill = diff>0)) + geom_col() + coord_flip() +
        scale_fill_manual(values=c(`TRUE`="firebrick",`FALSE`="steelblue"), guide="none") +
        labs(title="Cardiomyocyte TF activity KO-WT (top |diff|, descriptive)", x=NULL, y="KO - WT activity") + theme_bw())
dev.off()

cat("\n--- E2F-family regulon activity (KO - WT) ---\n")
print(as.data.frame(e2f_wide), row.names = FALSE)
cat("\n=== DONE tf_activity ===\n")
