suppressPackageStartupMessages(library(Matrix))
OUT <- "/out"
app <- readRDS("/in/app_data.rds")
cat("=== top-level slots ===\n"); print(names(app))
cat("built:", if (!is.null(app$built)) app$built else "NA", "\n\n")

w <- function(df, name) {
  utils::write.csv(df, file.path(OUT, name), row.names = FALSE, na = "")
  cat(sprintf("  wrote %-34s %6d rows x %3d cols\n", name, nrow(df), ncol(df)))
}

# ---- per-cell metadata -----------------------------------------------------
if (!is.null(app$meta)) {
  m <- app$meta
  if (!"cell" %in% names(m)) m$cell <- rownames(m)
  cat("=== app$meta columns ===\n"); print(names(m))
  w(m, "meta.csv")
  for (col in c("celltype","genotype","timepoint","Phase","cycling","orig.ident")) {
    if (col %in% names(m)) { cat("--", col, ":\n"); print(table(m[[col]], useNA="ifany")) }
  }
}
if (!is.null(app$cm) && !is.null(app$cm$meta)) {
  cm <- app$cm$meta
  if (!"cell" %in% names(cm)) cm$cell <- rownames(cm)
  cat("\n=== app$cm$meta columns ===\n"); print(names(cm))
  w(cm, "cm_meta.csv")
}

# ---- gene coverage: what can the model actually be scored on? ---------------
MODEL_NODES <- c("E2f1","E2f2","E2f3","E2f4","E2f5","E2f6","E2f7","E2f8",
  "Rb1","Rbl1","Rbl2","Ccnd1","Ccnd2","Ccne1","Ccne2","Ccna2","Ccnb1",
  "Cdk1","Cdk2","Cdk4","Cdk6","Cdkn1a","Cdkn1b","Cdt1","Gmnn","Ccng1",
  "Chek1","Wee1","Pkmyt1","Cdc25a","Cdc25b","Cdc20","Fzr1","Skp2","Fbxo5",
  "Ect2","Rhoa","Anln","Aurkb","Racgap1","Kif23","Cep55","Mki67","Top2a",
  "Hmgcs2","Fabp3","Pdk4","Ucp2","Erbb2","Amotl1","Nisch","Adra1b","Adrb1",
  "Mapk12","Mapk14","Slc2a1","Ldha","Pkm","Hk1","Hk2","Timeless","Ung","Msh6")
if (!is.null(app$genes)) {
  g <- as.character(app$genes)
  cov <- data.frame(gene = MODEL_NODES, in_curated_panel = MODEL_NODES %in% g)
  w(cov, "model_gene_coverage.csv")
  cat(sprintf("\ncurated panel: %d genes; model nodes covered: %d / %d\n",
              length(g), sum(cov$in_curated_panel), nrow(cov)))
  cat("MISSING from panel:", paste(cov$gene[!cov$in_curated_panel], collapse=", "), "\n")
  writeLines(g, file.path(OUT, "curated_genes.txt"))
}

# ---- the phase-calling gene lists, and whether we can rescore ---------------
S_GENES <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
 "CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76",
 "SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1",
 "TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
G2M_GENES <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
 "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2",
 "AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20",
 "TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR",
 "AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
if (!is.null(app$genes)) {
  up <- toupper(as.character(app$genes))
  covS <- sum(S_GENES %in% up); covG <- sum(G2M_GENES %in% up)
  cat(sprintf("\nphase-list coverage in curated panel: S %d/%d, G2M %d/%d\n",
              covS, length(S_GENES), covG, length(G2M_GENES)))
  pl <- data.frame(
    list_name = c(rep("S", length(S_GENES)), rep("G2M", length(G2M_GENES))),
    gene      = c(S_GENES, G2M_GENES),
    in_panel  = c(S_GENES %in% up, G2M_GENES %in% up),
    is_ko_target = c(S_GENES, G2M_GENES) %in% c("E2F7","E2F8"))
  w(pl, "phase_list_coverage.csv")
}

# ---- expression for the model nodes, from the curated panel -----------------
if (!is.null(app$expr) && !is.null(app$genes)) {
  ex <- app$expr; rownames(ex) <- as.character(app$genes)
  keep <- intersect(MODEL_NODES, rownames(ex))
  sub <- as.matrix(ex[keep, , drop = FALSE])
  df <- data.frame(gene = rownames(sub), sub, check.names = FALSE)
  utils::write.csv(df, file.path(OUT, "model_node_expr.csv"), row.names = FALSE)
  cat(sprintf("  wrote %-34s %6d genes x %d cells\n", "model_node_expr.csv",
              nrow(sub), ncol(sub)))
  # also the two phase lists' panel members, for the rescore
  ph <- intersect(union(S_GENES, G2M_GENES), toupper(rownames(ex)))
  idx <- which(toupper(rownames(ex)) %in% ph)
  sub2 <- as.matrix(ex[idx, , drop = FALSE])
  df2 <- data.frame(gene = rownames(sub2), sub2, check.names = FALSE)
  utils::write.csv(df2, file.path(OUT, "phase_gene_expr.csv"), row.names = FALSE)
  cat(sprintf("  wrote %-34s %6d genes x %d cells\n", "phase_gene_expr.csv",
              nrow(sub2), ncol(sub2)))
}

# ---- DE tables -------------------------------------------------------------
if (!is.null(app$tables)) {
  cat("\n=== app$tables ===\n"); print(names(app$tables))
  ct <- app$tables$ct_DE
  if (!is.null(ct)) {
    cat("ct_DE keys:", paste(names(ct), collapse=", "), "\n")
    for (k in names(ct)) if (grepl("Cardiomyocyte", k)) w(ct[[k]], paste0("de_", k, ".csv"))
  }
  for (nm in c("P0_cardiac_DE","P7_cardiac_DE","e2f_regulon","cellcycle_fraction",
               "composition","sub_summary","sub_subtype")) {
    x <- app$tables[[nm]]
    if (is.data.frame(x)) w(x, paste0("tbl_", nm, ".csv"))
  }
  sd_ <- app$tables$sub_DE
  if (!is.null(sd_) && !is.null(sd_[["res0.2"]])) {
    for (k in names(sd_[["res0.2"]])) w(sd_[["res0.2"]][[k]], paste0("subDE_res0.2_", k, ".csv"))
  }
}
cat("\nDONE\n")
