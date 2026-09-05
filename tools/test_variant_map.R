#!/usr/bin/env Rscript
# test_variant_map.R
# ---------------------------------------------------------------------------
# The Variant explorer's UMAP must show the clustering it claims to show.
#
# What this guards, and why it is not obvious. The map does NOT read the variant's own
# label vector -- that covers only the 21,598 cells in the DE matrix. It colours
# app$pcdims$cm$percell (all 42,416 cells, at each of dims 10/30/50) by the column
# SCT_snn_res.<r>. That is only correct because those columns were verified to reproduce
# the registry's cluster counts for all nine variants. If the sweep is ever rebuilt with a
# different column convention, the map would silently colour by the WRONG clustering --
# it would still render, still look plausible, and be wrong. So this asserts the cluster
# counts per dims against app$clusterings$registry, read straight from the bundle.
#
#   docker run --rm --network none -v "$PWD:/repo" -u "$(id -u):$(id -g)" -e HOME=/tmp \
#     lab-server-e2f-heart-scrna-dev:latest Rscript /repo/tools/test_variant_map.R
# ---------------------------------------------------------------------------
library(shiny)
.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
setwd(normalizePath(file.path(dirname(.this), "..", "shiny_app")))
FAIL <- 0; ok <- function(l,c){ if(!isTRUE(c)) FAIL<<-FAIL+1
  cat(sprintf("  [%s] %s\n", if(isTRUE(c))"PASS" else "FAIL", l)) }
# The registry is static bundle data, so read it straight from the file rather than
# reaching into the app's environment -- fewer moving parts, and it checks the map
# against the SAME source of truth the app claims to follow.
REG <- readRDS("app_data.rds")$clusterings$registry
want_for <- function(res) { r <- REG[abs(REG$resolution - as.numeric(res)) < 1e-9, ]
                            as.integer(r$n_clusters[order(r$dims)]) }

testServer(shinyAppFile("app.R"), {
  cat("== the map shows the selected variant's clustering, at all three cuts ==\n")
  for (v in c("cm_dims10_res0.2","cm_dims30_res0.2","cm_dims50_res0.3")) {
    res <- REG$resolution[REG$variant_id == v]
    session$setInputs(clu_var = v, clu_map_col = "cluster", clu_map_scope = "all")
    d <- clu_map_g()$percell
    per <- as.integer(sapply(sort(unique(d$dims)), function(x) length(unique(d$cluster[d$dims==x]))))
    cat(sprintf("    %-18s res %-4s clusters per dims %-10s registry %s\n", v, res,
        paste(per, collapse="/"), paste(want_for(res), collapse="/")))
    ok(paste("matches registry:", v),   identical(per, want_for(res)))
    ok(paste("all three cuts:", v),     identical(sort(unique(d$dims)), c(10L,30L,50L)))
    ok(paste("full compartment:", v),   nrow(d) == 3*42416)
  }
  cat("\n== 'only the selected variant' narrows to one cut ==\n")
  session$setInputs(clu_var = "cm_dims50_res0.3", clu_map_scope = "one")
  d <- clu_map_g()$percell
  ok("single dims panel", identical(sort(unique(d$dims)), 50L))
  ok("42,416 cells",      nrow(d) == 42416)

  cat("\n== every colour-by option actually DRAWS ==\n")
  # This must force the plot to render to a device. testServer's output$clu_map does NOT
  # draw -- it returns a descriptor -- so an earlier version of this test passed while the
  # panel errored on first paint for every user. Build the ggplot and print it to a null
  # PNG device: that is the step that executes pcdims_gg and would have caught it.
  session$setInputs(clu_map_scope = "all")
  draws <- function() { f <- tempfile(fileext=".png")
    r <- try({ grDevices::png(f, 900, 400); on.exit({grDevices::dev.off(); unlink(f)}, add=TRUE)
               print(clu_map_p()); TRUE }, silent = TRUE)
    isTRUE(r) }
  for (cb in c("cluster","genotype","timepoint","Phase")) {
    session$setInputs(clu_map_col = cb)
    ok(paste("colour by", cb, "draws"), draws())
  }
  cat("\n== the empty-selection window that broke it ==\n")
  # selectInput's value is "" before its observe populates the choices. req() must hold
  # the plot back rather than passing "" into pcdims_gg.
  session$setInputs(clu_map_col = "")
  r <- try(clu_map_p(), silent = TRUE)
  ok("empty colour-by is held by req(), not passed through",
     inherits(r, "try-error") && grepl("silent|argument", conditionMessage(attr(r,"condition"))) ||
     inherits(attr(r,"condition"), "shiny.silent.error"))
  session$setInputs(clu_map_col = "cluster")
  cat("\n== the note ==\n")
  h <- output$clu_map_note$html
  ok("renders",                       nchar(h) > 300)
  ok("warns colours not comparable",  grepl("not comparable between panels", h))
  ok("states cluster counts per cut", grepl("clusters at dims 10", h))
})
cat(sprintf("\n%s\n", if(FAIL==0) "ALL PASS" else paste(FAIL,"FAILURES")))
if (FAIL) quit(status=1)
