#!/usr/bin/env Rscript
# Build a shareable HTML report + PowerPoint deck from the our_analysis results.
# HTML is self-contained (figures base64-embedded). PPTX via officer/flextable.
# All KO-vs-WT content is framed DESCRIPTIVE (n=1, sex-confounded).
#
# The report CONTENT (helpers + the ordered `sections` list) now lives in
# report_content.R so it is shared, unchanged, with interactive/report.Rmd.

this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(this), "_common.R"))                 # defines OUT/OUTFIG/OUTTAB
source(file.path(dirname(this), "report_content.R"))          # CONFOUND, fig/tab, png_dims, top_de, rnd, sections
suppressWarnings(suppressMessages({ library(base64enc) }))

# ============================ HTML =========================================
build_html <- function() {
  esc <- function(x) gsub("<","&lt;", gsub("&","&amp;", x))
  img_tag <- function(p) if (is.na(p) || !file.exists(p)) "" else
    sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;height:auto;border:1px solid #ddd;margin:6px 0;">', base64encode(p))
  df_html <- function(df) {
    if (is.null(df) || !nrow(df)) return("")
    h <- paste0("<th>", esc(names(df)), "</th>", collapse = "")
    rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", esc(as.character(r)), "</td>", collapse = ""), "</tr>"))
    paste0("<table><thead><tr>", h, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table>")
  }
  css <- "body{font-family:Segoe UI,Arial,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222;line-height:1.5}
h1{color:#7b1fa2}h2{color:#7b1fa2;border-bottom:2px solid #eee;padding-top:18px}
table{border-collapse:collapse;margin:8px 0;font-size:13px}th,td{border:1px solid #ccc;padding:3px 7px;text-align:right}th{background:#f3e5f5}
.cap{color:#666;font-size:13px;font-weight:600;margin:2px 0 2px}.warn{background:#fff3e0;border-left:4px solid #e65100;padding:8px 12px;margin:8px 0}
.desc{color:#333;font-size:13.5px;margin:0 0 16px;line-height:1.45;background:#faf6fc;border-left:3px solid #ce93d8;padding:6px 12px}"
  body <- ""
  for (i in seq_along(sections)) {
    s <- sections[[i]]
    body <- paste0(body, if (i == 1) sprintf("<h1>%s</h1>", esc(s$title)) else sprintf("<h2>%s</h2>", esc(s$title)))
    cls <- if (grepl("caveat", s$title, ignore.case = TRUE)) ' class="warn"' else ""
    if (length(s$narr)) body <- paste0(body, sprintf("<ul%s>", cls), paste0("<li>", esc(s$narr), "</li>", collapse = ""), "</ul>")
    for (f in s$figs) if (!is.null(f)) {
      body <- paste0(body, img_tag(f[1]), sprintf('<div class="cap">%s</div>', esc(f[2])))
      if (length(f) >= 3 && !is.na(f[3]) && nzchar(f[3]))
        body <- paste0(body, sprintf('<div class="desc">%s</div>', esc(f[3])))
    }
    for (t in s$tabs) if (!is.null(t$df)) body <- paste0(body, df_html(rnd(t$df)), sprintf('<div class="cap">%s</div>', esc(t$cap)))
  }
  html <- sprintf("<!DOCTYPE html><html><head><meta charset='utf-8'><title>E2F7/8 scRNA pilot</title><style>%s</style></head><body>%s<hr><p class='cap'>Generated from our_analysis/. Descriptive pilot (n=1).</p></body></html>", css, body)
  out <- file.path(OUT, "E2F_scRNA_pilot_report.html")
  writeLines(html, out); cat("Wrote", out, "\n")
}

# ============================ PPTX =========================================
build_pptx <- function() {
  if (!requireNamespace("officer", quietly = TRUE)) { message("officer missing; skipping pptx"); return(invisible()) }
  library(officer); has_ft <- requireNamespace("flextable", quietly = TRUE)
  ppt <- read_pptx()
  SW <- 10; SH <- 7.5  # default 4:3 template (inches)
  place <- function(ppt, path, top = 1.5, maxw = 9, maxh = 5.4) {
    d <- png_dims(path); w <- maxw; h <- w * d["h"]/d["w"]
    if (h > maxh) { h <- maxh; w <- h * d["w"]/d["h"] }
    ph_with(ppt, external_img(path, width = w, height = h),
            location = ph_location(left = (SW-w)/2, top = top, width = w, height = h))
  }
  for (i in seq_along(sections)) {
    s <- sections[[i]]
    if (i == 1) {
      ppt <- add_slide(ppt, "Title Slide", "Office Theme")
      ppt <- ph_with(ppt, s$title, ph_location_type(type = "ctrTitle"))
      ppt <- ph_with(ppt, paste(s$narr, collapse = "\n"), ph_location_type(type = "subTitle"))
      next
    }
    # text/table slide
    ppt <- add_slide(ppt, "Title and Content", "Office Theme")
    ppt <- ph_with(ppt, s$title, ph_location_type(type = "title"))
    if (length(s$narr))
      ppt <- ph_with(ppt, unlist(s$narr), ph_location(left = 0.5, top = 1.4, width = 9, height = 5.5))
    if (length(s$tabs) && has_ft) for (t in s$tabs) if (!is.null(t$df)) {
      ppt <- add_slide(ppt, "Title and Content", "Office Theme")
      ppt <- ph_with(ppt, paste0(s$title, " — ", t$cap), ph_location_type(type = "title"))
      ft <- flextable::fontsize(flextable::autofit(flextable::flextable(rnd(t$df))), size = 9, part = "all")
      ppt <- ph_with(ppt, ft, ph_location(left = 0.5, top = 1.4))
    }
    # one slide per figure
    for (f in s$figs) if (!is.null(f) && !is.na(f[1]) && file.exists(f[1])) {
      ppt <- add_slide(ppt, "Title and Content", "Office Theme")
      ppt <- ph_with(ppt, paste0(s$title, " — ", f[2]), ph_location_type(type = "title"))
      ppt <- place(ppt, f[1], maxh = 4.6)
      if (length(f) >= 3 && !is.na(f[3]) && nzchar(f[3]))
        ppt <- ph_with(ppt, block_list(fpar(ftext(f[3], fp_text(font.size = 12)))),
                       location = ph_location(left = 0.5, top = 6.2, width = 9, height = 1.1))
    }
  }
  out <- file.path(OUT, "E2F_scRNA_pilot_report.pptx")
  print(ppt, target = out); cat("Wrote", out, "(", length(ppt), "slides )\n")
}

build_html()
tryCatch(build_pptx(), error = function(e) message("pptx build failed: ", conditionMessage(e)))
cat("=== DONE make_report ===\n")
