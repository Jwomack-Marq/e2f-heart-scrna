# Methods book

A chapter-per-analysis companion to the Shiny app: what each tab of the browser
is showing, how it was computed, with which parameters, **why that way**, and what
it cannot be used to claim.

Written for a reader who wants to review the analysis rather than operate the app —
the app shows results, this explains them. Every chapter names the app tabs it
documents, so the two can be read side by side.

## Read it

```bash
docs/render.sh                       # -> docs/_site/index.html
python3 -m http.server -d docs/_site 8080
```

`docs/_site/` is a build artifact and is git-ignored, like `shiny_app/deck.html`.

## How it is built, and why that way

**Rendering never executes code.** Every chapter declares `engine: markdown` and the
project sets `execute: enabled: false`, so `quarto render` needs only the Quarto CLI —
no R, no Bioconductor, and no `app_data.rds`. `docs/Dockerfile` is a 250 MB
`debian:bookworm-slim` + Quarto image with no R installed at all, which makes the
guarantee structural: a stray executable chunk fails the render instead of quietly
working on one machine and not another.

**Figures and tables are exported once, offline.**

```bash
docs/export.sh                       # all of it (needs shiny_app/app_data.rds)
docs/export.sh --list                # print the export ids without running
docs/export.sh --only=fourgroup      # re-export just the ids matching a regex
```

`docs/export_assets.R` runs in `e2f-enrich:latest` with the bundle mounted
**read-only** and writes:

| path | what | in git |
|---|---|---|
| `docs/assets/*.svg` | vector figures | yes |
| `docs/assets/qc-*.png` | the five upstream QC figures, decoded from `app$figs` | yes |
| `docs/_generated/*.md` | markdown table fragments the chapters `{{< include >}}` | yes |
| `docs/_generated/_manifest.md` | every export, its chapter, and what it shows | yes |

Two rules the exporter exists to enforce:

1. **A figure in the book is the figure in the app.** The exporter loads the app's own
   plotting functions out of `shiny_app/app.R` — everything above the `# ---- UI ----`
   marker — and calls them directly, so `go_dotplot_gg()` here and in the running app
   are the same closure over the same bundle. Nothing is redrawn from a description.
   Where a panel is a server *reactive* rather than a callable function (there is no way
   to call `comp_plot` from outside the server), the logic is mirrored in the exporter's
   `MIRRORS` section, and each mirror carries the `app.R` line range it mirrors so the
   two can be diffed by eye. The UMAP is the one figure whose rendering genuinely
   differs: the app draws it with WebGL plotly, which has no vector export, so the book
   shows a ggplot of the same coordinates and the same palette.
2. **No number in the book is typed by hand.** Load-bearing tables are generated
   fragments, so a rebuilt bundle changes the book. This is the discipline
   `shiny_app/build_deck.R` already follows, for the same reason: prose drifts from
   data, and a methods document is where that drift becomes a wrong claim in front of
   an audience. Numbers that do appear inline in prose carry a `file:line` citation to
   the source that sets them.

If you would rather not run the exporter, every plot in the app already has PDF/SVG/PNG
download buttons (`register_fig()`), and every table a CSV/XLSX one (`register_dl()`) —
hand-export into the same paths and the chapters are unchanged.

## Checks

```bash
python3 tools/check_docs_coverage.py --list
```

Asserts that all 18 top-level app tabs are claimed by exactly one chapter, that no
chapter claims a tab that no longer exists, that every referenced asset is present, and
that no chapter has escaped `engine: markdown`. `docs/render.sh` runs it before
rendering. It is Python, not R, for the same reason the render image has no R: it must
run on a machine with neither R nor the bundle.

## Adding or changing a chapter

1. Add the `.qmd` to `docs/` with `engine: markdown` and a `tabs:` list in the front
   matter naming the app tabs it documents (exact titles from `app.R`'s `nav_panel()`).
2. Add it to `book.chapters` in `_quarto.yml`.
3. Add any new figure or table to `docs/export_assets.R` as a `step(...)` and re-run
   `docs/export.sh --only=<id>`.
4. `python3 tools/check_docs_coverage.py && docs/render.sh`.

## What is deliberately not here

The mechanistic cardiomyocyte cell-cycle model in [`model/`](../model/) is the
quantitative half of the project and has its own documentation (`model/MODEL.md`,
`model/RESULTS.md`) and its own deck. This book covers the descriptive half — the
scRNA-seq analysis behind the browser — and links to the model rather than restating it.
