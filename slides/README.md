# Meeting decks

Quarto → PowerPoint. Rendering executes nothing (`engine: markdown`), so it needs only
the Quarto CLI — the same `e2f-docs:latest` image the methods book uses.

```bash
docs/export.sh                  # (re)generate slides/assets/ from the data bundle
slides/render.sh                # -> slides/qc-and-processing-review.pptx (+ .pdf preview)
slides/render.sh other.qmd      # render a different deck
```

## Decks

| File | For |
|---|---|
| `qc-and-processing-review.qmd` | 2026-08-26 meeting: processing, QC and thresholds, and the PIP-seq vs 10x platform question. Speaker notes on nearly every slide. |

## Two things to know before editing

**Figures must be PNG, and they come from the exporter.** Pandoc's pptx writer cannot be
trusted with SVG, and no SVG→PNG converter exists on this machine. So
`docs/export_assets.R` writes a raster twin of every methods-book figure into
`slides/assets/` (via `SLIDE_PNG_DIR`) — rendered from the same ggplot object through a
second device, not converted after the fact. One plotting source, two formats, nothing to
drift. `render.sh` fails before calling Quarto if a referenced asset is missing.

**One content block per slide, or pandoc splits it.** A slide that mixes a table or an
image with a paragraph becomes two slides, the second untitled and orphaned. Put prose in
`::: {.notes}`, and where text and a figure belong together use two columns:

```markdown
:::: {.columns}
::: {.column width="42%"}
- the point
:::
::: {.column width="58%"}
![](assets/figure.png)
:::
::::
```

Check for orphans after rendering — 22 slides in, 22 slides out:

```bash
unzip -l slides/qc-and-processing-review.pptx | grep -c 'ppt/slides/slide[0-9]*\.xml'
```

The `.pptx` and `.pdf` are build artifacts and are git-ignored; the `.qmd` is the source.
