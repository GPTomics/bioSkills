# ggplot2 Fundamentals - Usage Guide

## Overview

ggplot2 is the R reference for publication-quality static graphics, based on Wilkinson's grammar of graphics. A figure is composed as data + aesthetic mappings + one or more geometries + scales + facets + theme. The modern defaults for publication are theme_classic + panel.grid.off + Okabe-Ito categorical palette + cairo_pdf TrueType-embedded export.

## Prerequisites

```r
install.packages(c('ggplot2', 'scales', 'ggrepel', 'ggtext', 'viridis', 'scico', 'ggrastr'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a scatter with theme_classic, no grid, Okabe-Ito colors"
- "Faceted boxplot per tissue with free y-axis"
- "Add log10 transform to the y-axis with formatted breaks"
- "Programmatic plot from variable name passed as a string"
- "Save as 89mm cairo_pdf for Nature submission"

## Example Prompts

### Publication scatter

> "ggplot2 scatter with theme_classic baseline. Color by group using Okabe-Ito (#0072B2/#D55E00). Top/right spines removed. Save as 89mm cairo_pdf."

### Faceted box+jitter

> "Faceted boxplot of expression per tissue with overlaid jitter. Suppress boxplot outliers because jitter shows all points. Use Sheather-Jones bandwidth on violins."

### Programmatic from variable name

> "Build a plotting function that accepts the x-axis variable as a string and uses tidy evaluation (`aes(x = .data[[var]])`) rather than deprecated aes_string."

### Rich-text labels

> "Use ggtext to render axis title 'log<sub>2</sub> fold change' with proper subscripts and italics."

### Rasterized large-N scatter

> "Wrap geom_point in ggrastr::rasterise so 50000 points export as raster while axes stay vector."

## What the Agent Will Do

1. Set up base ggplot with explicit `aes()` mappings.
2. Add one or more geoms (geom_point, geom_boxplot, etc.) with explicit constant vs mapped parameters.
3. Configure scales - explicit color/fill palettes from CVD-safe sources; log/sqrt transforms; explicit limits.
4. Apply facets only when grouping requires it; use `scales = 'fixed'` unless free is justified.
5. Apply publication theme - `theme_classic` + panel grid removed + black axis/tick lines.
6. Add labels via `labs()`; rich text via ggtext if needed.
7. Save with `ggsave(..., device = cairo_pdf, units = 'mm')` at journal-spec width (89 or 183).

## Tips

- **Use `theme_classic` + `panel.grid = element_blank()`** as publication baseline. Remove top/right axis lines if Nature style: `theme(axis.line = element_line())`.

- **`device = cairo_pdf` for TrueType embedding.** Default ggsave can produce non-embedded fonts; journals reject.

- **`linewidth` replaces `size` for lines** in ggplot2 3.4+. `size` still works for points.

- **`aes_string` is deprecated.** Use `aes(x = .data[[var]])` for programmatic plots in ggplot2 3.0+.

- **Constants outside aes, mappings inside.** `geom_point(color = 'red')` for constant; `geom_point(aes(color = group))` for mapping. `aes(color = 'red')` makes 'red' a 1-level category.

- **`max.overlaps = Inf`** on ggrepel - default 10 silently drops labels.

- **`outlier.shape = NA`** on boxplot when jittering raw points; otherwise outliers render twice.

- **`facet_wrap(scales = 'fixed')` for cross-panel comparison** (default). Use `'free_y'` only when scales legitimately differ.

- **`units = 'mm'` explicit on ggsave.** Default is inches; Nature single column = 89 mm; double = 183 mm.

- **`ggrastr::rasterise(geom_point(...))`** for large N - vector axes + raster scatter is the publication compromise.

- **`ggtext::element_markdown()`** enables HTML in axis titles for subscripts and italics - much cleaner than `expression()`.

## Related Skills

- data-visualization/color-palettes - Palette selection
- data-visualization/multipanel-figures - patchwork composition
- data-visualization/distribution-plots - Box / violin / raincloud
- data-visualization/volcano-and-ma-plots - Volcano with ggrepel
- data-visualization/heatmaps-clustering - geom_tile and ComplexHeatmap
