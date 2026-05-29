# Interactive Visualization - Usage Guide

## Overview

Interactive HTML/web plots support exploration via zoom, pan, hover, and brushing. plotly is the modern Python and R default (ggplotly converts ggplot2 directly); bokeh suits server-side dashboards; gganimate / plotly frames cover animation. The static-export pipeline changed materially in 2025 - orca is EOL; Kaleido v1+ is the standard but no longer bundles Chrome and dropped EPS. Always produce static alongside interactive for journal submission.

## Prerequisites

```bash
pip install plotly bokeh kaleido altair
# Static export requires installed Chrome/Chromium
```

```r
install.packages(c('plotly', 'gganimate', 'htmlwidgets', 'DT', 'leaflet'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Make this scatter interactive with plotly - hover shows sample ID"
- "Convert this ggplot to plotly via ggplotly"
- "WebGL acceleration for 20000-point scatter"
- "Static PDF export via Kaleido v1; verify Chrome is installed"
- "Animate time-course UMAP with plotly frames"

## Example Prompts

### Interactive volcano

> "Build interactive volcano with plotly Express; hover_data = ['gene', 'baseMean', 'padj']. Reference lines at fdr=0.05 and lfc=1. WebGL for >5000 points."

### ggplotly conversion

> "Take this ggplot and pass to ggplotly with tooltip = c('text', 'x', 'y'). Save as self-contained HTML via htmlwidgets."

### Static export

> "Export the plotly figure as a 300 DPI PDF for the paper. Use Kaleido v1; verify Chrome is installed; produce a fallback PNG if PDF fails."

### Animated UMAP time-course

> "Animate UMAP across 5 timepoints using plotly animation_frame='timepoint'. Constant axis range across frames."

### Dashboard

> "Build a Bokeh server dashboard with linked-brushing scatter + bar chart for cohort exploration."

## What the Agent Will Do

1. Decide interactive vs static: notebook exploration / supplement / dashboard -> interactive; journal figure -> static (alongside).
2. For Python: plotly Express for fast onboarding; graph_objects for fine control; Scattergl for >5000 points.
3. For R: prefer ggplotly for ggplot conversion; htmlwidgets ecosystem for tables/maps/networks.
4. Configure hover_data, tooltips, color/size aesthetics.
5. For static export: confirm Kaleido v1 + Chrome installed; use `fig.write_image()` with no `engine=`.
6. For animation: limit to ≤100 frames; pre-aggregate per-frame data; constant axis range.
7. Save HTML as self-contained (`htmlwidgets::saveWidget(..., selfcontained = TRUE)`) for portability.
8. Always produce static PDF alongside HTML for journal submission.

## Tips

- **Always produce static alongside interactive.** Interactive HTML cannot be cited as a paper figure.

- **Kaleido v1 is the current static-export engine.** No `engine=` argument needed; orca is EOL.

- **Kaleido v1 requires installed Chrome/Chromium.** No longer bundled. Verify with `which chrome` or equivalent.

- **EPS export dropped in Kaleido v1.** Export PDF then convert via `pdf2ps` (ghostscript). For pure-vector EPS, may need a separate path.

- **WebGL (Scattergl) for >5000 points.** Default Scatter is slow at scale.

- **HTML > 10 MB is a warning sign.** Pre-aggregate (Datashader) or filter before exporting.

- **ggplotly is the lowest-friction R interactive path** - write ggplot, get plotly.

- **Self-contained HTML** for portability (`selfcontained = TRUE` in saveWidget). Default CDN-linked breaks offline.

- **For animation, cap at ~100 frames.** File size and viewer attention both degrade.

- **bokeh for dashboards / streaming;** plotly for static-figure-equivalent interactivity.

- **Don't cite interactive in figure captions.** Treat as supplement; static for the paper.

## Related Skills

- reporting/quarto-reports - Embed in scientific reports
- reporting/rmarkdown-reports - htmlwidgets in Rmd
- data-visualization/ggplot2-fundamentals - ggplot for ggplotly conversion
- data-visualization/dimensionality-reduction-plots - Interactive UMAP/PCA
- data-visualization/network-visualization - PyVis interactive networks
