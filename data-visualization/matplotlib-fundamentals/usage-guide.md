# matplotlib Fundamentals - Usage Guide

## Overview

matplotlib is the foundational Python plotting library. For publication-quality figures, use the **object-oriented Figure/Axes API** (not pyplot state-machine), `constrained_layout` for axis alignment, `pdf.fonttype=42` for TrueType font embedding (mandatory for Nature/IEEE/ACM), and CVD-safe palettes. seaborn integrates by sharing the same Figure/Axes objects.

## Prerequisites

```bash
pip install matplotlib seaborn numpy pandas
# CVD-safe Crameri scientific colormaps:
pip install cmcrameri
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a scatter plot using the Figure/Axes API with rasterized points and TrueType fonts"
- "Make a 2x3 subplot grid with constrained_layout"
- "Set up rcParams for Nature single-column figures (89mm width, 7pt body text)"
- "Use Okabe-Ito categorical palette and Crameri batlow colormap"
- "Embed TrueType fonts so the PDF passes journal submission"

## Example Prompts

### Publication-ready scatter

> "Volcano-style scatter at 89mm width with rasterized points, Okabe-Ito Up/Down/NS colors, top/right spines removed, TrueType-embedded PDF."

### Multi-panel layout

> "2x2 subplot grid with constrained_layout, panel labels a/b/c/d in bold upper-left of each subplot, shared y-axis for the left column."

### Sequential heatmap

> "imshow heatmap with Crameri batlow colormap, colorbar shrunk to 60% of axes height, vector axis text + raster cells."

### seaborn integration

> "Use sns.scatterplot with ax= argument to render into a pre-created matplotlib Axes."

### Migrate from pyplot

> "Refactor this pyplot script (plt.scatter, plt.xlabel) to the OO Figure/Axes API."

## What the Agent Will Do

1. Set rcParams at the top of the script for `pdf.fonttype=42`, font sizes (Nature 5-7pt), spine widths.
2. Build figures via `fig, ax = plt.subplots(figsize=(mm/25.4, mm/25.4), constrained_layout=True)`.
3. Render via Axes methods: `ax.scatter`, `ax.plot`, `ax.bar`, `ax.imshow`.
4. Add `rasterized=True` to scatter / imshow for >1000 elements.
5. Apply CVD-safe palettes from Okabe-Ito (categorical) or Crameri / viridis (continuous).
6. Set explicit limits, ticks, spines; remove top/right by default.
7. Save with `fig.savefig('out.pdf', dpi=300, bbox_inches='tight')`.

## Tips

- **Use the OO Figure/Axes API**, not pyplot state-machine. `fig, ax = plt.subplots(); ax.scatter(...)` not `plt.scatter(...)`.

- **`pdf.fonttype=42` is mandatory** for many journals. Type-3 fonts (default) are not searchable and rejected at submission.

- **`constrained_layout=True` instead of tight_layout** in matplotlib 3.6+. Handles colorbars and shared axes correctly.

- **`rasterized=True` for >1000 points** in scatter / imshow. Keep axes and text vector.

- **figsize in inches.** For Nature 89mm: `figsize=(89/25.4, 70/25.4)`. For double-column 183mm: `figsize=(183/25.4, height/25.4)`.

- **Top/right spines off** is the modern convention: `ax.spines[['top','right']].set_visible(False)`.

- **seaborn returns FacetGrid for figure-level functions** (`displot`, `catplot`, `relplot`). Use `.set_axis_labels()` not `.set_xlabel()`. Axes-level functions return Axes (use `ax=` argument).

- **Colorbar control**: `plt.colorbar(im, ax=ax, shrink=0.6, aspect=20)` to avoid colorbar dominating subplot.

- **Okabe-Ito for categorical** (CVD-safe 8 colors). Crameri batlow / viridis for sequential. RdBu_r or vik for diverging.

- **Save format**: PDF for vector axes + raster scatter (best journal compromise); PNG at 300 DPI for raster only; SVG for editable.

## Related Skills

- data-visualization/color-palettes - Palette selection
- data-visualization/multipanel-figures - GridSpec and complex layouts
- data-visualization/distribution-plots - seaborn box/violin/raincloud
- data-visualization/heatmaps-clustering - seaborn.clustermap
- data-visualization/volcano-and-ma-plots - matplotlib scatter for volcano
- reporting/figure-export - DPI / format / journal specs
