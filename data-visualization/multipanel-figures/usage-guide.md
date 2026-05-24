# Multi-Panel Figures - Usage Guide

## Overview

Multi-panel figures combine 2+ subplots into a single composed figure with shared elements (legends, axes) and panel labels (a, b, c) in journal convention. patchwork is the modern R default - version 1.2.0+ added `axes='collect'` for shared-axis collection. matplotlib subfigures (3.4+) is the modern Python approach. cowplot remains useful for edge cases where patchwork's alignment struggles.

## Prerequisites

```r
install.packages(c('patchwork', 'cowplot', 'gridExtra'))
# patchwork >= 1.2.0 required for axes='collect'
```

```bash
pip install matplotlib
# subfigures stable since matplotlib 3.4
```

## Quick Start

Tell your AI agent what you want to do:
- "Combine 4 ggplot objects in a 2x2 grid with shared legend, collected axes, lowercase panel labels"
- "matplotlib GridSpec layout: top row 2 panels, bottom row 1 spanning panel"
- "Sized at Nature single-column (89mm) with cairo_pdf TrueType embedding"
- "Inset a small plot in the upper-right corner of a larger plot"
- "Custom layout via patchwork design string"

## Example Prompts

### 2x2 grid with shared legend (patchwork)

> "Combine p1-p4 into a 2x2 grid using patchwork. Collect guides and axes; lowercase bold panel labels in upper-left. Save as 180mm double-column cairo_pdf."

### Top-2, bottom-1 layout

> "Use patchwork `(p1 | p2) / p3` for a top row of two panels and a single bottom panel spanning the full width."

### matplotlib GridSpec

> "matplotlib GridSpec with 2 rows, 3 columns. ax1 in [0,0], ax2 spans [0, 1:], ax3 spans [1, :]. constrained_layout for axis alignment."

### Inset

> "Add a small PCA plot as inset in the upper-right corner of the main scatter using patchwork::inset_element."

### Subfigures with colorbar

> "matplotlib subfigures: left subfigure has 2 stacked scatter panels, right subfigure has one heatmap with colorbar. constrained_layout for each subfigure independently."

## What the Agent Will Do

1. Identify subplot composition (rows, columns, spans, ratios).
2. For R: use patchwork operators (`+`, `/`, `|`) for simple grids; `plot_layout(design=...)` for complex.
3. For Python: GridSpec for fine control; subfigures for nested layouts.
4. Add panel labels in Nature (lowercase bold) or Cell (uppercase bold) convention.
5. Collect shared legends with `guides='collect'`; share axes with `axes='collect'` (patchwork ≥ 1.2.0).
6. Set sizing at journal specs (89 mm single col, 183 mm double; mm units explicit).
7. Save with `device = cairo_pdf` (R) or `pdf.fonttype=42` (Python) for TrueType.

## Tips

- **patchwork ≥ 1.2.0 for `axes='collect'`.** Verify with `packageVersion('patchwork')`. Older versions silently ignore.

- **`device = cairo_pdf` on ggsave** for TrueType font embedding. Default device can produce non-portable PDF.

- **`units = 'mm'` explicit on ggsave.** Nature single col = 89; double = 183; max height = 247.

- **Lowercase bold for Nature, uppercase bold for Cell** panel labels. `plot_annotation(tag_levels = 'a')` for lowercase; `'A'` for uppercase.

- **Panel label position consistency** requires standardized y-axis label widths across subplots. If they differ, either pad with whitespace or position tags inside the plotting area.

- **Shared legend requires identical scales** across subplots. Different `scale_color_manual` calls produce duplicate (not merged) legends.

- **For complex layouts**, patchwork's design string (`"AAB\nAAB\nCCC"`) is more readable than nested operators.

- **matplotlib `constrained_layout=True`** is the modern default (3.6+); handles colorbars correctly. `tight_layout()` does not.

- **subfigures (matplotlib 3.4+)** are stronger than nested GridSpec for complex composition because each subfigure has its own layout engine.

- **Inset with `patchwork::inset_element(p_inset, left, bottom, right, top)`** in normalized [0, 1] coordinates of the parent plot.

## Related Skills

- data-visualization/ggplot2-fundamentals - Individual ggplot subpanels
- data-visualization/matplotlib-fundamentals - Python equivalent
- reporting/figure-export - DPI / format / journal-spec
- data-visualization/color-palettes - Consistent palette across panels
