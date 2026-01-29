# Multi-Panel Figures - Usage Guide

## Overview
Combine individual ggplot2 plots into publication-ready multi-panel figures with consistent styling, shared legends, and proper annotations.

## Prerequisites
```r
install.packages(c('patchwork', 'cowplot', 'gridExtra'))
```

```bash
pip install matplotlib numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Combine my volcano plot and PCA into one figure"
- "Create a 2x2 panel figure with labels A, B, C, D"
- "Share the legend across all panels"

## Example Prompts
### Basic Assembly
> "Combine these three plots into a horizontal layout"

> "Stack my plots vertically with panel labels"

### Complex Layouts
> "Create a figure with a large plot on left and two small plots stacked on right"

> "Make a 2x2 grid where the top row spans both columns"

### Shared Elements
> "Share the legend across all panels"

> "Add an inset plot to my main figure"

## What the Agent Will Do
1. Arrange individual plots using patchwork operators
2. Apply panel labels (A, B, C, D)
3. Adjust relative widths/heights as needed
4. Collect and position shared legends
5. Export at publication resolution

## Tool Comparison

| Package | Language | Strengths |
|---------|----------|-----------|
| patchwork | R | Intuitive operators, ggplot2 native |
| cowplot | R | Simple API, good for basic layouts |
| gridExtra | R | Fine control, complex arrangements |
| matplotlib GridSpec | Python | Full control, Python workflows |

## Common Layouts
```r
# 1x3 horizontal
p1 + p2 + p3

# 3x1 vertical
p1 / p2 / p3

# 2x2 grid
(p1 + p2) / (p3 + p4)

# Wide left, narrow right
p1 + p2 + plot_layout(widths = c(2, 1))
```

## Tips
- Use patchwork for most R layouts (simplest syntax)
- Add panel labels with `plot_annotation(tag_levels = 'A')` in R
- Control relative sizes with `plot_layout(widths = ..., heights = ...)` in R
- Collect legends with `plot_layout(guides = 'collect')` in R
- Use GridSpec for Python workflows with full layout control
- Standard sizes: 7" single column, 14" double column for journals

## Related Skills
- **data-visualization/ggplot2-fundamentals** - Create individual plots
- **reporting/rmarkdown-reports** - Embed figures
