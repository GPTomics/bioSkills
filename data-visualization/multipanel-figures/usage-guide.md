# Multi-Panel Figure Assembly Usage Guide

## Overview

Combine individual ggplot2 plots into publication-ready multi-panel figures with consistent styling, shared legends, and proper annotations.

## Quick Start Prompts

- "Combine my volcano plot and PCA into one figure"
- "Create a 2x2 panel figure with labels A, B, C, D"
- "Share the legend across all panels"
- "Add an inset plot to my main figure"

## Tool Comparison

| Package | Strengths |
|---------|-----------|
| patchwork | Intuitive operators, ggplot2 native |
| cowplot | Simple API, good for basic layouts |
| gridExtra | Fine control, complex arrangements |

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

## Requirements

```r
install.packages(c('patchwork', 'cowplot', 'gridExtra'))
```

## Related Skills

- **data-visualization/ggplot2-fundamentals** - Create individual plots
- **reporting/rmarkdown-reports** - Embed figures
