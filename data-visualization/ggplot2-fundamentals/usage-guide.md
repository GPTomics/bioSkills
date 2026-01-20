# ggplot2 Fundamentals - Usage Guide

## Overview
ggplot2 is a declarative visualization system based on the Grammar of Graphics. Build publication-quality figures layer by layer.

## Prerequisites
```r
install.packages(c('ggplot2', 'patchwork', 'ggrepel', 'scales', 'RColorBrewer', 'viridis'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a volcano plot with labeled significant genes"
- "Make a multi-panel figure with panels A, B, C"
- "Set up a consistent theme for all my figures"

## Example Prompts
### Basic Plots
> "Create a scatter plot of gene expression vs significance"

> "Make a bar chart showing sample counts by condition"

### Customization
> "Apply a Nature-style theme to my plot"

> "Add a regression line with confidence interval"

### Publication Export
> "Export my figure at 300 DPI for journal submission"

> "Save as vector PDF for publication"

## What the Agent Will Do
1. Define data mapping with aesthetics (x, y, color, size)
2. Add appropriate geoms (points, lines, bars)
3. Customize scales, labels, and theme
4. Save at publication-quality resolution

## Grammar of Graphics

| Component | Description |
|-----------|-------------|
| Data | The dataset to visualize |
| Aesthetics | Mappings (x, y, color, size) |
| Geoms | Visual elements (points, lines, bars) |
| Scales | How data maps to aesthetics |
| Facets | Small multiples |
| Theme | Visual styling |

## Tips
- Always include axis labels with units
- Use colorblind-friendly palettes (viridis)
- Export at 300+ DPI for publication
- Use vector format (PDF) when possible
- Keep font size readable (8pt minimum)

## Related Skills
- **differential-expression/de-visualization** - Specialized plots
- **reporting/rmarkdown-reports** - Embed in reports
