# ggplot2 Fundamentals Usage Guide

## Overview

ggplot2 is a declarative visualization system based on the Grammar of Graphics. Build publication-quality figures layer by layer.

## Quick Start Prompts

- "Create a volcano plot with labeled significant genes"
- "Make a multi-panel figure with panels A, B, C"
- "Set up a consistent theme for all my figures"
- "Export my figure at 300 DPI for publication"

## Grammar of Graphics

| Component | Description |
|-----------|-------------|
| Data | The dataset to visualize |
| Aesthetics | Mappings (x, y, color, size) |
| Geoms | Visual elements (points, lines, bars) |
| Scales | How data maps to aesthetics |
| Facets | Small multiples |
| Theme | Visual styling |

## Workflow

1. **Define data** - Clean dataframe
2. **Map aesthetics** - x, y, color, fill
3. **Add geoms** - Points, lines, bars
4. **Customize** - Scales, labels, theme
5. **Save** - ggsave at appropriate DPI

## Requirements

```r
install.packages(c('ggplot2', 'patchwork', 'ggrepel', 'scales', 'RColorBrewer', 'viridis'))
```

## Publication Checklist

- [ ] Font size readable (≥8pt)
- [ ] Axis labels with units
- [ ] Legend positioned appropriately
- [ ] Colors accessible (colorblind-friendly)
- [ ] Resolution ≥300 DPI
- [ ] Vector format (PDF) when possible

## Related Skills

- **differential-expression/de-visualization** - Specialized plots
- **reporting/rmarkdown-reports** - Embed in reports
