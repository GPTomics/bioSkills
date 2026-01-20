# Color Palettes Usage Guide

## Overview

Proper color selection ensures figures are accessible, publication-ready, and effectively communicate data patterns.

## Palette Selection Guide

| Data Type | Palette Type | Examples |
|-----------|--------------|----------|
| Continuous (0 to max) | Sequential | viridis, Blues |
| Centered (-x to +x) | Diverging | RdBu, coolwarm |
| Categories | Qualitative | Set1, tab10, npg |

## Quick Start Prompts

- "Apply a colorblind-friendly palette to my plot"
- "Use a diverging color scheme for my heatmap"
- "Get the Nature journal color palette"
- "Create a custom palette matching my lab colors"

## Requirements

```r
# R
install.packages(c('viridis', 'RColorBrewer', 'ggsci', 'colorspace'))
```

```bash
# Python
pip install matplotlib seaborn
```

## Colorblind Considerations

- **Avoid**: Red-green combinations
- **Prefer**: viridis, cividis, Blue-Orange
- **Test**: Use colorblind simulators

## Journal Guidelines

Many journals have specific color requirements:
- Some require CMYK for print
- Ensure sufficient contrast
- Use consistent colors throughout figures

## Related Skills

- **data-visualization/ggplot2-fundamentals** - Apply palettes
- **data-visualization/heatmaps-clustering** - Heatmap colors
