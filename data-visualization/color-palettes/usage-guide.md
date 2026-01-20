# Color Palettes - Usage Guide

## Overview
Proper color selection ensures figures are accessible, publication-ready, and effectively communicate data patterns.

## Prerequisites
```r
# R
install.packages(c('viridis', 'RColorBrewer', 'ggsci', 'colorspace'))
```

```bash
# Python
pip install matplotlib seaborn
```

## Quick Start
Tell your AI agent what you want to do:
- "Apply a colorblind-friendly palette to my plot"
- "Use a diverging color scheme for my heatmap"
- "Get the Nature journal color palette"

## Example Prompts
### Colorblind-Friendly
> "Apply viridis colors to my heatmap"

> "Use a colorblind-safe palette for my categorical data"

### Journal Styles
> "Use the Nature journal color palette for my figure"

> "Apply Cell-style colors to my plot"

### Custom Palettes
> "Create a custom palette matching my lab colors (#1A5276, #F39C12, #27AE60)"

> "Generate a gradient from blue to red for my heatmap"

### Data-Appropriate
> "Use a diverging palette centered at zero for my log2FC data"

> "Apply a sequential palette for my p-value heatmap"

## What the Agent Will Do
1. Identify the data type (continuous, diverging, categorical)
2. Select appropriate palette type
3. Apply palette to the visualization
4. Ensure accessibility (colorblind-friendly)
5. Test visual contrast

## Palette Selection Guide

| Data Type | Palette Type | Examples |
|-----------|--------------|----------|
| Continuous (0 to max) | Sequential | viridis, Blues |
| Centered (-x to +x) | Diverging | RdBu, coolwarm |
| Categories | Qualitative | Set1, tab10, npg |

## Tips
- Avoid red-green combinations (colorblind unfriendly)
- Use viridis or cividis for safe continuous palettes
- Test with colorblind simulators
- Use consistent colors throughout all figures in a paper
- Check journal requirements (some require CMYK for print)

## Related Skills
- **data-visualization/ggplot2-fundamentals** - Apply palettes
- **data-visualization/heatmaps-clustering** - Heatmap colors
