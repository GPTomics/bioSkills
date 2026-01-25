# Figure Export Usage Guide

## Overview

This guide covers exporting publication-ready figures with proper resolution, sizing, and formatting.

## Prerequisites

```bash
# Python
pip install matplotlib seaborn

# R
install.packages(c('ggplot2', 'cowplot', 'patchwork'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Export my plot as a publication-ready PDF"
- "Save my figure at 300 DPI for journal submission"
- "Create a multi-panel figure with proper sizing"
- "Set up consistent styling for all my figures"

## Example Prompts

### Basic Export

> "Save my matplotlib figure as PDF and PNG at 300 DPI"

> "Export my ggplot as a 3.5 inch wide figure for a single-column journal"

### Journal Requirements

> "Format my figure for Nature submission (89mm single column width)"

> "Export my figure as TIFF with LZW compression for print publication"

### Styling

> "Set up Arial font and 8pt text for all my figures"

> "Create a consistent theme for my paper's figures"

### Multi-panel

> "Combine these 4 plots into a 2x2 figure with panel labels A, B, C, D"

> "Create a complex figure layout with different panel sizes"

## What the Agent Will Do

1. Set up figure dimensions and resolution
2. Configure fonts and styling
3. Export in requested format(s)
4. Add panel labels if multi-panel
5. Ensure colorblind-friendly palettes if requested

## Tips

- Vector formats (PDF, SVG) are preferred for line plots
- Use 300 DPI minimum for raster images
- Check journal guidelines for exact requirements
- Arial or Helvetica fonts are widely accepted
- Test figures at intended print size (not screen size)
- Use viridis or other colorblind-safe palettes
