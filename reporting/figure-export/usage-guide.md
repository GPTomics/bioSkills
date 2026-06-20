# Figure Export Usage Guide

## Overview

This guide covers exporting figures that survive journal production. The central idea is that a figure is four layers - vector structure, a dense raster data layer, editable type, and a color encoding - and export means giving each the representation it needs. The decisions that actually matter are: rasterize the dense layer (so a million-point scatter does not produce an unopenable vector file), keep fonts as editable text, choose colors that survive the print RGB-to-CMYK conversion and grayscale, and design at the journal's exact width.

## Prerequisites

```bash
# Python
pip install matplotlib numpy

# R
install.packages(c('ggplot2', 'ggrastr', 'ragg', 'svglite'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Export this UMAP without making a 200 MB PDF the typesetter can't open"
- "Save my figure for Nature at 89 mm with editable fonts"
- "My PDF text is uneditable in Illustrator - fix the font type"
- "Export a TIFF with LZW compression for Cell"
- "Pick a colorblind-safe palette that also survives grayscale printing"

## Example Prompts

### Hybrid (Vector + Raster)

> "Rasterize the scatter points but keep the axes and labels vector so the file opens"

> "Export my single-cell UMAP as a PDF that production can actually open"

### Fonts and Formats

> "Make the text in my matplotlib PDF selectable and editable"

> "Export my ggplot as a 174 mm two-column figure with embedded fonts"

### Journal Requirements

> "Format my figure for Nature submission at 89 mm single-column width, RGB"

> "Export as flattened 8-bit TIFF with LZW, no alpha, for PLOS"

### Color and Reproducibility

> "Replace the rainbow heatmap with a perceptually uniform map that reads in grayscale"

> "Make the exported PDF byte-stable so it does not show up as changed in git every run"

## What the Agent Will Do

1. Set publication rcParams/theme: TrueType fonts (`pdf.fonttype=42`), constrained layout, target font sizes
2. Design the figure at the journal's exact physical width in mm/inches
3. Rasterize dense data layers (`rasterized=True` / `ggrastr::rasterise`) while keeping structure and text vector
4. Choose a perceptually uniform or CVD-safe palette and add redundant encoding
5. Export the right format for the target (vector PDF/EPS for Nature/Science, LZW TIFF for Cell/PLOS), with byte-stable metadata

## Tips

- DPI is meaningless for pure vector; physical figure size is the real control. Design to size, do not rescale.
- The default matplotlib `pdf.fonttype` is Type 3 (uneditable). Set 42 once, globally, to avoid revise-and-resubmit over fonts.
- EPS cannot hold transparency - it flattens alpha to opaque. Use PDF/SVG, or rasterize the transparent layer.
- Saturated RGB blues/greens shift muddy after the journal's CMYK conversion; pull saturation back and submit RGB.
- Avoid `bbox_inches='tight'` for camera-ready figures - it makes output size depend on font metrics. Use `constrained_layout` instead.
- Never JPEG a figure with text or lines; block compression rings the edges.
- Inspect the saved file at 100%, not the notebook inline preview - the exported file is the deliverable.

## Related Skills

- data-visualization/ggplot2-fundamentals - Building the plots in R
- data-visualization/multipanel-figures - Composing multi-panel layouts
- data-visualization/color-palettes - Choosing perceptual and CVD-safe palettes
- reporting/publication-tables - The table counterpart to figure export
