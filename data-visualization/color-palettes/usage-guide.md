# Color Palettes for Scientific Visualization - Usage Guide

## Overview

Selecting a color palette is a perceptual decision, not an aesthetic one. The three modern criteria are perceptual uniformity (equal data steps look equal), color-vision-deficiency safety (~6% of males have red-green deficiency), and grayscale monotonicity (luminance increases monotonically along the colormap). The Crameri scientific colormaps, viridis family, cividis, and Okabe-Ito categorical palette satisfy these by construction. Rainbow / jet and most red-green pairings do not.

## Prerequisites

```r
install.packages(c('viridis', 'RColorBrewer', 'scico', 'khroma', 'ggsci', 'colorspace'))
```

```bash
pip install matplotlib colorcet cmcrameri colorspacious
```

## Quick Start

Tell your AI agent what you want to do:
- "Pick a sequential colormap for an expression heatmap"
- "Pick a diverging colormap for log-fold-change"
- "Pick a categorical palette for 6 cell types that is CVD-safe"
- "Replace this jet colormap with a perceptually-uniform alternative"
- "Verify this palette stays readable under deuteranopia"

## Example Prompts

### Heatmap palette (sequential)

> "Pick a sequential perceptually-uniform colormap for an expression heatmap. Verify it remains interpretable when converted to grayscale."

### Log-fold-change (diverging)

> "Use a diverging Crameri vik palette for a log-fold-change heatmap with symmetric bounds at +/- the 99th percentile of |LFC|. Zero must map to pure white."

### Cell-type categorical (CVD-safe)

> "Assign one CVD-safe Okabe-Ito color to each of these 7 cell types. Reserve grey for ambient/unassigned."

### Migration from jet

> "Find every plot in this notebook that uses cmap='jet' or 'rainbow' and replace with viridis (sequential), turbo (jet-like but perceptually uniform), or batlow (Crameri)."

### CVD audit

> "Run colorspace::cvd_emulator on the current palette and report whether the categories remain distinguishable under deuteranopia and protanopia."

## What the Agent Will Do

1. Identify whether the data axis is sequential (single direction), diverging (signed around midpoint), cyclic, or categorical.
2. Select from a CVD-safe shortlist: viridis/cividis/magma/batlow for sequential, vik/roma/RdBu for diverging, romaO for cyclic, Okabe-Ito for ≤8 categorical.
3. For diverging data, set symmetric bounds (`vmin = -max(|data|), vmax = +max(|data|)`) and a pure-white midpoint.
4. For categorical data, cap at 8 colors before switching to facet/shape/aggregation.
5. Apply the palette via the package-native scale function (`scale_color_scico`, `scale_color_viridis_c`, `scale_color_manual`).
6. Perform the grayscale test (convert PNG to grayscale; verify data order is still readable) before publication.
7. Run CVD simulation (colorspace::cvd_emulator or colorspacious) and report.

## Tips

- **Default to viridis or cividis for sequential, vik or RdBu for diverging.** All four are perceptually uniform and CVD-safe.

- **cividis is CVD-optimal** - designed to look near-identical to normal and CVD viewers (Nuñez 2018). Use for any figure where accessibility matters.

- **Okabe-Ito is the canonical 8-color CVD-safe categorical palette** (Wong 2011). Memorize the hexes: `#E69F00 #56B4E9 #009E73 #F0E442 #0072B2 #D55E00 #CC79A7 #000000`.

- **Crameri scientific colormaps** (batlow, lipari, vik, roma, bam, plus cyclic *O variants) are the modern publication standard (Crameri 2020 *Nat Commun*). Available via R `scico`, Python `cmcrameri`.

- **Symmetric bounds for diverging data.** `vmin = -vmax` ALWAYS. Asymmetric bounds misalign zero with the white midpoint and mis-encode signed data.

- **Pure white at the diverging midpoint, not light gray.** Gray reads as "weak signal"; white reads as "exactly zero."

- **The grayscale test is the most actionable CVD check.** Save figure as PNG, desaturate. If the data order is still readable, the colormap is luminance-monotonic and accessible.

- **Rainbow / jet is harmful** (Borland 2007, Crameri 2020). Non-monotonic luminance creates artifactual bands. Migrate to viridis or turbo.

- **Categorical palettes top out at ~8-10 distinguishable hues.** Above this, switch to faceting, aggregation, or color+shape combinations.

- **Journal-brand palettes (ggsci npg/aaas/lancet/jama)** are stylistic, not accessible. Use for compliance with house style only.

- **Custom palette construction**: pick endpoints with similar luminance; pass through pure white for diverging; verify with grayscale test before committing.

## Related Skills

- data-visualization/heatmaps-clustering - Robust diverging bounds for heatmaps
- data-visualization/volcano-and-ma-plots - Okabe-Ito Up/Down/NS conventions
- data-visualization/ggplot2-fundamentals - Applying palettes in ggplot2 scales
- data-visualization/dimensionality-reduction-plots - Categorical palette for cluster labels
