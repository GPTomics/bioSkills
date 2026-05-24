# Circos Plots - Usage Guide

## Overview

Circular genome plots (Krzywinski 2009) render chromosomes around a circle with stacked data tracks and inter-chromosomal arc/chord links. circlize (R; Gu 2014) is the modern R standard; pyCirclize (Python) is actively developed; Circos (Perl) remains the most powerful. The first decision is whether the circular layout adds information - Cleveland-McGill 1984 (replicated by Heer-Bostock 2010 on Mechanical Turk) establishes that position-on-common-scale (Cartesian) is the most accurate visual channel; circular position requires mental unwrapping and impairs value comparison. Use circular only when chromosome adjacency or chord-diagram pairwise relationships matter.

## Prerequisites

```r
install.packages('circlize')
```

```bash
pip install pycirclize
# CLI (Perl):
brew install circos    # or conda install -c bioconda circos
```

## Quick Start

Tell your AI agent what you want to do:
- "Circos plot of CNV log2 across all autosomes with red gain / blue loss / grey neutral"
- "Inter-chromosomal SV link diagram with chord arcs colored by variant type"
- "Karyotype with cytoband ideograms and a variant density histogram outer track"
- "Hi-C inter-chromosomal contact summary as circos heatmap"

## Example Prompts

### CNV summary circos

> "Circos with hg38 ideograms (start at 12 o'clock), gene density histogram, variant density heatmap, CNV log2 scatter (red >0.3, blue <-0.3, grey otherwise). End with `circos.clear()`."

### SV link diagram

> "Inter-chromosomal SV chord plot. Color links by SV type (deletion, duplication, inversion, translocation). Filter to >100 supporting reads."

### Hi-C summary

> "Chromosome-level Hi-C interaction summary using circos heatmap track. Color by log2 contact frequency."

### Pre-flight critique

> "Before building this circos, is the data better suited to a heatmap? Apply the Cleveland-McGill 1984 effectiveness ranking."

## What the Agent Will Do

1. Confirm circular layout is appropriate (chromosome adjacency OR chord diagram). If not, suggest heatmap/linear alternative.
2. Initialize ideogram via `circos.initializeWithIdeogram` (matching species/build).
3. Set explicit chromosome.index for correct ordering (chr1...22, X, Y).
4. Build data tracks: histogram, heatmap, scatter, lines, in inside-to-outside order.
5. Add chord links with color/width encoding interaction strength.
6. Use CVD-safe palette (Crameri vik for diverging tracks, Okabe-Ito for categorical).
7. End every plot block with `circos.clear()` to reset global state.
8. Export PDF / PNG; verify chromosome rotation matches publication convention.

## Tips

- **Justify circular layout.** Heatmap usually wins (Cleveland-McGill 1984). Use circular only for adjacency-meaningful contexts.

- **Always end with `circos.clear()`** - circos.par settings are global; persist into next plot if not cleared.

- **Explicit `chromosome.index`** for correct order (chr1...22, X, Y); default is input-file order.

- **Match species to data build.** `species = 'hg38'` for hg38 coordinates; mismatch produces wrong cytoband display.

- **Inside-to-outside track building.** First `circos.genomicTrack` is outermost; subsequent calls add inward.

- **Filter dense link sets** - >5000 chord arcs produce a black blob. Aggregate or threshold.

- **`gap.degree` for visual breaks** - bigger gap before chr1 (or between autosomes and sex chromosomes) helps reader orient.

- **`start.degree = 90`** puts chr1 at 12 o'clock (most common convention).

- **For per-chromosome CNV** consider linear karyoploteR (`copy-number/cnv-visualization`) - often clearer than circos.

- **pyCirclize (Python) is younger but actively developed** - Pythonic API; production-ready as of 2024.

## Related Skills

- copy-number/cnv-visualization - karyoploteR linear alternative
- variant-calling/structural-variant-calling - SV data for links
- hi-c-analysis/hic-visualization - Hi-C circular contact display
- data-visualization/genome-tracks - Linear genome track alternative
- data-visualization/color-palettes - Track and link palettes
