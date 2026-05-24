# Distribution Plots - Usage Guide

## Overview

Distribution plots compare continuous data across a small number of groups. The default `bar of mean ± SEM` collapses distributional information and misleads (Weissgerber 2015). Modern alternatives are raincloud (Allen 2019) for N 30-200, beeswarm/quasirandom for N 10-30, raw points for N < 10, and letter-value plots (Hofmann 2017) or violin with explicit bandwidth for large N.

## Prerequisites

```r
install.packages(c('ggplot2', 'ggbeeswarm', 'ggdist', 'gghalves', 'lvplot'))
```

```bash
pip install seaborn matplotlib ptitprince
```

## Quick Start

Tell your AI agent what you want to do:
- "Replace this bar plot with a raincloud showing each individual point"
- "Plot expression per cluster with quasirandom jitter for N=15 per group"
- "Make a split violin comparing condition within each cell type"
- "Use Sheather-Jones bandwidth so bimodality is visible"
- "Add N to each x-axis tick label"

## Example Prompts

### Raincloud for N 30-200

> "Build a horizontal raincloud plot comparing biomarker level across 3 treatment arms (n=80 each). Show every point on one side, box+median in the middle, half-violin on the other side."

### Beeswarm for small N

> "Quasirandom beeswarm of expression across 6 cell types (n=15 each) with a horizontal crossbar for the median per group."

### Single-cell split violin

> "Split violin of gene X expression per cluster, with control vs treatment on the two halves. Use Sheather-Jones bandwidth so bimodal expression is visible."

### Letter-value plot for large N

> "Letter-value (boxen) plot for N=2000 per group instead of standard boxplot - preserves tail structure."

### Migration from bar plot

> "Convert this bar-of-mean ± SEM figure to a raincloud. Keep the same color scheme. Add N per group."

## What the Agent Will Do

1. Inspect N per group; pick encoding from the N decision table.
2. Default to raw-point overlay for N < 30; raincloud for N 30-200; letter-value for N > 200.
3. For violin: use explicit `bw = 'SJ'` (Sheather-Jones); not the default Silverman.
4. For box + jitter: set `outlier.shape = NA` to avoid double-rendering outliers.
5. Annotate N per group (x-axis tick label or text stratum).
6. Color via Okabe-Ito or other CVD-safe palette.
7. Apply statistical annotation via `ggpubr::stat_compare_means` or `statannotations` if requested (see `statistical-annotation` skill).

## Tips

- **Replace bar plots with raincloud or beeswarm.** Weissgerber 2015 documented that bars dominate top physiology journals despite collapsing distributional information.

- **Show every point for N < 30.** Aggregate summaries hide too much at small N.

- **Use Sheather-Jones bandwidth on violins** (`bw = 'SJ'`). Default Silverman oversmooths bimodality.

- **`trim = FALSE` on violin** so KDE doesn't appear to extend beyond data range.

- **`outlier.shape = NA` on boxplot when jittering** to avoid double-rendering.

- **Notched boxplots need N ≥ 15.** Below this notches can extend beyond Q1/Q3 and look broken.

- **Annotate N.** Add to x-tick label (`Control (n=12)`) or as text stratum.

- **For paired comparisons in 2 conditions, use split violin** (gghalves, introdataviz) - more compact than dodged violins.

- **Letter-value plot (`geom_lv` / `boxenplot`) for large N** - preserves tail structure that boxplot collapses (Hofmann 2017).

- **Quasirandom jitter is deterministic** (van der Corput sequence). Reruns produce identical layouts; standard `geom_jitter` is stochastic.

## Related Skills

- data-visualization/statistical-annotation - p-value brackets between groups
- data-visualization/color-palettes - CVD-safe categorical palettes
- data-visualization/ggplot2-fundamentals - Underlying grammar
- single-cell/markers-annotation - Stacked / split violin for scRNA
- clinical-biostatistics/effect-measures - Effect size to report alongside distribution
