# Repertoire Visualization - Usage Guide

## Overview

Renders publication figures for TCR/BCR repertoires and encodes how to read them so they are not misinterpreted. Every figure inherits two upstream choices: the clonotype definition (nucleotide vs amino-acid CDR3, with or without V/J) and the sequencing depth. Diversity, rarefaction, and overlap figures are only comparable across samples after depth normalization (downsample to a common depth, or read rarefaction curves at a shared x); a bare Shannon bar across unequal-depth libraries is misleading. Amino-acid CDR3 inflates apparent sharing via convergent recombination, so the clonotype definition must be held constant across every figure in a comparison and stated on it. The skill covers V-J chord/circos, CDR3 spectratype, clonal-space stratification, clonal tracking, rarefaction/extrapolation, overlap heatmaps, and clonotype-similarity networks.

## Prerequisites

```bash
pip install matplotlib seaborn pandas numpy networkx rapidfuzz
# R plotting routes
Rscript -e "install.packages(c('circlize', 'iNEXT'))"
# VDJtools (Java) provides PlotFancyVJUsage and RarefactionPlot; run RInstall once for its R deps
```

## Quick Start

Tell your AI agent what to plot:
- "Plot a V-J pairing chord diagram for this sample"
- "Draw a CDR3 spectratype and tell me if it looks polyclonal"
- "Compare diversity with rarefaction curves, not a raw Shannon bar"
- "Make an overlap heatmap using a depth-robust metric"
- "Track my top 10 clones across timepoints"
- "Build a clonotype-similarity network and show the threshold effect"

## Example Prompts

### Gene Usage

> "Create a V-J pairing chord diagram for sample S1, and a V-by-J frequency heatmap as an alternative"

> "Show V gene usage across samples, and warn me if the differences could be primer bias"

### Clonal Structure

> "Plot a frequency-weighted CDR3 spectratype and interpret whether it is Gaussian (polyclonal) or spiked (expanded)"

> "Show clonal-space homeostasis as a stacked bar of rare/small/medium/large/hyperexpanded strata"

### Diversity Done Right

> "Compare diversity between treatment and control with rarefaction/extrapolation curves read at a common depth"

> "I have samples at very different depths - make a diversity figure that is not confounded by depth"

### Comparison and Dynamics

> "Generate an overlap heatmap with Morisita-Horn on depth-normalized samples"

> "Track the top clones over the vaccination timecourse and flag sampling-zero contractions"

> "Build a CDR3-similarity network and run a threshold sensitivity sweep"

## What the Agent Will Do

1. Confirm the clonotype definition (nt vs aa, +/- V/J) and hold it constant across figures.
2. Check whether depth normalization is needed and, for diversity/overlap, downsample or use rarefaction.
3. Compute the metric behind the figure (frequency strata, Hill numbers, overlap matrix, pairwise distances).
4. Render the chosen figure with matplotlib/seaborn, R circlize/iNEXT, or VDJtools plotting.
5. Annotate the figure with the metric, weighting, and threshold so it reads unambiguously.

## Tips

- Rarefaction curves, not a Shannon bar, are the correct diversity comparison across unequal depth; read all curves at a shared x and extrapolate at most 2-3x observed depth.
- For overlap heatmaps prefer Morisita-Horn (abundance-weighted, depth-robust) over Jaccard (presence/absence, depth-biased), and normalize depth first.
- A similarity network's structure depends entirely on the distance threshold - report the metric and cutoff, and sweep it.
- Spectratype: weight by frequency to see expansions, by unique clonotypes to see diversity; label which is plotted.
- Clonal tracking: downsample timepoints to common depth before calling a clone "contracted"; absence is often a sampling zero.
- Never compare V/J usage across different primer sets or platforms - primer-amplification bias mimics biology.
- Use colorblind-friendly palettes (viridis, Set2) and keep clone-tracking to the top 10-20 clones for legibility.

## Related Skills

- vdjtools-analysis - Compute the diversity/overlap inputs (depth-normalized)
- mixcr-analysis - Produce clonotype tables
- scirpy-analysis - Single-cell clonal-expansion overlays
- data-visualization/ggplot2-fundamentals - General ggplot2 grammar
- data-visualization/heatmaps-clustering - Overlap/usage heatmap techniques
