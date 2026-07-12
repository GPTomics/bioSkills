# Sashimi Plots - Usage Guide

## Overview
Visualize splicing events as sashimi-style plots showing per-sample read coverage and splice junction arcs labeled by read counts. Tools differ in input handling, group aggregation, and customization: ggsashimi (publication overlays), rmats2sashimiplot (rMATS-aware), MAJIQ-VOILA (LSV posteriors interactive HTML), leafviz (leafcutter Shiny), Jutils (tool-agnostic), pyGenomeTracks (multi-track figures).

## Prerequisites
```bash
# Python / CLI
pip install pandas
conda install -c bioconda ggsashimi rmats2sashimiplot pygenometracks

# R for leafviz (leafcutter is on GitHub, not Bioconductor)
# devtools::install_github('davidaknowles/leafcutter/leafcutter')

# MAJIQ-VOILA bundled with MAJIQ V3 (majiq.biociphers.org)

# Jutils
conda install -c bioconda jutils
```

## Quick Start
Tell your AI agent what you want to do:
- "Create sashimi plots for the top 20 differential splicing events from rMATS"
- "Visualize a specific exon-skipping event with samples grouped by condition"
- "Generate publication-quality sashimi with intron shrinking and per-condition aggregation"
- "Make MAJIQ VOILA interactive HTML for browsing LSV posteriors"
- "Plot splicing alongside ChIP-seq tracks for the same locus"

## Example Prompts

### Single Event Visualization
> "Plot a sashimi for chr17:43094000-43125000 (BRCA1 region) with control vs treatment groups, intron shrinking, and matched y-axis scales."

### Batch Plotting from Differential Output
> "Iterate ggsashimi over the top 25 significant rMATS events; output per-event PDFs with flanking 500nt context."

### MAJIQ VOILA
> "Run voila on MAJIQ deltapsi output to generate interactive HTML browser of LSV posteriors and splice graphs."

### leafviz
> "Generate leafcutter Shiny app for interactive cluster-level browsing with NMD annotation."

### Multi-Track Figures
> "Use pyGenomeTracks to combine RNA-seq coverage tracks with H3K4me3 ChIP-seq for the same locus."

### Tool-Agnostic Heatmaps
> "Use Jutils to create a unified heatmap of significant events across rMATS, leafcutter, and MAJIQ output."

## What the Agent Will Do
1. Choose visualization tool based on the upstream differential analysis (rMATS, leafcutter, MAJIQ) or general region request
2. Configure sample grouping with consistent colors per condition
3. Apply intron shrinking, fixed y-scales, and replicate aggregation flags
4. Iterate over significant events for batch plotting
5. Generate publication-ready PDFs/SVGs with annotation tracks

## Tips
- Use `--shrink` for genes with large introns (TTN, brain genes) to keep exons visible
- `--fix-y-scale` is essential for cross-group comparisons; otherwise auto-rescaling exaggerates differences
- Aggregate replicates with `-O 3 -A mean_j` to reduce clutter while preserving variance via alpha
- Limit to 3-4 groups per figure; more becomes hard to read
- Include 200-500 nt flanking exons for splicing context
- For MXE events, plot both alternative exons; otherwise only half of the event is visible
- VOILA requires MAJIQ build output (splicegraph.zarr in V3; the V2 splicegraph.sql is deprecated) plus the .voila file
- IGV sashimi is for interactive ad-hoc inspection, not figures

## Related Skills

- differential-splicing - Identify events to plot
- splicing-quantification - Context for PSI values
- data-visualization/genome-tracks - Multi-track figure design
- data-visualization/ggplot2-fundamentals - ggsashimi customization
