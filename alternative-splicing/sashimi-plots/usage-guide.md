# Sashimi Plots - Usage Guide

## Overview
Create sashimi plots to visualize splicing events with read coverage and splice junction counts. These plots show coverage across exons and arcs representing junction-spanning reads, essential for validating differential splicing results.

## Prerequisites
```bash
# ggsashimi (recommended)
pip install ggsashimi

# rmats2sashimiplot (for rMATS output)
pip install rmats2sashimiplot

# Dependencies
pip install matplotlib pysam pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Create sashimi plots for my significant splicing events"
- "Visualize the exon skipping event at gene X"
- "Generate sashimi plots comparing control vs treatment"
- "Plot junction coverage for my top differential splicing results"

## Example Prompts

### Single Event
> "Create a sashimi plot for the region chr1:1000000-1010000 showing my control and treatment samples."

> "Visualize the exon skipping event in TP53 with junction read counts."

### Batch Processing
> "Generate sashimi plots for my top 20 significant differential splicing events from rMATS."

> "Create publication-ready sashimi plots with shrunk introns for all significant events."

### Customization
> "Make sashimi plots with samples grouped by condition and color-coded."

> "Generate high-resolution PDF sashimi plots with consistent y-axis scaling."

## What the Agent Will Do
1. Create sample grouping file with BAM paths and conditions
2. Define genomic regions for events of interest
3. Generate sashimi plots with coverage and junction arcs
4. Apply visualization options (shrunk introns, fixed y-scale)
5. Export plots in requested format (PDF, PNG, SVG)

## Tips
- Use `--shrink` option for genes with large introns
- Set `--fix-y-scale` when comparing between groups
- Aggregate replicates with `-A mean` to reduce clutter
- Minimum junction reads (-M) of 5-10 filters noise
- Limit to 3-4 groups per plot for readability
- Include flanking exons to show full splicing context

## Related Skills
- differential-splicing - Identify events to plot
- splicing-quantification - Context for PSI values
- data-visualization/ggplot2-fundamentals - Further customization
