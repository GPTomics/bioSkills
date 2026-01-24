# Spatial Statistics - Usage Guide

## Overview

This skill covers computing spatial statistics for spatial transcriptomics data using Squidpy, including Moran's I, co-occurrence analysis, and neighborhood enrichment.

## Prerequisites

```bash
pip install squidpy scanpy pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Find spatially variable genes"
- "Calculate Moran's I for my spatial data"

## Example Prompts

### Spatial Autocorrelation
> "Compute Moran's I for all genes"

> "Find genes with significant spatial patterns"

### Co-localization
> "Run co-occurrence analysis on my clusters"

> "Which cell types co-localize?"

### Neighborhood Analysis
> "Compute neighborhood enrichment"

> "Are these two clusters spatially associated?"

## What the Agent Will Do

1. Ensure spatial neighbor graph exists
2. Compute requested spatial statistic
3. Return results with significance values
4. Optionally visualize results

## Tips

- **Moran's I** - Positive values indicate spatial clustering
- **Build neighbors first** - All statistics require spatial neighbor graph
- **Gene subset** - Computing for all genes can be slow; subset to HVGs
- **Multiple testing** - Apply FDR correction for many genes
