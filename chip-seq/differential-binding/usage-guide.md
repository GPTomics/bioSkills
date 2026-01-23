# Differential Binding Analysis - Usage Guide

## Overview

Identify differentially bound regions between ChIP-seq conditions using DiffBind, which applies normalization and statistical testing to read counts in consensus peaks.

## Prerequisites

```r
BiocManager::install('DiffBind')
```

## Quick Start

Tell your AI agent what you want to do:
- "Find peaks that change between treatment and control"
- "Identify differential binding sites with FDR < 0.05"
- "Compare H3K27ac peaks between wild-type and knockout"

## Example Prompts

### Basic Analysis
> "Run differential binding analysis comparing treated vs untreated samples"

> "Find regions with increased binding in the drug treatment condition"

> "Identify peaks lost upon knockdown of my transcription factor"

### Sample Setup
> "Create a DiffBind sample sheet from my BAM and peak files"

> "Set up a contrast comparing tumor vs normal samples"

### Results Interpretation
> "Generate an MA plot of differential binding results"

> "Export the significant differential peaks to a BED file"

## What the Agent Will Do

1. Create a sample sheet with BAM files, peak files, and condition labels
2. Load samples and count reads in consensus peak regions
3. Normalize counts and set up the experimental contrast
4. Run differential binding analysis using DESeq2 or edgeR
5. Generate diagnostic plots (PCA, MA plot, volcano plot)
6. Export significant differential peaks with fold changes and FDR values

## Tips

- Requires at least 2 replicates per condition (3+ recommended)
- Use `summits = 250` to re-center peaks on summit and standardize width
- Positive fold change means gained in treatment, negative means lost
- Check PCA plot to verify samples cluster by condition, not batch
- If few peaks are significant, check replicate concordance and peak quality
