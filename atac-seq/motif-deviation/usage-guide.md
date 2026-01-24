# Motif Deviation Analysis

## Overview

Analyze transcription factor motif accessibility variability across samples using chromVAR. This tool identifies TFs whose binding sites show differential accessibility between conditions, complementing peak-level differential analysis.

## Prerequisites

```r
BiocManager::install(c('chromVAR', 'motifmatchr', 'JASPAR2020', 'TFBSTools'))
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')  # or your genome
```

## Quick Start

Tell your AI agent what you want to do:
- "Run chromVAR on my ATAC-seq peaks and counts"
- "Find which TF motifs have variable accessibility across my samples"
- "Identify differential motif accessibility between treatment and control"
- "Generate a heatmap of top variable TF motifs"

## Example Prompts

### Basic Analysis
> "Run chromVAR on my ATAC-seq data to find variable TF motifs"

> "Compute motif deviation scores for my samples"

> "Which transcription factors show the highest variability in accessibility?"

### Differential Analysis
> "Compare TF motif accessibility between my two conditions"

> "Find motifs that are more accessible in treated vs control samples"

> "Run differential motif analysis with limma on chromVAR deviations"

### Visualization
> "Create a heatmap of top 50 variable TF motifs"

> "Plot PCA of chromVAR deviation scores colored by condition"

> "Make a volcano plot of differential motif accessibility"

### Integration
> "Combine chromVAR results with my DiffBind differential peaks"

> "Correlate TF motif deviations with gene expression"

## What the Agent Will Do

1. Load peak counts and genomic ranges
2. Add GC bias correction for accurate deviation estimation
3. Filter low-quality peaks and samples
4. Match TF motifs from JASPAR to peaks
5. Compute per-sample deviation scores
6. Calculate variability across samples
7. Identify differentially accessible motifs between conditions
8. Generate visualizations (heatmaps, PCA, volcano plots)

## Tips

- chromVAR works best with 10,000+ peaks and >1500 reads per sample
- GC bias correction is essential for accurate deviation estimates
- Use JASPAR vertebrate core collection for comprehensive TF coverage
- Variability > 2.0 indicates highly variable motifs worth investigating
- Combine with footprinting (TOBIAS) for TF binding confirmation
- For scATAC-seq, aggregate cells by cluster before running chromVAR
- Background peak selection affects results; use getBackgroundPeaks for custom control
