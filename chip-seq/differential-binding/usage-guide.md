# Differential Binding Analysis - Usage Guide

## Overview

Identify differentially bound regions between ChIP-seq conditions. Supports DiffBind (from BAM/peak files), DESeq2 (R from count matrix), and PyDESeq2 (Python from count matrix). All approaches apply normalization and statistical testing to detect significant binding changes.

## Prerequisites

```r
BiocManager::install(c('DiffBind', 'DESeq2', 'apeglm'))
```

```bash
pip install pydeseq2 pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Find peaks that change between treatment and control"
- "Identify differential binding sites with FDR < 0.05"
- "Compare H3K27ac peaks between wild-type and knockout"
- "Run differential binding on my count matrix"
- "Find differentially bound regions from this peaks-by-samples count table"

## Example Prompts

### From Count Matrix
> "I have a ChIP-seq count matrix with treated and control columns. Find differentially bound peaks."

> "Run DESeq2 on my peak count matrix to identify differential binding between conditions"

> "Use Python to find differentially bound regions from counts.tsv"

### From BAM Files (DiffBind)
> "Run differential binding analysis comparing treated vs untreated samples using DiffBind"

> "Create a DiffBind sample sheet from my BAM and peak files"

> "Find regions with increased binding in the drug treatment condition"

### Normalization
> "My experiment involves global H3K27me3 loss — what normalization should I use for differential binding?"

> "Compare binding changes using full library size normalization instead of default"

### Results Interpretation
> "Generate an MA plot of differential binding results"

> "Export the significant differential peaks to a BED file"

## What the Agent Will Do

### From Count Matrix
1. Load the count matrix and construct sample metadata from column names
2. Create a DESeqDataSet (R) or DeseqDataSet (Python) with appropriate design formula
3. Apply minimal pre-filtering appropriate for ChIP-seq (less aggressive than RNA-seq)
4. Run the DESeq2 pipeline (size factor estimation, dispersion estimation, statistical testing)
5. Extract results at the specified significance threshold (default padj < 0.05)
6. Export results with peak IDs, log2 fold changes, p-values, and significance calls

### From BAM Files (DiffBind)
1. Create a sample sheet with BAM files, peak files, and condition labels
2. Load samples and count reads in consensus peak regions
3. Normalize counts and set up the experimental contrast
4. Run differential binding analysis using DESeq2 or edgeR
5. Generate diagnostic plots (PCA, MA plot, volcano plot)
6. Export significant differential peaks with fold changes and FDR values

## Tips

- Requires at least 2 replicates per condition (3+ recommended for adequate statistical power)
- Positive fold change means higher binding in treatment; negative means lower
- Check PCA plot to verify samples cluster by condition, not batch
- Normalization matters most: DESeq2 default RLE normalization on a count matrix assumes most peaks are stable between conditions. This works for typical TF ChIP but fails when global binding changes occur (e.g., EZH2 inhibitor, BET inhibitor). Use full library size normalization or spike-in normalization in those cases.
- Pre-filter conservatively for ChIP-seq: peaks are already enriched regions. RNA-seq-style filtering (e.g., rowSums >= 10) can remove truly differential peaks that have condition-specific signal.
- Set `alpha` in `results()` to match the intended significance threshold (e.g., `alpha = 0.05` for padj < 0.05) for optimal independent filtering
- LFC shrinkage (apeglm) improves fold-change estimates for ranking but does not change p-values or padj
- If few peaks are significant, check replicate concordance, peak quality, and normalization choice

## Related Skills

- peak-calling - Generate input peak files for DiffBind
- peak-annotation - Annotate differential peaks to genomic features
- differential-expression/deseq2-basics - DESeq2 fundamentals
- chipseq-visualization - Heatmaps and profile plots of differential regions
