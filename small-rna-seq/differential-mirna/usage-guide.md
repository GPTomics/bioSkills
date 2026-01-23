# Differential miRNA Expression - Usage Guide

## Overview

Identify differentially expressed miRNAs between conditions using DESeq2 or edgeR with considerations specific to small RNA data.

## Prerequisites

```r
BiocManager::install(c('DESeq2', 'edgeR', 'apeglm', 'EnhancedVolcano', 'pheatmap'))
```

## Quick Start

Tell your AI agent:
- "Find differentially expressed miRNAs between treatment and control"
- "Run DESeq2 on my miRNA count matrix"
- "Create a volcano plot of DE miRNAs"
- "Filter for significant miRNAs with |log2FC| > 1"

## Example Prompts

### Basic DE Analysis

> "Run differential expression on my miRge3 counts"

> "Compare miRNA expression between tumor and normal samples"

> "Find miRNAs with adjusted p-value < 0.05"

### Visualization

> "Create a volcano plot of differentially expressed miRNAs"

> "Make a heatmap of significant miRNAs"

> "Plot MA of miRNA expression changes"

### Export and Filtering

> "Export significant miRNAs to CSV"

> "Filter for miRNAs with at least 2-fold change"

> "List the top 20 upregulated miRNAs"

## What the Agent Will Do

1. Load miRNA count matrix and sample metadata
2. Create DESeq2 dataset with appropriate design
3. Filter low-expressed miRNAs (< 10 total reads)
4. Run DESeq2 and apply apeglm shrinkage
5. Filter by significance (padj < 0.05, |log2FC| > 1)
6. Generate visualizations and export results

## Tips

- **Filter low counts** - miRNAs with < 10 total reads are unreliable
- **Use apeglm shrinkage** - improves log2FC estimates for low-count miRNAs
- **Check normalization** - miRNA libraries can have different compositions
- **Multiple testing** - always use adjusted p-values
- **|log2FC| > 1** is standard for biologically meaningful changes
