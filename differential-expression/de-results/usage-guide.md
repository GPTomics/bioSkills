# DE Results - Usage Guide

## Overview

This skill covers extracting, filtering, annotating, and exporting differential expression results for downstream analysis. Works with both DESeq2 and edgeR output.

## Prerequisites

```r
# For annotation
BiocManager::install(c('org.Hs.eg.db', 'biomaRt'))
# For Excel export
install.packages('openxlsx')
```

## Quick Start

Tell your AI agent what you want to do:
- "Extract significant genes with padj < 0.05 and |log2FC| > 1"
- "Add gene symbols to my DESeq2 results"
- "Export results to Excel with separate sheets for up and down genes"

## Example Prompts

### Filtering
> "Get the top 100 most significant genes"

> "Filter for up-regulated genes with at least 2-fold change"

> "Find genes significant in both DESeq2 and edgeR"

### Annotation
> "Add gene symbols and descriptions from Ensembl"

> "Convert Ensembl IDs to gene symbols"

> "Annotate results with my custom annotation file"

### Export
> "Save significant genes to CSV"

> "Create an Excel file with all results and significant genes"

> "Prepare a ranked gene list for GSEA"

## What the Agent Will Do

1. Extract results from DESeq2/edgeR object
2. Apply requested filters (padj, log2FC, baseMean)
3. Add gene annotations if requested
4. Format and export results

## Common Thresholds

| Filter | Typical Value | Description |
|--------|--------------|-------------|
| padj | < 0.05 | Standard significance |
| padj | < 0.01 | Stringent |
| \|log2FC\| | > 1 | 2-fold change |
| \|log2FC\| | > 0.585 | 1.5-fold change |
| baseMean | > 10 | Minimum expression |

## Tips

- Always use adjusted p-values (padj/FDR), not raw p-values
- Check for NA values in padj (indicates filtering or outliers)
- Use lfcShrink() results for ranking genes
- Export both all results and filtered significant genes
- Include multiple test corrections in reports
