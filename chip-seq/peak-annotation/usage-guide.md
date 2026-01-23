# Peak Annotation - Usage Guide

## Overview

Annotate ChIP-seq peaks to genomic features (promoters, exons, introns, intergenic) and associated genes using ChIPseeker.

## Prerequisites

```r
BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))

# For functional enrichment
BiocManager::install('clusterProfiler')
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my peaks to the nearest genes"
- "Show the distribution of peaks across genomic features"
- "Find which genes have peaks in their promoters"

## Example Prompts

### Basic Annotation
> "Annotate my narrowPeak file to genes and genomic features"

> "Get gene symbols for all peaks in my BED file"

> "Find the nearest TSS for each peak"

### Visualization
> "Create a pie chart showing peak distribution across genomic features"

> "Plot the distance of peaks to the nearest TSS"

> "Compare peak annotation between two samples"

### Functional Analysis
> "Run GO enrichment on genes with promoter peaks"

> "Export annotated peaks with gene symbols to CSV"

## What the Agent Will Do

1. Load peak files (narrowPeak, broadPeak, or BED format)
2. Match chromosome naming style (UCSC vs Ensembl)
3. Annotate peaks to genomic features using the appropriate TxDb
4. Add gene symbols and Entrez IDs from the annotation database
5. Generate visualization plots (pie chart, bar plot, TSS distance plot)
6. Export annotated peaks with all gene information

## Tips

- Ensure chromosome names match between peaks and TxDb (use `seqlevelsStyle()`)
- Include `annoDb` parameter to get gene symbols, not just Entrez IDs
- Promoter regions are typically defined as +/- 3kb from TSS
- Use `distanceToTSS` to filter for proximal vs distal peaks
- Multiple peaks can annotate to the same gene
