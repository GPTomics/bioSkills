# differential-expression

## Overview

Differential expression analysis using R/Bioconductor packages DESeq2 and edgeR for RNA-seq count data. Covers the complete workflow from count matrix to visualizations and significant gene lists, including decision guidance for method selection, result interpretation, and prokaryotic organisms.

**Tool type:** r | **Primary tools:** DESeq2, edgeR, ggplot2, pheatmap

## Skills

| Skill | Description |
|-------|-------------|
| deseq2-basics | DESeq2 workflow with decision guidance, prokaryotic support, PyDESeq2 alternative |
| edger-basics | edgeR quasi-likelihood workflow with test selection and v4 guidance |
| batch-correction | Batch correction with ComBat, limma, SVA; when to use each method |
| de-visualization | MA plots, volcano plots, PCA, heatmaps with diagnostic interpretation |
| de-results | Filter, annotate, export results; interpretation guidance, GSEA/ORA preparation |
| timeseries-de | Time-series DE with limma splines, maSigPro, ImpulseDE2 |

## Example Prompts

- "Run DESeq2 on my count matrix with treated vs control"
- "Analyze differential expression controlling for batch"
- "Apply log fold change shrinkage with apeglm"
- "Run edgeR quasi-likelihood analysis on my RNA-seq data"
- "Create contrasts to compare all treatments against control"
- "Use TMM normalization and glmQLFit"
- "Create a volcano plot highlighting significant genes"
- "Make a heatmap of top 50 DE genes"
- "Generate PCA plot colored by treatment group"
- "Check p-value histogram for analysis quality"
- "Extract genes with padj < 0.05 and |log2FC| > 1"
- "Add gene symbols to my results"
- "Prepare a ranked gene list for GSEA from my DE results"
- "Run DESeq2 on my bacterial RNA-seq data with appropriate KEGG annotation"
- "Run DE analysis using PyDESeq2 in Python"

## Requirements

```r
BiocManager::install(c('DESeq2', 'edgeR', 'apeglm'))
install.packages(c('ggplot2', 'pheatmap', 'ggrepel'))
```

```bash
# Python alternative
pip install pydeseq2
```

## Related Skills

- **rna-quantification** - Generate count matrices from BAM or FASTQ
- **pathway-analysis** - GO/KEGG enrichment of DE genes
- **single-cell** - Single-cell differential expression
- **alignment-files** - Process BAM files for counting
- **expression-matrix** - Count matrix loading and gene ID mapping
