# differential-expression

## Overview

Differential expression analysis using R/Bioconductor packages DESeq2 and edgeR for RNA-seq count data. Covers the complete workflow from count matrix to visualizations and significant gene lists.

**Tool type:** r | **Primary tools:** DESeq2, edgeR, ggplot2, pheatmap

## Skills

| Skill | Description |
|-------|-------------|
| deseq2-basics | DESeq2 workflow: DESeqDataSet, normalization, testing, lfcShrink |
| edger-basics | edgeR workflow: DGEList, TMM normalization, glmQLFit, testing |
| batch-correction | Batch correction with ComBat, limma, SVA |
| de-visualization | MA plots, volcano plots, PCA, heatmaps with ggplot2/pheatmap |
| de-results | Filter significant genes, add annotations, export results |
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
- "Create an MA plot showing DE genes"
- "Extract genes with padj < 0.05 and |log2FC| > 1"
- "Add gene symbols to my results"
- "Export significant genes to Excel"
- "How many genes are up vs down regulated?"

## Requirements

```r
BiocManager::install(c('DESeq2', 'edgeR', 'apeglm'))
install.packages(c('ggplot2', 'pheatmap', 'ggrepel'))
```

## Related Skills

- **rna-quantification** - Generate count matrices from BAM or FASTQ
- **pathway-analysis** - GO/KEGG enrichment of DE genes
- **single-cell** - Single-cell differential expression
- **alignment-files** - Process BAM files for counting
