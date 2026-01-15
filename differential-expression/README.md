# Differential Expression

Differential expression analysis using R/Bioconductor packages DESeq2 and edgeR for RNA-seq count data.

## Overview

This category covers the complete differential expression workflow: analyzing RNA-seq count data with DESeq2 or edgeR, visualizing results with MA plots, volcano plots, PCA, and heatmaps, and extracting significant genes for downstream analysis.

**Tool type:** `r`
**Primary tools:** DESeq2, edgeR, ggplot2, pheatmap

## Skills

| Skill | Description |
|-------|-------------|
| [deseq2-basics](deseq2-basics/) | DESeq2 workflow: DESeqDataSet, normalization, testing, lfcShrink |
| [edger-basics](edger-basics/) | edgeR workflow: DGEList, TMM normalization, glmQLFit, testing |
| [de-visualization](de-visualization/) | MA plots, volcano plots, PCA, heatmaps with ggplot2/pheatmap |
| [de-results](de-results/) | Filter significant genes, add annotations, export results |

## Workflow

```
Count Matrix + Sample Metadata
    |
    v
[deseq2-basics] or [edger-basics]
    |
    +---> Normalization
    |         |
    |         v
    |     Dispersion estimation
    |         |
    |         v
    |     Statistical testing
    |
    v
DE Results (log2FC, p-values)
    |
    +---> [de-visualization] - QC plots, publication figures
    |
    +---> [de-results] - Filter, annotate, export
              |
              v
          Gene lists for pathway analysis
```

## DESeq2 vs edgeR

| Aspect | DESeq2 | edgeR |
|--------|--------|-------|
| Model | Negative binomial | Negative binomial + QL |
| Normalization | Median of ratios | TMM |
| Shrinkage | LFC shrinkage (apeglm) | Empirical Bayes on dispersions |
| Best for | General use | Complex designs, small samples |
| Output | Wald statistic, padj | F-statistic, FDR |

## Example Prompts

### DESeq2
- "Run DESeq2 on my count matrix with treated vs control"
- "Analyze differential expression controlling for batch"
- "Apply log fold change shrinkage with apeglm"

### edgeR
- "Run edgeR quasi-likelihood analysis on my RNA-seq data"
- "Create contrasts to compare all treatments against control"
- "Use TMM normalization and glmQLFit"

### Visualization
- "Create a volcano plot highlighting significant genes"
- "Make a heatmap of top 50 DE genes"
- "Generate PCA plot colored by treatment group"

### Results
- "Extract genes with padj < 0.05 and |log2FC| > 1"
- "Add gene symbols to my results"
- "Export significant genes to Excel"

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

# Core packages
BiocManager::install(c('DESeq2', 'edgeR', 'apeglm'))

# Visualization
install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer', 'ggrepel'))
BiocManager::install('EnhancedVolcano')

# Annotation
BiocManager::install(c('org.Hs.eg.db', 'biomaRt'))

# Export
install.packages('openxlsx')
```

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Count matrix | Integer matrix | Genes (rows) x Samples (columns) |
| Sample metadata | Data frame | Condition, batch, etc. |
| Design formula | R formula | ~condition or ~batch + condition |

## Notes

- **Biological replicates required** - DESeq2 removed support for no-replicate designs in v1.22
- **Use lfcShrink()** - DESeq2's betaPrior is deprecated; use lfcShrink() with type='apeglm'
- **glmQLFit preferred** - edgeR's quasi-likelihood framework is more robust than exact test
- **decidetestsDGE() removed** - Use decideTests() instead in edgeR v4.4+
- **estimateDisp() optional** - In edgeR v4+, glmQLFit() estimates dispersions internally

## Related Skills

- **sequence-io** - Prepare input sequences
- **alignment-files** - Process BAM files for counting
- **database-access** - Fetch gene annotations from NCBI
- **pathway-analysis** (planned) - GO/KEGG enrichment of DE genes

## References

- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
- [RNA-seq workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
