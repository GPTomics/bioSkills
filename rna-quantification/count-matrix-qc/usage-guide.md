# Count Matrix QC - Usage Guide

## Overview
Quality control of count matrices is essential before differential expression analysis to identify outliers, batch effects, and sample quality issues.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('DESeq2')
install.packages('pheatmap')
```

## Quick Start
Tell your AI agent what you want to do:
- "Check the quality of my count matrix before DE analysis"
- "Identify outlier samples in my RNA-seq data"
- "Make a PCA plot to see if my samples cluster by condition"

## Example Prompts
### Basic QC
> "Run QC on my count matrix and identify any problematic samples"

> "Check library sizes and gene detection rates across my samples"

### Visualization
> "Create a sample correlation heatmap from my count data"

> "Make a PCA plot colored by condition and batch"

### Outlier Detection
> "Check if any samples are outliers based on PCA and correlation"

> "Should I remove sample X based on the QC metrics?"

### Batch Effects
> "Check if my samples show batch effects in PCA"

> "Add batch correction to my DESeq2 design"

## What the Agent Will Do
1. Load the count matrix and sample metadata
2. Calculate library sizes and gene detection rates
3. Normalize counts for visualization (VST or rlog)
4. Generate sample correlation heatmap
5. Create PCA plot to assess sample clustering
6. Identify any outliers or batch effects
7. Recommend next steps based on QC results

## Key QC Checks

1. **Library sizes** - Total counts per sample
2. **Gene detection** - Number of genes with counts
3. **Sample correlation** - Replicates should cluster together
4. **PCA** - Samples should separate by condition, not batch
5. **Outlier detection** - Identify problematic samples

## Quick QC in R

```r
library(DESeq2)
library(pheatmap)

# Load data
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Normalize
vsd <- vst(dds, blind = TRUE)

# Sample correlation
pheatmap(cor(assay(vsd)))

# PCA
plotPCA(vsd, intgroup = 'condition')
```

## What to Look For

### Good Signs
- Replicates cluster together in PCA
- High correlation (>0.9) between replicates
- Similar library sizes across samples
- Clear separation between conditions

### Warning Signs
- Sample doesn't cluster with its group
- Low correlation with replicates
- Very different library size
- PCA driven by batch, not condition

## Handling Problems

### Outlier Samples
```r
# Option 1: Remove outlier
dds <- dds[, colnames(dds) != 'outlier_sample']

# Option 2: Flag for sensitivity analysis
# Run with and without outlier
```

### Batch Effects
```r
# Add batch to design
design(dds) <- ~ batch + condition
```

### Low Library Size
- Consider excluding samples with <1M reads
- Or use weighted analysis

## Recommended Thresholds

| Metric | Good | Concerning |
|--------|------|------------|
| Library size | >5M | <1M |
| Genes detected | >12,000 | <8,000 |
| Replicate correlation | >0.95 | <0.85 |
| Mapping rate | >70% | <50% |

## Tips
- Always run QC before differential expression analysis
- Use blind=TRUE for vst/rlog when doing unsupervised QC
- Document any samples removed and justify the decision
- Run sensitivity analysis if unsure about removing a sample
- Check for batch effects even if not explicitly part of the design
