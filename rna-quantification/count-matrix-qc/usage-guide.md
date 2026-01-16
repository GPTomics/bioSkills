# Count Matrix QC Usage Guide

Quality control of count matrices is essential before differential expression analysis to identify outliers, batch effects, and sample quality issues.

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
