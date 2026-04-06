# Expression Matrix Normalization - Usage Guide

## Overview

Normalize and transform RNA-seq count matrices using appropriate methods for each downstream task. Covers between-sample normalization (TMM, RLE), within-sample normalization (TPM, FPKM), variance-stabilizing transformations (VST, rlog, log-CPM), GC/length bias correction, and single-cell normalization (scran).

## Prerequisites

```bash
pip install pandas numpy scipy scanpy anndata
```

```r
BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'scran', 'EDASeq', 'cqn'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Normalize my count matrix for PCA and heatmap visualization"
- "Apply TMM normalization to my edgeR DGEList"
- "Transform my counts using VST for clustering"
- "Should I use FPKM or TPM for my cross-sample comparison?"

## Example Prompts
### Choosing a Method
> "I have raw counts from featureCounts and want to do PCA followed by DESeq2 analysis. What normalization should I use for each step?"

> "My samples have very different library sizes (5M to 50M reads). Which transformation handles this best?"

> "I need to provide normalized expression values for a collaborator to use in WGCNA. What format should I export?"

### Applying Normalization
> "Apply DESeq2 VST to my count matrix with blind=FALSE since I know the experimental design"

> "Compute TMM normalization factors and log-CPM values for my edgeR analysis"

> "Normalize my single-cell counts using scran deconvolution with pre-clustering"

### Troubleshooting
> "My PCA shows strong batch effects even after normalization. Should I use ComBat or include batch in the model?"

> "DESeq2 size factors for two samples are very different from the rest. What does this mean?"

> "My GSEA results seem driven by gene length rather than biology. How do I correct for this?"

## What the Agent Will Do

1. Identify the downstream task (DE, visualization, clustering, ML, reporting)
2. Select the appropriate normalization method based on the task and data type
3. Check whether the data is bulk or single-cell, and library size variability
4. Apply normalization using the correct tool and parameters
5. Validate results (check size factor distribution, PCA, mean-variance relationship)
6. Warn about common mistakes (pre-normalizing before DE, using FPKM for cross-sample comparison)

## Method Selection Flowchart

```
What is the downstream task?
|
+-- Differential expression (DESeq2/edgeR/limma) --> Use RAW COUNTS (tools normalize internally)
|
+-- Visualization (PCA, heatmaps, clustering)
|   |
|   +-- Using DESeq2 ecosystem --> VST (blind=FALSE for downstream; blind=TRUE for initial QC)
|   +-- Using edgeR ecosystem --> log-CPM with prior.count >= 2
|   +-- Library sizes vary >10-fold --> rlog (slower but handles unequal libraries)
|   +-- >100 samples --> VST (rlog too slow)
|
+-- WGCNA / gene network analysis --> VST (blind=FALSE)
|
+-- Machine learning / biomarkers --> VST (blind=FALSE)
|
+-- Gene set variation (GSVA/ssGSEA) --> log2(TPM+1) or VST
|
+-- Within-sample gene comparison --> TPM
|
+-- Reporting expression levels --> TPM or DESeq2 normalized counts
|
+-- Single-cell --> scran deconvolution (rigorous) or scanpy normalize_total (exploratory)
```

## Complete R Workflow

```r
library(DESeq2)
library(edgeR)

counts <- read.delim('counts.tsv', row.names=1)
coldata <- read.csv('metadata.csv', row.names=1)
coldata$condition <- factor(coldata$condition, levels=c('control', 'treated'))

# DESeq2 path: normalization + transformation
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- estimateSizeFactors(dds)
cat('Size factors:', sizeFactors(dds), '\n')

# VST for visualization (blind=FALSE for known design)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup='condition')

# For DE: let DESeq() handle normalization from raw counts
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)

# edgeR path: TMM + log-CPM
y <- DGEList(counts=counts, group=coldata$condition)
y <- calcNormFactors(y, method='TMM')
log_cpm <- cpm(y, log=TRUE, prior.count=2)
```

## Complete Python Workflow

```python
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)

# Simple log-CPM for exploratory visualization
lib_sizes = counts.sum(axis=0)
prior_count = 2
cpm_vals = (counts + prior_count) / (lib_sizes + 2 * prior_count) * 1e6
log_cpm = np.log2(cpm_vals)

# PCA on log-CPM
pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_cpm.T)

plt.figure(figsize=(8, 6))
for condition in metadata['condition'].unique():
    mask = metadata['condition'] == condition
    plt.scatter(pca_result[mask, 0], pca_result[mask, 1], label=condition, s=60)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
plt.legend()
plt.tight_layout()
plt.savefig('pca_logcpm.png', dpi=150)

# For rigorous VST, use PyDESeq2 or call R via rpy2
from pydeseq2.dds import DeseqDataSet
dds = DeseqDataSet(counts=counts.T, metadata=metadata, design='~condition')
dds.vst()
vst_counts = pd.DataFrame(dds.layers['vst_counts'], index=dds.obs_names, columns=dds.var_names)
```

## Validating Normalization

```r
# Check size factors are reasonable (expect 0.5 to 2.0 for most samples)
sf <- sizeFactors(dds)
cat('Size factor range:', range(sf), '\n')
if (any(sf < 0.1 | sf > 10)) {
    warning('Extreme size factors detected -- check for outlier samples or contamination')
}

# Check mean-variance relationship after VST (should be flat)
meanSdPlot(assay(vsd))
```

```python
# Check library sizes before and after normalization
print('Raw library sizes:')
print(counts.sum().describe())
print('\nNormalized library sizes (should be more uniform):')
norm_counts = counts.div(size_factors, axis=1)
print(norm_counts.sum().describe())
```

## Tips

- Never provide pre-normalized data to DESeq2, edgeR, or limma-voom -- these tools model raw count distributions
- Use VST with `blind=FALSE` for all downstream analysis except initial unbiased QC
- For edgeR visualization, set `prior.count >= 2` in `cpm(log=TRUE)` to reduce low-count noise
- DESeq2 size factors assume most genes are NOT differentially expressed -- use `controlGenes` or spike-ins when this assumption is violated
- If size factors vary >5-fold, investigate: this often indicates contamination, degradation, or sample mislabeling
- FPKM/RPKM should never be used for cross-sample differential expression -- use TPM for within-sample comparison only
- After log-transformation, check that sparsity decreased (zeros become non-zero); sparse matrix format may no longer save memory
- For gene set enrichment, consider GC/length bias correction (EDASeq, cqn) to avoid false positives driven by gene properties rather than biology

## Related Skills

- expression-matrix/counts-ingest - Load count data before normalization
- expression-matrix/sparse-handling - Sparse format considerations post-normalization
- differential-expression/deseq2-basics - DESeq2 normalization in context of DE
- differential-expression/edger-basics - edgeR TMM normalization
- rna-quantification/count-matrix-qc - QC before normalization
- single-cell/preprocessing - Single-cell normalization workflows
