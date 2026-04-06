---
name: bio-expression-matrix-normalization
description: Normalize and transform RNA-seq count matrices for differential expression, visualization, and clustering. Covers between-sample (TMM, RLE, upper quartile), within-sample (TPM, FPKM), variance-stabilizing (VST, rlog), and single-cell (scran) methods. Use when choosing or applying normalization to expression data.
tool_type: mixed
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, pandas 2.2+, numpy 1.26+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Expression Matrix Normalization

## Why Normalization Matters: Composition Bias

Simple library-size normalization (dividing by total counts) fails because RNA populations differ between samples. If a treatment massively upregulates a few genes, those genes consume a disproportionate share of sequencing reads, making every other gene appear downregulated even if unchanged. This composition bias is not a theoretical edge case -- it occurs routinely in stress responses, viral infection, and knockout experiments.

All robust normalization methods (TMM, RLE/DESeq2, upper quartile) address composition bias by estimating scaling factors from the majority of non-differentially-expressed genes.

## Normalization Decision Table

| Task | Method | Tool | Input |
|------|--------|------|-------|
| DE analysis (DESeq2) | RLE / median-of-ratios | `DESeq()` handles internally | Raw integer counts |
| DE analysis (edgeR) | TMM | `calcNormFactors()` | Raw integer counts |
| DE analysis (limma-voom) | TMM + voom | `calcNormFactors()` then `voom()` | Raw integer counts |
| PCA, heatmaps, clustering | VST or rlog | `vst(dds)` or `rlog(dds)` | DESeqDataSet |
| PCA/clustering (edgeR/limma) | log-CPM | `edgeR::cpm(y, log=TRUE)` | DGEList |
| WGCNA | VST or log-CPM | `vst(dds, blind=FALSE)` | DESeqDataSet |
| GSVA / ssGSEA | log2(TPM+1) or VST | Precomputed | TPM or DESeqDataSet |
| Machine learning / biomarkers | VST | `vst(dds, blind=FALSE)` | DESeqDataSet |
| Reporting expression levels | TPM (within-sample) | Quantification tool output | N/A |
| Cross-sample comparison | DESeq2 normalized counts | `counts(dds, normalized=TRUE)` | DESeqDataSet |
| Single-cell RNA-seq | scran deconvolution | `scran::computeSumFactors()` | SingleCellExperiment |

Critical rule: DE tools (DESeq2, edgeR, limma-voom) normalize internally from raw counts. Never provide pre-normalized data to these tools.

## Between-Sample Normalization

### RLE / Median of Ratios (DESeq2)

**Goal:** Estimate per-sample size factors that correct for library size and composition bias.

**Approach:** Compute a pseudo-reference sample (geometric mean per gene across all samples), then take the median of per-gene ratios between each sample and the reference. The median is robust to the minority of DE genes.

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Normalized counts for export/visualization
norm_counts <- counts(dds, normalized=TRUE)
```

Size factor interpretation: a factor of 1.2 means the sample has 20% more sequencing depth (after composition adjustment) than the reference. Genes with zero counts in any sample are excluded from the geometric mean calculation.

Size factor alternatives for special cases:

| Scenario | Solution |
|----------|----------|
| High zero-inflation (single-cell) | `estimateSizeFactors(dds, type='poscounts')` |
| Very small libraries | `estimateSizeFactors(dds, type='iterate')` |
| Known stable reference genes | `estimateSizeFactors(dds, controlGenes=stable_idx)` |
| Majority of genes DE (e.g., prokaryotic stress) | Spike-in normalization or `controlGenes` |

### TMM (edgeR)

**Goal:** Compute normalization factors that account for composition bias using trimmed mean of M-values.

**Approach:** Select a reference sample, compute gene-wise log-ratios (M-values) and average expression (A-values) relative to the reference, trim extreme values, and compute a weighted mean scaling factor.

```r
library(edgeR)

y <- DGEList(counts=counts, group=metadata$condition)
y <- calcNormFactors(y, method='TMM')

# TMM trims 30% of M-values and 5% of A-values by default
# Normalization factors are stored in y$samples$norm.factors
# Effective library size = lib.size * norm.factors
```

```python
# Python equivalent using edgepy or manual calculation
import numpy as np
import pandas as pd

def tmm_norm_factors(counts):
    '''Simplified TMM normalization factors.'''
    lib_sizes = counts.sum(axis=0)
    ref_idx = np.argmin(np.abs(np.log(lib_sizes) - np.mean(np.log(lib_sizes))))
    ref = counts.iloc[:, ref_idx]
    factors = np.ones(counts.shape[1])
    for i in range(counts.shape[1]):
        if i == ref_idx:
            continue
        sample = counts.iloc[:, i]
        keep = (ref > 0) & (sample > 0)
        m = np.log2(sample[keep] / lib_sizes[i]) - np.log2(ref[keep] / lib_sizes[ref_idx])
        a = 0.5 * (np.log2(sample[keep] / lib_sizes[i]) + np.log2(ref[keep] / lib_sizes[ref_idx]))
        m_lo, m_hi = np.percentile(m, [15, 85])
        a_lo, a_hi = np.percentile(a, [2.5, 97.5])
        keep2 = (m >= m_lo) & (m <= m_hi) & (a >= a_lo) & (a <= a_hi)
        factors[i] = 2 ** np.average(m[keep2])
    return factors / np.exp(np.mean(np.log(factors)))
```

### Upper Quartile

```r
y <- calcNormFactors(y, method='upperquartile')
```

Simpler than TMM; uses the 75th percentile of non-zero counts. Less robust to asymmetric DE but adequate when composition bias is mild.

## Within-Sample Normalization (TPM, FPKM)

### Why FPKM/RPKM Is Problematic

FPKM normalizes by total mapped reads AND gene length, enabling comparison of Gene A to Gene B within one sample. However, the per-million scaling factor depends on the RNA composition of that sample. If one sample has a few massively expressed genes consuming 20% of reads, all other genes' FPKM values are deflated relative to another sample. Wagner et al. (2012, Theory Biosci) showed that average FPKM varies between samples even for the same genome.

TPM (Transcripts Per Million) partially fixes this by normalizing to a constant sum (1 million per sample), but it still represents relative abundance within a sample's sequenced population, not absolute abundance. Neither TPM nor FPKM should be used for differential expression testing.

| Unit | Within-sample comparable | Between-sample comparable | Use for DE |
|------|-------------------------|--------------------------|------------|
| Raw counts | No | No | Yes (with DE tools) |
| CPM | No (no length correction) | Only with TMM/RLE factors | No |
| FPKM/RPKM | Yes | No | No |
| TPM | Yes | Partially (same composition caveat) | No |
| Normalized counts | No | Yes | Via DE tool |
| VST/rlog | No | Yes | No (use for viz) |

### Computing TPM

**Goal:** Convert raw counts to TPM for within-sample gene expression comparison.

**Approach:** Divide each gene's count by its effective length (in kilobases), then normalize so the column sums to 1 million.

```python
import pandas as pd
import numpy as np

def counts_to_tpm(counts, gene_lengths):
    '''Convert raw counts to TPM. gene_lengths in bp.'''
    rate = counts.div(gene_lengths / 1000, axis=0)  # RPK
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    return tpm
```

```r
counts_to_tpm <- function(counts, gene_lengths) {
    rate <- counts / (gene_lengths / 1000)  # RPK
    tpm <- t(t(rate) / colSums(rate)) * 1e6
    return(tpm)
}
```

## Variance-Stabilizing Transformations

### VST (DESeq2)

**Goal:** Transform counts so that variance is approximately constant across the mean, suitable for PCA, heatmaps, clustering, and distance-based analyses.

**Approach:** Fit the dispersion-mean trend and apply a variance-stabilizing function. The `blind` parameter controls whether the design is used during dispersion estimation.

```r
library(DESeq2)

# blind=TRUE: re-estimates dispersions ignoring design; use for unbiased QC
vsd <- vst(dds, blind=TRUE)

# blind=FALSE: uses dispersions from the full model; recommended for most analyses
# Prevents the transformation from treating design-based differences as noise
vsd <- vst(dds, blind=FALSE)

vst_matrix <- assay(vsd)
```

The `blind` parameter decision:
- `blind=TRUE`: for initial QC where the goal is unbiased sample assessment (does PCA structure match expected groups?)
- `blind=FALSE`: for all downstream use (heatmaps, WGCNA, machine learning) because it preserves biologically meaningful variance

### rlog (DESeq2)

**Goal:** Apply regularized log transformation that handles low counts and variable library sizes better than VST.

**Approach:** Fit a shrinkage model per gene per sample that pulls low-information estimates toward the mean.

```r
# Better than VST when library sizes vary widely (>10-fold range)
# Much slower -- O(genes x samples); impractical for >100 samples
rld <- rlog(dds, blind=FALSE)
rlog_matrix <- assay(rld)
```

| Criterion | VST | rlog |
|-----------|-----|------|
| Speed | Fast | Slow (O(genes x samples)) |
| >100 samples | Recommended | Not recommended |
| Unequal library sizes | Adequate | Better |
| Low-count genes | May be noisy | Better shrinkage |
| Default choice | Yes | Only when library sizes vary >10-fold |

### log-CPM (edgeR)

**Goal:** Produce log-transformed CPM values suitable for visualization, with dampened low-count noise.

**Approach:** Add a scaled prior count before log transformation to prevent log(0) and reduce noise from low-count genes.

```r
library(edgeR)

# Default prior.count=0.25; use prior.count >= 2 for visualization
log_cpm <- cpm(y, log=TRUE, prior.count=2)
```

```python
import numpy as np

def log_cpm(counts, prior_count=2):
    '''Log2 CPM with prior count for noise dampening.'''
    lib_sizes = counts.sum(axis=0)
    # Scale prior proportional to library size (edgeR convention)
    cpm_vals = (counts + prior_count) / (lib_sizes + 2 * prior_count) * 1e6
    return np.log2(cpm_vals)
```

The `prior.count` parameter matters: the default (0.25) is adequate for statistical use, but for heatmaps and PCA, a larger prior (2-5) dampens noise from genes with very few counts, reducing visual artifacts from low-expression noise.

## GC Content and Gene Length Bias

Standard normalization (TMM, RLE) does not correct for sample-specific GC content bias or gene length bias. These biases arise from library preparation (fragmentation, PCR) and create systematic differences in read coverage that correlate with gene properties.

When this matters: gene set enrichment analysis is particularly susceptible. A 2019 PLOS Biology study showed that uncorrected sample-specific gene-length bias causes recurrent false positives in GSEA.

### EDASeq

**Goal:** Remove GC content and gene length bias using within-lane normalization.

**Approach:** Two rounds of full-quantile normalization -- first within GC-content bins per sample, then between samples.

```r
library(EDASeq)

# Requires gene-level GC content
feature_data <- data.frame(gc=gene_gc_content, length=gene_lengths, row.names=rownames(counts))
data <- newSeqExpressionSet(counts=as.matrix(counts), featureData=feature_data, phenoData=coldata)

# Within-lane GC normalization, then between-lane normalization
data_norm <- withinLaneNormalization(data, 'gc', which='full')
data_norm <- betweenLaneNormalization(data_norm, which='full')
normalized_counts <- counts(data_norm)
```

### cqn (Conditional Quantile Normalization)

```r
library(cqn)

# Combines GC and length correction with quantile normalization
cqn_result <- cqn(counts, x=gene_gc_content, lengths=gene_lengths)

# Use offset in edgeR
y$offset <- cqn_result$glm.offset
```

## Single-Cell Normalization

### scran Deconvolution

**Goal:** Estimate per-cell size factors that handle the high zero-inflation of single-cell data.

**Approach:** Pool cells in overlapping windows (sorted by library size), compute pool-level size factors, then deconvolve back to individual cells. Rough pre-clustering prevents mixing distinct cell types with very different transcriptome sizes.

```r
library(scran)
library(scater)

# Pre-cluster to avoid mixing cell types
clusters <- quickCluster(sce)

# Deconvolution-based size factors
sce <- computeSumFactors(sce, clusters=clusters)

# Apply normalization
sce <- logNormCounts(sce)
```

```python
import scanpy as sc

# scanpy's normalize_total is simpler but less robust than scran
# For rigorous analysis, use scran via rpy2 or anndata2ri
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

Standard bulk normalization (TMM, RLE) fails on single-cell data because many genes have zero counts due to dropout, violating the assumption that most genes are detected in most samples. scanpy's `normalize_total` is a simple CPM-like normalization that works for exploratory analysis but does not handle composition bias as well as scran.

## Pre-Filtering Before Normalization

**Goal:** Remove genes with insufficient expression to reduce noise and multiple testing burden.

**Approach:** Apply expression thresholds that account for experimental design (group sizes).

### edgeR: filterByExpr (Design-Aware)

```r
library(edgeR)

# Automatically determines thresholds based on design and library sizes
# min.count=10, min.total.count=15 by default
# Uses smallest group size from design matrix
keep <- filterByExpr(y, design=model.matrix(~condition, data=metadata))
y <- y[keep, , keep.lib.sizes=FALSE]
```

`filterByExpr` is design-aware: it requires a gene to have CPM above a threshold in at least n samples, where n is the smallest group size. This is more principled than arbitrary thresholds.

### DESeq2: Independent Filtering

DESeq2 performs automatic independent filtering at the `results()` step, optimizing a mean-expression cutoff to maximize adjusted p-value rejections. Manual pre-filtering (`rowSums(counts(dds)) >= 10`) is purely for computational speed and memory -- it does not affect statistical results because independent filtering handles it downstream.

```r
# Speed-only pre-filter (does not affect statistical results)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# The alpha in results() should match the planned FDR cutoff
# because independent filtering is optimized at this threshold
res <- results(dds, alpha=0.05)
```

Key distinction: edgeR requires explicit filtering (via `filterByExpr`) because its quasi-likelihood pipeline does not have DESeq2's independent filtering. Forgetting `filterByExpr` in edgeR inflates the multiple testing burden.

## Common Mistakes

| Mistake | Why it fails | Correct approach |
|---------|-------------|-----------------|
| Providing normalized data to DESeq2/edgeR | DE tools model raw count distributions; normalized input violates assumptions | Always use raw integer counts |
| Using FPKM for cross-sample DE | Composition bias makes FPKM incomparable between samples | Use TMM/RLE-normalized counts via DE tools |
| Using VST/rlog values for DE testing | Transformed values lose the count distribution properties DE tools need | VST/rlog for visualization only |
| Applying log2(counts+1) globally | Pseudocount of 1 is arbitrary; does not stabilize variance properly | Use VST, rlog, or edgeR's log-CPM with prior.count |
| Normalizing single-cell data with TMM | High dropout violates TMM assumptions | Use scran deconvolution or scanpy normalize_total |
| Forgetting `filterByExpr` in edgeR | Inflated multiple testing burden | Always filter before `estimateDisp` |
| Using `blind=TRUE` for downstream VST | Over-shrinks transformed values by ignoring known design | Use `blind=FALSE` for heatmaps, WGCNA, ML |

## Related Skills

- expression-matrix/counts-ingest - Load count data before normalization
- differential-expression/deseq2-basics - DESeq2 handles normalization internally
- differential-expression/edger-basics - edgeR TMM normalization
- single-cell/preprocessing - Single-cell normalization workflows
- rna-quantification/count-matrix-qc - QC before normalization
