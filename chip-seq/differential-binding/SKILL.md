---
name: bio-chipseq-differential-binding
description: Identifies differentially bound ChIP-seq regions between conditions using DiffBind (from BAMs), DESeq2, or PyDESeq2 (from count matrices). Handles normalization, statistical testing, and fold-change estimation with ChIP-seq-specific considerations. Use when comparing ChIP-seq binding between experimental conditions.
tool_type: mixed
primary_tool: DiffBind
---

## Version Compatibility

Reference examples tested with: DiffBind 3.12+, DESeq2 1.42+, edgeR 4.0+, PyDESeq2 0.4+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show pydeseq2` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Binding Analysis

**"Compare ChIP-seq binding between conditions"** → Identify genomic regions with statistically significant differences in transcription factor or histone mark occupancy between experimental groups.
- R (from BAMs): `DiffBind::dba()` → `dba.count()` → `dba.normalize()` → `dba.analyze()`
- R (from count matrix): `DESeq2::DESeqDataSetFromMatrix()` → `DESeq()` → `results()`
- Python (from count matrix): `pydeseq2.DeseqDataSet()` → `.deseq2()` → `DeseqStats()`

## Choosing an Approach

| Scenario | Recommended |
|----------|-------------|
| Have BAM + peak files | DiffBind (handles consensus peaks, counting, normalization) |
| Have a count matrix (peaks x samples) | DESeq2 directly (R) or PyDESeq2 (Python) |
| Need background normalization | DiffBind with `background=TRUE` or csaw |
| Suspect global binding changes | Custom size factors or spike-in normalization |
| Pre-defined regions (promoters, enhancers) | DESeq2/PyDESeq2 on counted regions |

DiffBind wraps DESeq2/edgeR internally. When a count matrix already exists, DiffBind adds no statistical value — its main contributions (consensus peak generation, summit re-centering, read counting) are already done.

## ChIP-seq Normalization

Normalization is the single most consequential analytical decision for differential binding. The choice of reference reads matters far more than algorithm (RLE vs TMM produce essentially identical results on the same reference reads).

### Reference Reads

| Reference | How applied | Assumption | Risk |
|-----------|-------------|------------|------|
| Reads in peaks | DESeq2 default `estimateSizeFactors()` on count matrix | Most peaks NOT differentially bound | Fails with global changes; can reverse conclusions |
| Full library | Total aligned reads from BAM | Sequencing depth is main technical variable | Conservative; DiffBind 3.0+ default |
| Background bins | 10-15 kb genomic bins via csaw | Background is stable across conditions | Robust to composition bias; requires BAMs |
| Spike-in | Exogenous chromatin (e.g., Drosophila) | Spike-in ratio is constant | Gold standard for global changes; requires spike-in protocol |

### When Default DESeq2 Normalization Works

Running `estimateSizeFactors()` on a ChIP-seq count matrix computes RLE normalization on reads-in-peaks. This works when:
- Most peaks are expected to be stable (typical TF ChIP-seq with localized changes)
- Changes are roughly symmetric (similar numbers of gained and lost peaks)
- Differential peaks are a minority of the total peak set

### When It Fails

- **Global binding changes**: EZH2 inhibitor reducing H3K27me3, BET inhibitor reducing BRD4, knockdown of the ChIPped factor
- **Asymmetric changes**: predominantly gained or lost, not balanced
- **Small peak sets** where differential peaks are a large fraction of total

**Diagnostic**: if the MA plot loess curve deviates substantially from y=0, normalization is not centering correctly. Verify biological plausibility of the gain/loss ratio.

### Custom Size Factors

When default normalization is inappropriate:

```r
# Full library size normalization (when total aligned reads are known)
library_sizes <- c(15e6, 14e6, 16e6, 15e6, 13e6, 14e6)
names(library_sizes) <- colnames(counts(dds))
sizeFactors(dds) <- library_sizes / mean(library_sizes)

# Use known stable peaks as normalization reference
stable_idx <- which(rownames(dds) %in% known_stable_peaks)
dds <- estimateSizeFactors(dds, controlGenes = stable_idx)
```

## From Count Matrix (R - DESeq2)

**Goal:** Identify differentially bound peaks from a pre-computed count matrix using DESeq2's negative binomial framework.

**Approach:** Load counts into DESeqDataSet, apply minimal ChIP-seq-appropriate filtering, run the DESeq2 pipeline, extract results at the desired significance threshold.

```r
library(DESeq2)

counts <- read.delim('counts.tsv', row.names = 1, check.names = FALSE)
coldata <- data.frame(
    condition = factor(c(rep('treated', 3), rep('control', 3))),
    row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# ChIP-seq peaks are already enriched regions; filter less aggressively than RNA-seq
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = 'control')
dds <- DESeq(dds)

# alpha matches intended significance threshold (optimizes independent filtering)
res <- results(dds, alpha = 0.05)
```

### Pre-filtering for ChIP-seq

Unlike RNA-seq where thousands of genes have zero expression, ChIP-seq peaks were called because they have signal. Aggressive filtering can remove truly differential peaks — a peak present in treatment but absent in control will have near-zero counts in control samples.

| Feature set size | Recommended filter |
|-----------------|-------------------|
| < 500 peaks | `rowSums(counts(dds)) >= 1` (remove only all-zero rows) |
| 500-5,000 peaks | `rowSums(counts(dds)) >= 5` |
| > 5,000 peaks | `rowSums(counts(dds)) >= 10` or standard RNA-seq thresholds |

### The Alpha Parameter

The `alpha` argument in `results()` controls independent filtering optimization. Setting it to match the intended significance threshold (e.g., `alpha = 0.05` when filtering at padj < 0.05) maximizes discoveries at that cutoff. Using the default `alpha = 0.1` while filtering at padj < 0.05 is slightly suboptimal but the effect is usually small.

### Log Fold Change Shrinkage

LFC shrinkage (apeglm, ashr) improves fold-change estimates for ranking and visualization but does NOT change p-values or padj. Always use padj from `results()` for significance calls.

```r
library(apeglm)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
```

## From Count Matrix (Python - PyDESeq2)

**Goal:** Run differential binding analysis in Python using PyDESeq2's implementation of the DESeq2 negative binomial framework.

**Approach:** Load count data into a DeseqDataSet (samples as rows, features as columns), run the pipeline, extract results with BH-adjusted p-values.

```python
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
metadata = pd.DataFrame(
    {'condition': ['treated'] * 3 + ['control'] * 3},
    index=counts.columns
)

inference = DefaultInference(n_cpus=4)
dds = DeseqDataSet(counts=counts.T, metadata=metadata, design='~ condition',
                   refit_cooks=True, inference=inference)
dds.deseq2()

ds = DeseqStats(dds, contrast=['condition', 'treated', 'control'], inference=inference)
ds.summary()
results = ds.results_df
```

PyDESeq2 expects counts with samples as rows and features as columns — transpose a peaks-by-samples matrix with `.T`. The `contrast` argument takes `['variable', 'numerator', 'denominator']`, so positive log2FC means higher in the numerator condition.

## DiffBind Workflow (From BAMs)

**Goal:** Run the complete DiffBind pipeline from BAM and peak files to differential binding results.

**Approach:** Create a sample sheet linking BAMs, peaks, and metadata; load into DiffBind; count reads in consensus peaks; normalize; run statistical testing.

### Sample Sheet

```r
samples <- data.frame(
    SampleID = c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'),
    Tissue = c('cell', 'cell', 'cell', 'cell'),
    Factor = c('H3K4me3', 'H3K4me3', 'H3K4me3', 'H3K4me3'),
    Condition = c('control', 'control', 'treatment', 'treatment'),
    Replicate = c(1, 2, 1, 2),
    bamReads = c('ctrl1.bam', 'ctrl2.bam', 'treat1.bam', 'treat2.bam'),
    Peaks = c('ctrl1_peaks.narrowPeak', 'ctrl2_peaks.narrowPeak',
              'treat1_peaks.narrowPeak', 'treat2_peaks.narrowPeak'),
    PeakCaller = c('macs', 'macs', 'macs', 'macs')
)
write.csv(samples, 'samples.csv', row.names = FALSE)
```

### Load, Count, Normalize, Analyze

```r
library(DiffBind)

dba_obj <- dba(sampleSheet = 'samples.csv')
dba_obj <- dba.count(dba_obj, summits = 250, minOverlap = 2)
dba_obj <- dba.normalize(dba_obj)
dba_obj <- dba.contrast(dba_obj, design = '~ Condition')
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)
```

### DiffBind Normalization Options

```r
# Default: full library size (conservative, recommended starting point)
dba_obj <- dba.normalize(dba_obj)

# Background normalization (robust to composition bias; uses csaw bins)
dba_obj <- dba.normalize(dba_obj, background = TRUE)

# Inspect applied normalization
dba.normalize(dba_obj, bRetrieve = TRUE)
```

| Method | When to use |
|--------|------------|
| `DBA_NORM_LIB` + full library (default) | Most analyses; conservative baseline |
| `DBA_NORM_RLE` + reads in peaks | Most peaks expected stable; typical TF ChIP |
| Background (`background=TRUE`) | Suspected composition bias or global changes |
| Spike-in (`spikein=TRUE`) | ChIP-Rx experiments with exogenous reference |
| Loess offsets (`offsets=TRUE`) | Abundance-dependent efficiency biases (use cautiously — can over-normalize) |

### Results and Export

```r
db_all <- dba.report(dba_obj, th = 1)
results_df <- as.data.frame(db_all)
write.csv(results_df, 'differential_binding.csv', row.names = FALSE)

# Significant peaks at FDR < 0.05
db_sig <- dba.report(dba_obj, th = 0.05)

# Export to BED
library(rtracklayer)
export(db_sig, 'diff_peaks.bed', format = 'BED')
```

### Visualization

```r
dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)
dba.plotMA(dba_obj)
dba.plotVolcano(dba_obj)
dba.plotHeatmap(dba_obj, contrast = 1, correlations = FALSE)
```

## Multi-Factor Design

```r
# DiffBind with blocking factor
dba_obj <- dba.contrast(dba_obj, design = '~ Batch + Condition')
dba_obj <- dba.analyze(dba_obj)

# DESeq2 directly with batch correction
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                               design = ~ batch + condition)
dds <- DESeq(dds)
res <- results(dds, name = 'condition_treated_vs_control')
```

## DiffBind 3.0+ Notes

Key defaults changed in DiffBind 3.0:
- `summits=200` recenters peaks (was `FALSE`); set `summits=FALSE` for broad histone marks
- `dba.normalize()` is now required before analysis
- Blacklist filtering applied by default
- Full library size normalization replaces reads-in-peaks (more conservative)
- Use design formulas instead of `group1`/`group2` for contrasts

## Sample Sheet Columns

| Column | Required | Description |
|--------|----------|-------------|
| SampleID | Yes | Unique identifier |
| Condition | Yes | Experimental condition |
| Replicate | Yes | Replicate number |
| bamReads | Yes | Path to BAM file |
| Peaks | Yes | Path to peak file |
| PeakCaller | Yes | macs, bed, narrow |
| bamControl | No | Path to input BAM |
| Tissue | No | Tissue/cell type |
| Factor | No | ChIP target |

## Related Skills

- peak-calling - Generate input peak files
- peak-annotation - Annotate differential peaks to genes
- differential-expression/deseq2-basics - DESeq2 fundamentals for count-based testing
- differential-expression/edger-basics - Alternative statistical framework
- pathway-analysis/go-enrichment - Functional enrichment of peak-associated genes
