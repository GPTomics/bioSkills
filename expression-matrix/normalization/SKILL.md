---
name: bio-expression-matrix-normalization
description: Normalizes and transforms RNA-seq count matrices for DE, visualization, clustering, and ML. Covers between-sample (TMM, TMMwsp, RLE/median-of-ratios, upper quartile), within-sample (TPM, FPKM/RPKM), variance-stabilizing (VST, rlog, log-CPM), GC-content correction (cqn, EDASeq), and single-cell (scran deconvolution, scanpy normalize_total). Encodes the composition-bias rationale, the "most genes not DE" assumption and its catastrophic failure modes (MYC amplification, apoptosis, viral host shutoff, prokaryotic stress), the "lengthScaledTPM is not TPM" naming trap, the "TPM is not for DE" rule, the blind=TRUE vs FALSE decision, ERCC spike-in normalization (SBN), and the single-cell zero-inflation breakdown of TMM/RLE. Use when choosing or applying normalization, debugging shifted-MA-plot diagnostics, handling zero-heavy single-cell data, or correcting GC bias.
tool_type: mixed
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, limma 3.58+, pandas 2.2+, numpy 1.26+, scanpy 1.10+, scran 1.30+, scater 1.30+, EDASeq 2.36+, cqn 1.48+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Expression Matrix Normalization

**"Normalize my counts for X"** -> Pick a method that matches the downstream task (raw counts for DE; TMM/RLE-scaled CPM for cross-sample comparison; VST/rlog for visualization/ML; scran for single-cell) and the data structure (bulk vs single-cell, zero-heavy or not, length-biased or not).

## The Single Most Important Modern Insight -- TMM/RLE assume most genes are NOT DE; that assumption is wrong in MYC, apoptosis, viral host shutoff, and prokaryotic stress

TMM (Robinson & Oshlack 2010 *Genome Biol* 11:R25) and RLE/median-of-ratios (Anders & Huber 2010 *Genome Biol* 11:R106) both rely on the assumption that the MAJORITY of genes are unchanged between samples. The scaling factor is computed from a trimmed reference. When that assumption holds (most bulk RNA-seq), both methods are excellent.

When the assumption fails catastrophically:

| Biology | Mechanism | Symptom |
|---------|-----------|---------|
| MYC amplification | MYC drives global 2-3x transcriptional amplification | Scaling factors absorb the global shift; reported LFCs muted |
| Apoptosis / cell death | Massive transcriptional shutdown | Surviving (often mitochondrial) transcripts appear up-regulated |
| Viral host shutoff (HSV-1, vaccinia) | Host mRNA degraded by viral nucleases | Apparent up-regulation of non-degraded host genes |
| Transcription / splicing inhibitors (DRB, flavopiridol) | Pol II elongation blocked | Similar global shutdown signature |
| Cell-cycle synchronization | G1 vs S vs M differ in total RNA per cell | Per-cell RNA shifts; library-size-only normalization mis-corrects |
| Prokaryotic stress | Bacteria rewire large fractions of transcriptome | >50% of genes truly DE; trimmed mean is no longer "background" |

Detection: MA plot shows the bulk cloud clearly shifted off zero; reported fold changes don't match qPCR or Western for known-DE genes.

Fix: ERCC spike-in normalization (Jiang 2011 *Genome Res* 21:1543) -- 96 synthetic RNAs spiked proportional to cell number; use spike-in counts as the size factor. Or: `controlGenes=` with a curated stable housekeeping set. The "most genes not DE" assumption is checkable; check it on every dataset that might violate it.

A second insight: **`lengthScaledTPM` from tximport is NOT TPM**. It is a count-scale matrix (sums to library size) with length-bias removed, designed as input to limma-voom (which cannot accept offsets). The "TPM" in the option name has misled many users into reporting these values as normalized abundance. See `expression-matrix/counts-ingest` for the full `countsFromAbundance` decision tree.

A third insight: **VST/rlog values are NEVER input to DE**. They are for PCA, heatmaps, clustering, ML features. DE tools model the count distribution directly; passing VST/rlog values silently violates their assumptions.

## Normalization Decision Table

| Task | Method | Tool | Input |
|------|--------|------|-------|
| DE testing (DESeq2) | RLE / median-of-ratios | `DESeq()` internally | Raw integer counts |
| DE testing (edgeR) | TMM / TMMwsp | `normLibSizes()` internally | Raw integer counts |
| DE testing (limma-voom) | TMM + voom precision weights | `normLibSizes()` then `voom()` | Raw integer counts |
| PCA, heatmaps, clustering (n > 30) | VST | `vst(dds, blind = FALSE)` | DESeqDataSet |
| PCA, heatmaps, clustering (n < 30, library sizes vary >4x) | rlog | `rlog(dds, blind = FALSE)` | DESeqDataSet |
| PCA, heatmaps (edgeR/limma) | log-CPM | `cpm(y, log=TRUE, prior.count=2)` | DGEList |
| WGCNA | VST or log-CPM | `vst(blind=FALSE)` | DESeqDataSet |
| GSVA / ssGSEA | log2(TPM+1) or VST | precomputed | TPM or DESeqDataSet |
| ML biomarker model | VST | `vst(blind=FALSE)` | DESeqDataSet |
| Cross-sample expression reporting | DESeq2 normalized counts | `counts(dds, normalized=TRUE)` | DESeqDataSet |
| Within-sample gene ranking | TPM | quantification tool output | -- |
| Single-cell normalization | scran deconvolution | `computeSumFactors` + `logNormCounts` | SingleCellExperiment |
| Single-cell, simple/exploratory | log1p(CPM-like) | `sc.pp.normalize_total + sc.pp.log1p` | AnnData |
| Cross-sample but composition-shifted (majority-DE biology) | Spike-in (ERCC) or `controlGenes=` | DESeq2/edgeR with offset | Raw counts + controls |
| GC-content bias (cross-platform integration) | cqn or EDASeq | `cqn()` returns offset | Counts + per-gene GC + length |

## Between-Sample Normalization

### RLE / Median of Ratios (DESeq2)

**Goal:** Estimate per-sample size factors that correct for library size and composition bias.

**Approach:** Geometric mean per gene across samples forms a pseudo-reference; per-sample size factor is the median ratio of (sample count / reference count). Median is robust to a minority of DE genes.

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)
```

Size factor interpretation: 1.2 means the sample has 20% more sequencing depth (after composition adjustment) than the reference.

The geometric mean is undefined for any gene with a zero in any sample, so RLE `type='ratio'` (default) silently EXCLUDES genes with any zero from the reference. For single-cell or zero-heavy data, use:

```r
dds <- estimateSizeFactors(dds, type = 'poscounts')
```

`poscounts` uses only positive entries per gene, salvaging the geometric mean.

| Scenario | Use |
|----------|-----|
| Standard bulk RNA-seq | Default `type='ratio'` |
| Zero-heavy (single-cell, sparse) | `type='poscounts'` |
| Very small libraries | `type='iterate'` |
| Known stable reference genes | `controlGenes=stable_idx` |
| Majority-DE biology (stress, MYC, viral) | `controlGenes=` or spike-in SBN |

### TMM / TMMwsp (edgeR)

**Goal:** Compute normalization factors that account for composition bias via trimmed mean of M-values.

**Approach:** Select reference; compute gene-wise log-ratios (M-values) and average expression (A-values); trim extremes; weighted mean scaling factor.

```r
library(edgeR)

y <- DGEList(counts = counts, group = coldata$condition)
y <- normLibSizes(y)
y$samples$norm.factors
```

`normLibSizes()` is the v4 canonical name (was `calcNormFactors()` in v3; same function). `method='TMM'` remains the documented default; `method='TMMwsp'` (TMM with singleton pairing) is an alternative for samples with many zeros and is the preferred choice for sparse / single-cell-pseudobulk data. Pass `method=` explicitly for reproducibility. Both old name and method still work.

TMM defaults: `logratioTrim = 0.3` (trim top 30% and bottom 30% of M values), `sumTrim = 0.05`. The factor enters the GLM as part of the offset; counts are NOT divided.

For visualization:

```r
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
cpm_vis <- cpm(y, normalized.lib.sizes = TRUE)
```

`prior.count = 2` is the modern edgeR default for `log=TRUE`. Smaller priors (0.25) leave low-count log values noisy; larger priors (5-10) shrink them further toward zero. limma-trend assumes the log-CPM uses cpm log=TRUE.

### Upper Quartile

```r
y <- normLibSizes(y, method = 'upperquartile')
```

Per-sample 75th percentile of non-zero counts (Bullard et al. 2010 *BMC Bioinformatics* 11:94). Robust to the "most genes unchanged" assumption -- useful when TMM/RLE fail. Less common in modern practice; usable for microRNA-seq or targeted panels where TMM doesn't have enough genes to trim.

### Spike-In Normalization (ERCC)

```r
ercc_idx <- grep('^ERCC-', rownames(counts))
y <- DGEList(counts = counts[-ercc_idx, ], group = coldata$condition)
ercc_lib <- colSums(counts[ercc_idx, ])
y$samples$norm.factors <- ercc_lib / mean(ercc_lib)
```

Spike-in based normalization (SBN) uses ERCC counts as the per-sample size factor proxy, decoupling normalization from endogenous transcript composition. Required for MYC, viral, or other majority-DE biology. Caveat: spike-ins have technical variance; the spike-in:endogenous ratio depends on input cell number being accurately measured (often violated).

## Within-Sample Normalization (TPM, FPKM)

| Unit | Definition | Within-sample comparable | Between-sample comparable | Use for DE |
|------|------------|-------------------------|--------------------------|------------|
| Raw counts | reads per gene | No | No | YES (via DE tools) |
| CPM | reads / library_size * 1e6 | No (no length correction) | Only with TMM/RLE factors | No |
| RPKM/FPKM (Mortazavi 2008) | reads / (length_kb * lib_size_in_M) | Yes | No (composition bias) | No |
| TPM (Wagner 2012) | (reads/length) / sum(reads/length) * 1e6 | Yes | Partially (same composition caveat) | No |
| Normalized counts | DESeq2 / edgeR scaled | No | Yes | Via DE tool |
| VST/rlog | DESeq2 stabilized | No | Yes | NO (visualization only) |

```python
import pandas as pd

def counts_to_tpm(counts, gene_lengths):
    '''Convert raw counts to TPM. gene_lengths in bp.'''
    rate = counts.div(gene_lengths / 1000, axis=0)
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    return tpm
```

```r
counts_to_tpm <- function(counts, gene_lengths) {
    rate <- counts / (gene_lengths / 1000)
    t(t(rate) / colSums(rate)) * 1e6
}
```

TPM's sum-to-1e6-per-sample makes cross-sample comparisons dimensionally coherent but does NOT fix composition shifts -- TPM remains a COMPOSITIONAL value. Under massive global changes, TPM destroys the absolute-scale signal.

Wagner GP, Kin K, Lynch VJ 2012 *Theory Biosci* 131:281 is the canonical "why TPM > RPKM" reference. Their headline: "average FPKM varies between samples even for the same genome."

CRITICAL: do NOT use TPM as input to DESeq2 / edgeR / limma-voom. They model the count distribution directly. TPM is for: within-sample ranking, gene-set comparison within sample, deconvolution (CIBERSORT, EPIC, MCP-counter all require TPM), and reporting expression levels.

## Variance-Stabilizing Transformations (VST, rlog)

**Goal:** Transform counts so variance is approximately constant across the mean, suitable for PCA, heatmaps, clustering, ML features.

**Approach:** Fit dispersion-mean trend and apply variance-stabilizing function.

```r
library(DESeq2)

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

vst_matrix <- assay(vsd)
```

| Criterion | VST | rlog |
|-----------|-----|------|
| Speed | Fast (1 sec for thousands of samples) | Slow (30+ sec for 100 samples) |
| n > 30 | Recommended | Impractical |
| Unequal library sizes (>4x) | Adequate | Better (more shrinkage) |
| Low-count genes | May be noisy | Better shrinkage |
| Default choice | YES | Only when n<30 AND lib sizes vary >4x |

`blind=TRUE` (vst default) re-estimates dispersions IGNORING the design -- appropriate for unbiased QC ("are samples consistent independent of design?"). `blind=FALSE` uses fitted dispersions -- appropriate for downstream visualization after the design is settled. Modern DESeq2 vignette recommends `blind=FALSE` for any plot AFTER the model is fit.

`vst()` uses 1000 most-variable genes by default to fit the dispersion trend. With <1000 genes after filtering, set `nsub` lower.

CRITICAL: VST/rlog values are for visualization only. NEVER pass them to `DESeq()` or `glmQLFit()`. DE tools require raw counts.

### log-CPM (edgeR / limma)

```r
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
```

```python
import numpy as np
def log_cpm(counts, prior_count=2):
    lib_sizes = counts.sum(axis=0)
    cpm_vals = (counts + prior_count) / (lib_sizes + 2 * prior_count) * 1e6
    return np.log2(cpm_vals)
```

`prior.count = 2` (edgeR modern default) shrinks low-count log values toward zero, reducing visual artifacts from low-expression noise. For statistical use (limma-trend), this default is appropriate; for heatmaps and PCA, larger priors (3-5) further dampen noise.

voom uses 0.5 added per cell with `(lib.size + 1)` denominator -- not strictly equivalent to `cpm(log=TRUE, prior.count=0.5)`. Its mean-variance trend handles low-count noise via per-observation weights. Do not confuse voom's internal pre-log handling with `cpm(log=TRUE)`.

## GC Content and Length Bias (cqn, EDASeq)

Standard TMM/RLE do NOT correct sample-specific GC content bias or gene length bias. These biases arise from library preparation (fragmentation, PCR) and create systematic differences in read coverage correlated with gene properties.

Most affected: GSEA. Sample-specific gene-length bias causes recurrent false positives.

EDASeq (Risso, Schwartz, Sherlock, Dudoit 2011 *BMC Bioinformatics* 12:480) applies within-lane GC normalization then between-lane normalization:

```r
library(EDASeq)

fd <- data.frame(gc = gene_gc, length = gene_lengths, row.names = rownames(counts))
data <- newSeqExpressionSet(as.matrix(counts), featureData = fd, phenoData = coldata)

data_norm <- withinLaneNormalization(data, 'gc', which = 'full')
data_norm <- betweenLaneNormalization(data_norm, which = 'full')
normalized_counts <- counts(data_norm)
```

cqn (Hansen, Irizarry, Wu 2012 *Biostatistics* 13:204) fits a smooth conditional-quantile model on log2 expression as a function of GC and length; the offset enters the GLM directly:

```r
library(cqn)
cqn_res <- cqn(counts, x = gene_gc, lengths = gene_lengths)
y$offset <- cqn_res$glm.offset
```

Use when comparing samples sequenced on different platforms or with different library prep chemistries. Run BEFORE TMM/RLE (or supply the cqn offset to edgeR/DESeq2 directly).

## Single-Cell Normalization

### scran Deconvolution

**Goal:** Per-cell size factors that handle the high zero-inflation of single-cell data without violating TMM/RLE assumptions.

**Approach:** Pool cells in overlapping windows; compute pool-level size factors; deconvolve back to individual cells. Pre-cluster to avoid mixing cell types with very different transcriptome sizes.

```r
library(scran)
library(scater)

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)
```

`quickCluster` provides a rough partition; without it, mixing a cell with 200 detected genes (resting T cell) and 5000 detected genes (activated macrophage) produces wrong size factors.

### scanpy normalize_total

```python
import scanpy as sc

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

Simpler than scran (CPM-like with target_sum=1e4) but less robust to composition bias. Adequate for exploratory analysis; for rigorous DE (especially pseudobulk), use scran-style or aggregate to bulk and use DESeq2.

Standard bulk TMM/RLE FAIL on single-cell data because many genes have zero counts due to dropout, violating "most genes detected in most samples." Use scran or `poscounts`.

## Pre-Filtering Before Normalization

```r
library(edgeR)
keep <- filterByExpr(y, design = model.matrix(~ condition, coldata))
y <- y[keep, , keep.lib.sizes = FALSE]
```

`filterByExpr` is design-aware: requires CPM above threshold in at least n samples, where n is the smallest group size from the design. See `differential-expression/edger-basics` for internals.

DESeq2 performs automatic independent filtering at `results()` time. Manual pre-filtering for DESeq2 is for SPEED only -- it does not affect statistical results.

```r
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
```

Key distinction: edgeR REQUIRES explicit `filterByExpr` (no automatic independent filtering); DESeq2 does both. Forgetting `filterByExpr` in edgeR inflates the multiple-testing burden -- the single biggest reason an edgeR analysis underperforms a comparable DESeq2 analysis.

## Per-Method Failure Modes

### TMM/RLE fail under MYC amplification

**Trigger:** Cancer vs normal with MYC amplification; MA plot shows the bulk cloud clearly above zero; reported fold changes 30-50% smaller than qPCR or Western.

**Mechanism:** MYC drives global ~2x mRNA increase. TMM/RLE assume most genes unchanged -- the assumption is wrong; trimmed reference is dominated by genuinely up-regulated genes; size factors absorb the global shift.

**Symptom:** Known up-regulated genes show muted LFC; spike-in or qPCR shows the real magnitude.

**Fix:** Use ERCC spike-ins for SBN; or supply curated stable housekeeping genes via `controlGenes=` or as `y$offset`.

### Used VST/rlog values as DE input

**Trigger:** PCA looked good on VST; user passed VST matrix to limma `lmFit`; many DE genes.

**Mechanism:** DE tools require raw counts. VST/rlog values are log-scale and homoskedastic by design -- the count-distribution assumptions are violated.

**Symptom:** P-value histogram non-uniform; gene lists don't replicate.

**Fix:** Always use raw counts as input to DE tools. VST/rlog is for visualization, clustering, ML -- not testing.

### Single-cell with TMM -- zero-inflation breaks it

**Trigger:** scRNA-seq data passed through `normLibSizes(y, method='TMM')`; many cells have NA size factors or extreme values.

**Mechanism:** TMM trims top and bottom 30% of M-values. With many zeros, the M-value distribution is dominated by undefined values.

**Symptom:** Size factor estimation fails or produces extreme outliers.

**Fix:** scran deconvolution (with `quickCluster` pre-clustering). For DESeq2 on pseudobulk, `type='poscounts'`.

### TPM used for DE

**Trigger:** Pipeline computed TPM in Python; user wrote `DESeqDataSetFromMatrix(round(tpm), coldata, ...)`.

**Mechanism:** TPM is compositional within sample; cross-sample comparisons are fractions, not abundances. DESeq2 models the count distribution; rounded TPM is neither raw counts nor a meaningful transformation.

**Symptom:** Bizarre dispersion estimates; DE list dominated by housekeeping shifts.

**Fix:** Always use raw counts as input. For Salmon/kallisto output, use tximport.

### blind=TRUE for downstream visualization

**Trigger:** `vst(dds, blind = TRUE)` for the results-figure PCA; clusters look weaker than expected.

**Mechanism:** `blind=TRUE` ignores the design when fitting dispersions, treating biological signal as noise -- appropriate for QC, suboptimal for results figures.

**Symptom:** PCA shows less separation than `blind=FALSE` would.

**Fix:** Use `blind=FALSE` for downstream visualization AFTER the design is settled.

### Quantile normalization across RNA-seq samples with global shift

**Trigger:** Multi-platform integration; user applied quantile normalization to harmonize; downstream DE shows fewer hits than expected.

**Mechanism:** Quantile normalization forces every sample's distribution to be identical by rank-replacing. Assumes the true expression distribution is the same across samples (valid for microarrays with bounded dynamic range; often false for RNA-seq with global shifts).

**Symptom:** Real biological signal erased; especially harmful when biological condition truly shifts the distribution.

**Fix:** TMM/RLE for between-sample composition correction; quantile only for cross-platform integration where dynamic range mismatch dominates.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| Size factors all NA | Zero-heavy data with default `type='ratio'` | `type='poscounts'` |
| MA plot bulk cloud shifted off zero | TMM/RLE assumption violated | Spike-in normalization; `controlGenes=` |
| `lengthScaledTPM` reported as TPM in methods | Misleading option name | Clarify: count-scale matrix with length bias removed; NOT TPM |
| DE list reproducibly different across normalization methods | Marginal cases sensitive to normalization choice | Report results using two methods; flag genes that change rank |
| `cpm(y, log=TRUE)` very noisy at low counts | Default `prior.count` too small | `prior.count = 2` (modern edgeR default) |
| VST produces NA matrix | Too few genes after filtering | Set `nsub` lower in `vst(dds, nsub=...)` |
| Quantile-normalized RNA-seq missing known DE | Quantile erases the global biological shift | Use TMM/RLE; reserve quantile for cross-platform integration |

## References

- Anders S, Huber W. 2010. Differential expression analysis for sequence count data. *Genome Biol* 11(10):R106. doi:10.1186/gb-2010-11-10-r106
- Robinson MD, Oshlack A. 2010. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol* 11(3):R25. doi:10.1186/gb-2010-11-3-r25
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Bullard JH, Purdom E, Hansen KD, Dudoit S. 2010. Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. *BMC Bioinformatics* 11:94. doi:10.1186/1471-2105-11-94
- Wagner GP, Kin K, Lynch VJ. 2012. Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory Biosci* 131(4):281-285. doi:10.1007/s12064-012-0162-3
- Mortazavi A, Williams BA, McCue K, Schaeffer L, Wold B. 2008. Mapping and quantifying mammalian transcriptomes by RNA-Seq. *Nat Methods* 5(7):621-628. doi:10.1038/nmeth.1226
- Bolstad BM, Irizarry RA, Astrand M, Speed TP. 2003. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. *Bioinformatics* 19(2):185-193. doi:10.1093/bioinformatics/19.2.185
- Smid M et al. 2018. Gene length corrected trimmed mean of M-values (GeTMM) processing of RNA-seq data performs similarly in intersample analyses while improving intrasample comparisons. *BMC Bioinformatics* 19:236. doi:10.1186/s12859-018-2246-7
- Lun ATL, Bach K, Marioni JC. 2016. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. *Genome Biol* 17:75. doi:10.1186/s13059-016-0947-7
- Jiang L et al. 2011. Synthetic spike-in standards for RNA-seq experiments. *Genome Res* 21(9):1543-1551. doi:10.1101/gr.121095.111
- Risso D, Schwartz K, Sherlock G, Dudoit S. 2011. GC-content normalization for RNA-Seq data. *BMC Bioinformatics* 12:480. doi:10.1186/1471-2105-12-480
- Hansen KD, Irizarry RA, Wu Z. 2012. Removing technical variability in RNA-seq data using conditional quantile normalization. *Biostatistics* 13(2):204-216. doi:10.1093/biostatistics/kxr054
- Chen Y et al. 2025. edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. *Nucleic Acids Res* 53(2):gkaf018. doi:10.1093/nar/gkaf018

## Related Skills

- counts-ingest - tximport offsets vs scaledTPM vs lengthScaledTPM; biotype pre-filtering before normalize
- gene-id-mapping - Filtering rRNA/Mt biotypes affects normalization
- metadata-joins - Sample alignment before normalization
- sparse-handling - Single-cell sparse matrix normalization patterns
- differential-expression/deseq2-basics - DESeq2 RLE / median-of-ratios internals
- differential-expression/edger-basics - edgeR TMM / TMMwsp internals
- differential-expression/batch-correction - Normalization vs batch correction distinction
- differential-expression/de-visualization - VST/rlog blind choice; log-CPM prior count
- single-cell/preprocessing - Single-cell normalization workflows
- rna-quantification/count-matrix-qc - QC before normalization
