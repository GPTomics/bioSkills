---
name: bio-small-rna-seq-differential-mirna
description: Tests miRNAs for differential expression with DESeq2 or edgeR using small-RNA-aware normalization and filtering. Use when deciding which normalization survives a library dominated by a few hyper-abundant miRNAs (compositional fragility); choosing DESeq2 vs edgeR vs a compositional method; setting a lower prefilter than mRNA; handling biofluid data with no endogenous normalizer; or remembering that RPM is for display and TDMD can make a miRNA drop without transcriptional repression.
tool_type: r
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, apeglm 1.24+, EnhancedVolcano 1.20+, pheatmap 1.0.12+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential miRNA Expression

**"Find differentially expressed miRNAs between my conditions"** -> Test a raw miRNA count matrix for expression changes, accounting for the compositional fragility that makes miRNA normalization harder than mRNA.
- R: `DESeq2::DESeq()` or `edgeR::glmQLFTest()` on RAW miRNA counts

## The governing principle: a few miRNAs dominate the library, so normalization is the dominant decision

A miRNA library is not a gently varying pool of thousands of features like an mRNA library. A handful of tissue-dominant miRNAs can be more than half of all reads, and the expressed repertoire is only hundreds to low-thousands of miRNAs. Two consequences follow, and they matter more than the choice of DE engine. First, global-scaling normalizers (DESeq2 median-of-ratios, edgeR TMM) assume most features are not differentially expressed and the count distribution is roughly symmetric; when one dominant miRNA shifts between conditions it absorbs the size factor and distorts every other miRNA's normalized value, manufacturing phantom changes. The normalization choice genuinely changes which miRNAs are called DE (Garmire 2012; Tam 2015) - so filter low-count noise FIRST, inspect whether a few miRNAs dominate, and report the normalizer. Second, empirical-Bayes dispersion shrinkage borrows strength across features, so with only hundreds of miRNAs the prior is estimated from a small, noisy population and is weaker than on ~20k genes; apeglm LFC shrinkage matters more for the many low-count miRNAs.

Two reframes prevent classic mistakes. RPM is for display and cross-sample viewing, never for testing - hand RAW counts to DESeq2/edgeR, which model the count distribution themselves. And a miRNA going DOWN does not necessarily mean transcriptional repression: target-directed miRNA degradation (TDMD, via ZSWIM8) lets a highly complementary target trigger decay of the miRNA itself (Han 2020; Shi 2020), so interpret a drop as a change in steady-state level, not automatically as reduced biogenesis.

A third decision is the level of testing. Mature-miRNA-level DE answers "which miRNAs changed" with good power; isomiR-level DE is sparser (more features and zeros, weaker per-feature power, heavier multiplicity), and 5' isomiRs shift the seed and can move OPPOSITE to the canonical mature form - so never silently sum 5' isomiRs into the mature count. Collapse to mature for the standard question; test at isomiR resolution only when isomiR identity is the biology.

## Decision: which normalization / method

| Method | Normalization assumption | Best when | Fails when |
|--------|--------------------------|-----------|------------|
| DESeq2 (median-of-ratios) | most features stable; symmetric | balanced designs, no single runaway miRNA | one miRNA dominates and shifts (compositional) |
| edgeR TMM (glmQLF) | most features stable; trimmed mean | similar to DESeq2; flexible GLM | strong composition shift; default 30%/5% trim built for thousands of mRNAs |
| upper-quartile / quantile / Lowess | rank/quantile-based | skewed miRNA distributions (often better-behaved per Garmire) | when the global shape itself is the biology |
| spike-in (cel-miR-39) | external technical scale | biofluids with no endogenous reference; controls extraction | does not correct ligation bias or biological composition |
| RUVg (RUVSeq) | unwanted variation from control miRNAs | hidden batch/technical structure global scaling misses | controls poorly chosen |
| CLR + ALDEx2 (compositional) | treat counts as compositional | as a sensitivity analysis when a few miRNAs dominate | still blind to a global pool shift; more conservative |

When a perturbation moves the WHOLE pool (e.g. Dicer/Drosha loss), every internal normalizer - including CLR - forces the average change to zero and is blind to it; only external spike-ins or cell-number normalization detect a global shift (Lovén 2012).

## Load the count matrix

**Goal:** Read raw miRNA counts and build sample metadata for testing.

**Approach:** Load the miRge3/miRDeep2 count CSV (raw, not RPM) and define the condition factor.

```r
library(DESeq2)

counts <- read.csv('miR.Counts.csv', row.names = 1)   # RAW counts, not RPM
coldata <- data.frame(
    condition = factor(c('control', 'control', 'treated', 'treated')),
    row.names = colnames(counts))
```

## DESeq2 analysis

**Goal:** Identify miRNAs that change between conditions with small-RNA-aware filtering and shrinkage.

**Approach:** Build a DESeqDataSet from rounded raw counts, prefilter at a lower threshold than mRNA, run DESeq2, then shrink LFCs with apeglm for the many low-count miRNAs.

```r
dds <- DESeqDataSetFromMatrix(
    countData = round(counts),     # DESeq2 needs integers
    colData = coldata,
    design = ~ condition)

# Lower prefilter than mRNA: miRNA libraries have fewer total counts, and most
# miRBase entries are near-zero noise. Justify the threshold; do not test everything.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Inspect for compositional risk: a size factor far from 1, or one miRNA that is a
# large fraction of reads, is a warning that median-of-ratios may be distorted.
sizeFactors(dds)

res <- results(dds, contrast = c('condition', 'treated', 'control'))
# apeglm shrinks via a named coef; for an arbitrary/multi-level contrast not expressible
# as one coef, use type = 'ashr' instead.
res_shrunk <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
res_shrunk <- res_shrunk[order(res_shrunk$padj), ]
```

## edgeR alternative

**Goal:** Test the same data with edgeR's quasi-likelihood GLM as a cross-check.

**Approach:** Build a DGEList, filter with filterByExpr, TMM-normalize, estimate dispersion, and run the QL F-test.

```r
library(edgeR)

dge <- DGEList(counts = round(counts), group = coldata$condition)
keep <- filterByExpr(dge, group = coldata$condition)   # pass group or it treats all samples as one
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)                            # TMM

design <- model.matrix(~ condition, data = coldata)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
res_edger <- topTags(qlf, n = Inf)$table               # edgeR uses $FDR, not $padj
```

## Report effect size and expression level, not just FDR

**Goal:** Avoid calling low-count miRNAs DE on the strength of unstable fold-changes.

**Approach:** Filter on shrunk LFC and FDR, but always inspect base mean / CPM, because a significant LFC on a ~5-count miRNA is almost always noise.

```r
sig <- subset(as.data.frame(res_shrunk), padj < 0.05 & abs(log2FoldChange) > 1)
sig$baseMean <- res_shrunk[rownames(sig), 'baseMean']  # keep expression level visible
sig <- sig[order(sig$padj), ]
```

## Visualize

**Goal:** Show the result with a volcano plot and a heatmap of significant miRNAs.

**Approach:** Use EnhancedVolcano on the shrunk results and a variance-stabilized, row-scaled pheatmap.

```r
library(EnhancedVolcano); library(pheatmap)

EnhancedVolcano(res_shrunk, lab = rownames(res_shrunk),
    x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 1,
    title = 'Differential miRNA expression')

# vst() subsets 1000 genes to fit the dispersion trend and ERRORS on miRNA-sized data
# (hundreds of features) - use the full varianceStabilizingTransformation instead.
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vsd)[rownames(sig), , drop = FALSE]
pheatmap(t(scale(t(mat))), annotation_col = coldata['condition'],
    show_rownames = nrow(mat) < 50)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Everything looks DE in one direction | One dominant miRNA shifted and distorted the size factors | Filter first; inspect sizeFactors; try upper-quartile/quantile or remove the runaway from size-factor estimation |
| Inflated significance on tiny miRNAs | RPM (or unfiltered low counts) fed to the test | Use RAW counts and a lower prefilter; report baseMean for every call |
| `filterByExpr` warns "all samples one group" | group/design not passed | `filterByExpr(dge, group = coldata$condition)` |
| edgeR results have no `padj` column | edgeR names the FDR column `FDR` | Use `topTags(...)$table$FDR`, not `$padj` |
| Biofluid DE driven by a few samples | hemolysis/batch confound; no endogenous normalizer | Add cel-miR-39 spike-in normalization; flag hemolysis (miR-451a:miR-23a-3p); model batch |
| A known miRNA "down" but its gene is unchanged | TDMD (target-driven degradation), not transcription | Interpret as steady-state change; check pri/pre-miRNA or ZSWIM8 context before claiming repression |
| `vst()` errors "less than 'nsub' rows" | vst() subsets 1000 genes; miRNA datasets have only hundreds | Use `varianceStabilizingTransformation(dds, blind=FALSE)` (full VST) instead of `vst()` |

## Related Skills

- mirge3-analysis - Produces the raw count matrix
- mirdeep2-analysis - Alternative quantification
- target-prediction - Predict and validate targets of DE miRNAs
- differential-expression/deseq2-basics - General DESeq2 mechanics
- differential-expression/edger-basics - General edgeR mechanics

## References

- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550. doi:10.1186/s13059-014-0550-8
- Robinson MD, McCarthy DJ, Smyth GK. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26:139-140. doi:10.1093/bioinformatics/btp616
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35:2084-2092. doi:10.1093/bioinformatics/bty895
- Garmire LX, Subramaniam S. 2012. Evaluation of normalization methods in mammalian microRNA-Seq data. *RNA* 18:1279-1288. doi:10.1261/rna.030916.111
- Tam S, Tsao MS, McPherson JD. 2015. Optimization of miRNA-seq data preprocessing. *Brief Bioinform* 16:950-963. doi:10.1093/bib/bbv019
- Han J, LaVigne CA, Jones BT, et al. 2020. A ubiquitin ligase mediates target-directed microRNA decay independently of tailing and trimming. *Science* 370:eabc9546. doi:10.1126/science.abc9546
- Shi CY, Kingston ER, Kleaveland B, et al. 2020. The ZSWIM8 ubiquitin ligase mediates target-directed microRNA degradation. *Science* 370:eabc9359. doi:10.1126/science.abc9359
- Lovén J, Orlando DA, Sigova AA, et al. 2012. Revisiting global gene expression analysis. *Cell* 151:476-482. doi:10.1016/j.cell.2012.10.012
