# Differential Abundance - Usage Guide

## Overview
Identify proteins with significantly different abundance between experimental conditions using statistical testing, multiple testing correction, and fold change shrinkage. Covers the full pipeline from raw intensities through preprocessing, statistical modeling, and accurate effect size estimation.

## Prerequisites
```bash
pip install numpy pandas scipy statsmodels
```
```r
BiocManager::install(c("limma", "DEqMS", "ashr", "proDA"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differentially abundant proteins between treatment and control in my intensity matrix"
- "Run limma analysis on my protein data with log2 transformation and median normalization"
- "Identify significant proteins with shrunk fold change estimates"
- "Perform differential abundance testing on my label-free proteomics data"

## Example Prompts

### Full Pipeline
> "I have raw protein intensities in a TSV file with samples as columns. Log2 transform, median normalize, and run differential abundance testing between case and control groups. Report fold changes with shrinkage applied."

> "Analyze my TMT proteomics data for differential abundance between treatment and control. I have PSM counts per protein, so use DEqMS for the analysis."

### Statistical Testing
> "Run limma differential analysis comparing treatment vs control groups on my normalized protein matrix"

> "Use proDA for differential testing on my label-free data that has about 30% missing values"

### Complex Designs
> "Set up a limma model with treatment and batch as covariates in the design matrix"

> "Run paired differential analysis for my before/after samples"

### Results and Visualization
> "Create a volcano plot of the differential abundance results"

> "Filter to proteins with adjusted p-value below 0.05 and absolute log2FC above 1"

## What the Agent Will Do
1. Load raw protein intensity matrix and sample metadata
2. Preprocess: log2 transform raw intensities, apply normalization (median centering, cyclic loess, or VSN depending on data characteristics)
3. Select statistical method based on sample size, data type, and available metadata (limma for small n, DEqMS with PSM counts, proDA for extensive missing values, Welch's t-test for large n in Python)
4. Define experimental design and contrasts
5. Fit statistical model with empirical Bayes moderation
6. Apply Benjamini-Hochberg multiple testing correction
7. Choose fold change reporting strategy based on downstream use (raw for GSEA/meta-analysis, ashr-shrunk for effect size recovery)
8. Generate results table and optional visualizations

## Statistical Methods

| Method | Best for | Key advantage |
|--------|----------|---------------|
| limma | Small n (3-5), general purpose | Borrows variance across proteins; ~10-20 extra effective df |
| DEqMS | When PSM/peptide counts available | Weights variance by quantification depth |
| proDA | Label-free with >20% missing values | Models dropout without imputation |
| Welch's t-test | Large n (>10), Python-only | Simple, reliable with sufficient samples |
| MSstats | Complex designs, technical replicates | Feature-level mixed models |

## Normalization Methods

| Method | When to use | Key assumption |
|--------|-------------|----------------|
| Median centering | Default starting point; robust to missing values | Majority of proteins unchanged |
| Cyclic loess | Unbalanced DE (asymmetric up/down regulation) | Majority of proteins unchanged |
| VSN | Heteroscedastic data; input is raw (not log2) | Parametric variance model holds |
| Quantile | TMT with complete data | Identical sample distributions |

Median normalization subtracts each sample's median log2 value and optionally re-centers to the global median, ensuring all samples share a common center. This removes systematic loading differences while preserving biological signal.

## Fold Change Reporting

Raw fold changes are noisy but unbiased estimates of the true biological effect. How to handle them depends on what comes next:

- **GSEA / pathway analysis**: Use raw fold changes for all proteins. These methods rank by effect size and rely on the full continuous distribution, including small non-significant effects. Do not zero or threshold FCs before GSEA.
- **Effect size recovery** (e.g., "which proteins truly changed and by how much?"): Apply ashr shrinkage in R to produce posterior mean estimates. ashr smoothly shrinks uncertain effects toward zero while preserving well-supported ones. This is the principled Bayesian approach.
- **Reporting tables**: Report raw FC with adjusted p-value and confidence interval. The p-value communicates uncertainty; the FC communicates magnitude. Downstream consumers can threshold as needed.
- **Cross-study comparison / meta-analysis**: Use raw FCs with standard errors as input. Shrinkage is study-specific and should not be applied before pooling.

Avoid hard-thresholding FCs at a p-value cutoff (e.g., zeroing non-significant FCs). This creates artificial discontinuities and discards information that downstream methods may need.

## Significance Thresholds

Typical thresholds for proteomics:
- **Adjusted p-value**: < 0.05 (or 0.01 for stringent)
- **Log2 fold change**: > 1 (2-fold) or > 0.58 (1.5-fold)
- Use `treat()` + `topTreat()` in limma for minimum-effect-size testing rather than post-hoc FC filtering, which can inflate FDR

## Tips
- Always log2-transform raw intensities before normalization and testing
- Always use adjusted p-values (not raw) for significance calls
- Use `equal_var=False` in `scipy.stats.ttest_ind` for Welch's t-test (the default is Student's)
- Pass `method='fdr_bh'` explicitly to `statsmodels.stats.multitest.multipletests` (the default is Holm-Sidak, not BH)
- Include batch as a covariate in the design matrix if samples were processed separately; do not use `removeBatchEffect()` before testing
- Consider ashr fold change shrinkage in R when effect size accuracy matters; for GSEA or meta-analysis, use raw FCs
- Check volcano plot symmetry. Strongly asymmetric patterns may indicate normalization issues
- Report: number tested, normalization method, statistical method, thresholds, and number significant

## Related Skills

- quantification - Protein-level abundance estimation and normalization before testing
- proteomics-qc - Quality control and batch effect assessment
- differential-expression/deseq2-basics - Analogous empirical Bayes concepts for RNA-seq
- data-visualization/specialized-omics-plots - Volcano plots, MA plots, heatmaps
