# Statistical Analysis - Usage Guide

## Overview

Statistical analysis identifies metabolites associated with biological conditions. Covers the full pipeline from raw intensities through preprocessing (log2 transformation, normalization), statistical modeling with limma or Welch's t-tests, fold change estimation, and multivariate methods for biomarker discovery.

## Prerequisites

```bash
# Python
pip install scipy statsmodels numpy pandas matplotlib scikit-learn
```

```r
# R packages
BiocManager::install(c("limma"))
install.packages(c("mixOmics", "ropls", "pROC", "ashr"))
```

## Quick Start

Tell your AI agent what you want to do:
- "Find differentially abundant metabolites between treatment and control"
- "Log2 transform and normalize my metabolomics data, then run limma for differential analysis"
- "Identify significant metabolites with shrunk fold change estimates"
- "Run PLS-DA and identify important features by VIP score"

## Example Prompts

### Full Pipeline
> "I have raw metabolite intensities in a TSV file. Log2 transform, PQN normalize, and run differential testing between case and control groups. Report fold changes and adjusted p-values."

> "Process my untargeted metabolomics feature table: handle missing values, normalize, run limma with empirical Bayes, and generate a volcano plot."

### Univariate Analysis
> "Run limma moderated t-tests on my log2-transformed metabolomics data with BH FDR correction"

> "Perform Welch's t-tests in Python on my metabolomics data with Benjamini-Hochberg correction"

> "Perform ANOVA across my three treatment groups and identify significant metabolites"

### Multivariate Analysis
> "Run PCA for exploratory analysis and check sample grouping"

> "Build a PLS-DA model with 10-fold cross-validation and calculate VIP scores"

> "Use OPLS-DA for biomarker discovery between disease and healthy groups"

### Biomarker Selection
> "Identify metabolites with VIP > 1 and FDR < 0.05"

> "Calculate ROC curves and AUC for top candidate biomarkers"

### Fold Change and Effect Size
> "Apply ashr fold change shrinkage to get more accurate effect size estimates from my limma results"

> "Use treat() to test for metabolites with at least 1.5-fold change"

## What the Agent Will Do

1. Load raw intensity matrix and sample metadata
2. Preprocess: handle zeros (distinguish technical vs biological), log2 transform, apply normalization (PQN for untargeted LC-MS, or skip if data is pre-normalized)
3. Select statistical method based on sample size and design (limma for small n, Welch's for large n in Python, Wilcoxon for non-normal data)
4. Define experimental design and contrasts
5. Fit statistical model with empirical Bayes moderation (limma) or per-feature testing (Python)
6. Apply Benjamini-Hochberg multiple testing correction
7. Choose fold change reporting strategy (raw for GSEA/pathway analysis, ashr-shrunk for effect size recovery)
8. Generate results table and optional visualizations (volcano plot, PCA, heatmap)
9. Optionally run multivariate models (PLS-DA, OPLS-DA) for classification and VIP ranking

## Statistical Methods

| Method | Best for | Key advantage |
|--------|----------|---------------|
| limma | Small n (3-10), general purpose | Borrows variance across features; ~10-20 extra effective df |
| Welch's t-test | Large n (>10), Python-only | Simple, reliable with sufficient samples |
| Wilcoxon rank-sum | Non-normal after log transform | Distribution-free; less powerful when normality holds |
| Two-part test | Zero-inflated features (BPMVs) | Separately tests presence/absence and abundance |
| PCA | Any | Unsupervised exploration, QC, outlier detection |
| PLS-DA | 2+ groups | Supervised classification with VIP ranking |
| OPLS-DA | 2 groups | Separates predictive from orthogonal variation |

## Normalization Methods

| Method | When to use | Key assumption |
|--------|-------------|----------------|
| PQN | Default for untargeted LC-MS metabolomics | Majority of features unchanged |
| QC-RSC (LOESS) | Multi-batch studies with pooled QC samples | QC samples represent technical drift |
| VSN | High zero rate; heteroscedastic data | Parametric variance model holds |
| TIC | NMR metabolomics; quick exploration | All features contribute equally |
| Cyclic loess | Asymmetric DE pattern | Majority of features unchanged |
| None | IS-corrected; single-batch balanced design | No systematic bias present |

PQN is preferred over median centering or TIC for untargeted metabolomics because it is more robust to dominant high-abundance features.

## Fold Change Reporting

Raw fold changes are noisy but unbiased. How to handle them depends on what comes next:

- **Pathway analysis (QEA/GSEA)**: Use raw fold changes for all features. Do not zero or threshold FCs before enrichment analysis.
- **Effect size recovery**: Apply ashr shrinkage in R for posterior mean estimates. Smoothly shrinks uncertain effects toward zero.
- **Reporting tables**: Report raw FC with adjusted p-value. The p-value communicates uncertainty; the FC communicates magnitude.
- **Cross-study / meta-analysis**: Use raw FCs with standard errors. Shrinkage is study-specific.

Avoid hard-thresholding FCs at a p-value cutoff (setting non-significant FCs to zero). Use `treat()` + `topTreat()` in limma for minimum-effect-size testing instead of post-hoc FC filtering.

## Significance Thresholds

| Metric | Threshold |
|--------|-----------|
| Adjusted p-value | < 0.05 (discovery: < 0.1) |
| Log2 fold change | > 1 (2-fold) or > 0.58 (1.5-fold) |
| VIP score | > 1 |
| AUC | > 0.7 (moderate), > 0.8 (good) |

## Tips

- Always log2-transform raw intensities before statistical testing; compute fold change as difference of means on the log2 scale (geometric mean ratio), not `log2(mean_ratio)` which uses arithmetic means
- Use limma with `eBayes(trend=TRUE, robust=TRUE)` as the default for metabolomics -- even with n > 5, it is never worse than a plain t-test and is often better
- In Python, set `equal_var=False` in `scipy.stats.ttest_ind()` for Welch's t-test (scipy defaults to Student's)
- Pass `method='fdr_bh'` explicitly to `statsmodels.multipletests()` (the default is Holm-Sidak, not BH)
- Include batch as a covariate in the design matrix if samples were processed separately; do not use `removeBatchEffect()` before testing
- MetaboAnalyst uses arithmetic mean FC and Student's t-test by default; be aware when comparing results
- Consider ashr fold change shrinkage in R when effect size accuracy matters; for GSEA or pathway analysis, use raw FCs
- Distinguish biological zeros (metabolite truly absent) from technical zeros (below detection limit) before imputation
- Validate PLS-DA with permutation testing (Q2 should exceed permuted values)
- VIP > 1 is a common threshold, but combine with FDR for confidence
- Check volcano plot symmetry -- strongly asymmetric patterns may indicate normalization issues
- Report: number tested, normalization method, statistical method, thresholds, and number significant

## Related Skills

- normalization-qc - Data preparation and batch correction
- pathway-mapping - Functional interpretation of differential metabolites
- multi-omics-integration/mixomics-analysis - Advanced multivariate methods
- proteomics/differential-abundance - Analogous empirical Bayes concepts for proteomics
- data-visualization/specialized-omics-plots - Volcano plots, MA plots, heatmaps

## References

- limma: doi:10.1093/nar/gkv007
- ashr: doi:10.1093/biostatistics/kxw041
- mixOmics: doi:10.1371/journal.pcbi.1005752
- ropls (OPLS-DA): doi:10.1021/acs.jproteome.5b00354
- PQN normalization: doi:10.1021/ac051632c
