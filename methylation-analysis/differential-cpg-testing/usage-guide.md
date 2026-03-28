# Differential CpG Testing - Usage Guide

## Overview

Per-CpG differential methylation testing from bisulfite sequencing count data or beta-value matrices. Covers the full pipeline from raw counts through beta/M-value computation, coverage filtering, statistical testing, multiple testing correction, and effect size reporting.

## Prerequisites

```bash
pip install numpy pandas scipy statsmodels
```

```r
BiocManager::install(c('limma', 'DSS'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Test each CpG for differential methylation between my case and control groups"
- "Run Welch's t-test on beta values from my bisulfite count data with BH FDR correction"
- "Use limma on M-values to find differentially methylated CpGs in my small-sample experiment"
- "Compute delta-beta effect sizes and identify significant CpGs"

## Example Prompts

### From Count Data
> "I have a TSV with per-CpG methylated and total read counts for 6 case and 6 control samples. Compute beta values, filter CpGs with any sample below 10x coverage, run Welch's t-test, and apply BH FDR correction."

> "Read my bisulfite count data, compute beta values as meth/total, and test each CpG for differential methylation between treatment groups."

### Method Selection
> "I have only 3 samples per group from RRBS. Use limma on M-values for per-CpG differential methylation testing."

> "My beta value distributions look non-normal at many CpGs. Use Mann-Whitney U test instead of Welch's t-test."

> "Analyze my BS-seq count data with DSS beta-binomial model for statistically rigorous per-CpG testing."

### Results and Effect Sizes
> "Filter my differential methylation results to CpGs with |delta_beta| > 0.2 and FDR < 0.05."

> "Generate a results table with cpg_id, mean_case_beta, mean_ctrl_beta, delta_beta, pvalue, padj, and significant columns."

## What the Agent Will Do

1. Read per-CpG methylation data (count matrix or beta-value matrix)
2. Compute beta values from counts if starting from raw data (beta = meth / total)
3. Apply coverage filtering (minimum 10x default, 99.9th percentile upper limit)
4. Select statistical test based on sample size and data characteristics
5. Run per-CpG testing between groups (Welch's t-test, Mann-Whitney, limma, or DSS)
6. Apply Benjamini-Hochberg FDR correction for multiple testing
7. Calculate delta_beta effect sizes (mean case beta minus mean control beta)
8. Output results table with statistics, adjusted p-values, and significance calls

## Method Selection Guide

| Scenario | Recommended | Key Advantage |
|----------|-------------|---------------|
| Large n (>10/group), Python | Welch's t-test + BH | Fast, simple, reliable with adequate samples |
| Large n, non-normal | Mann-Whitney U + BH | No distributional assumptions |
| Small n (3-5/group) | limma on M-values | Borrows variance across CpGs via empirical Bayes |
| Count-level modeling | DSS beta-binomial | Models biological + sampling variance jointly |
| Unreplicated | Fisher's exact | Only option without replicates (use with caution) |

## Tips

- For small samples (n < 5/group), use M-values (logit of beta) for statistical testing with limma -- M-values are homoscedastic unlike beta values. For large samples (n > 10/group), t-test on beta values directly is standard and acceptable
- Report delta_beta (difference in mean beta values) for effect sizes -- directly interpretable as percentage-point methylation difference
- Always pass `equal_var=False` to `scipy.stats.ttest_ind` -- the default is Student's t-test which assumes equal variance
- Always pass `method='fdr_bh'` to `statsmodels.stats.multitest.multipletests` -- the default is Holm-Sidak, not Benjamini-Hochberg
- For fewer than 5 samples per group, strongly prefer limma on M-values over a simple t-test
- Coverage filtering at 10x is standard; increase to 30x+ for targeted bisulfite sequencing
- Effect size thresholds: 0.10 for discovery/EWAS, 0.20 for standard analysis, 0.30 for stringent filtering
- Compute delta_beta from original beta values, not from model coefficients -- methylKit's meth.diff and limma's logFC on M-values do not equal mean(case) - mean(ctrl) on beta values
- Do not pool counts across biological replicates for Fisher's exact test -- this ignores biological variance and inflates false positives

## Related Skills

- methylkit-analysis - Chi-squared testing with methylKit in R
- methylation-calling - Generate per-CpG count data from Bismark BAM files
- dmr-detection - Region-level testing after per-CpG analysis
- differential-expression/deseq2-basics - Analogous empirical Bayes concepts for count data
- proteomics/differential-abundance - Similar Welch t-test + BH workflow for proteomics
