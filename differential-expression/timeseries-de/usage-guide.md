# Time-Series DE - Usage Guide

## Overview

Identify genes with significant temporal expression patterns across time-course experiments using spline models, polynomial regression, or likelihood ratio tests.

## Prerequisites

```r
BiocManager::install(c('limma', 'edgeR', 'maSigPro', 'DESeq2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Find genes with significant temporal patterns in my time-course data"
- "Identify genes that respond differently over time between conditions"
- "Cluster genes by their temporal expression profiles"

## Example Prompts

### limma with Splines
> "Fit a spline model to my time-course RNA-seq data"

> "Test for genes with significant time effects using limma"

> "Find genes with condition-specific temporal dynamics"

### maSigPro Analysis
> "Run maSigPro on my multi-condition time-series experiment"

> "Identify genes with polynomial time patterns"

> "Cluster significant genes by temporal profile"

### DESeq2 Approach
> "Use DESeq2 LRT to test for time effects"

> "Compare full model with time against reduced model"

### Visualization
> "Plot expression trajectories for my top time-varying genes"

> "Create a heatmap of genes clustered by temporal pattern"

## What the Agent Will Do

1. Set up appropriate design matrix with time variable (splines or polynomial)
2. Normalize counts and apply voom transformation
3. Fit linear model with time terms
4. Test for significant time effects or time:condition interactions
5. Cluster significant genes by temporal profile
6. Visualize expression trajectories

## Method Comparison

| Method | Best For |
|--------|----------|
| limma + splines | Smooth temporal patterns |
| maSigPro | Multiple conditions over time |
| ImpulseDE2 | Impulse-like responses |
| DESeq2 LRT | Discrete time comparisons |

## Tips

- Choose spline degrees based on number of timepoints (rule of thumb: df <= unique time points / 2)
- Include biological replicates at each timepoint for statistical power
- Test for group:time interaction to find condition-specific temporal dynamics
- Use ns() for natural splines or bs() for B-splines in design formulas
- maSigPro works well for experiments with multiple conditions and many timepoints; cite Nueda 2014 (RNA-seq update), NOT the 2006 Conesa paper (microarray original)
- For short series (<8 time points), naive pairwise DESeq2 + LRT typically OUTPERFORMS dedicated time-course tools per the Spies 2019 benchmark; reserve ImpulseDE2 for longer series with monotonic-then-saturating biology
- For repeated measures (same subject sampled over time), use DREAM (linear mixed model with `(1|subject)`) -- treating repeated observations as independent is pseudoreplication and inflates type-I error
- ImpulseDE2 was removed from Bioconductor at the 3.13 release (May 2021); install from BiocArchive or the YosefLab GitHub mirror

## Related Skills

- differential-expression/deseq2-basics - Standard DE analysis
- differential-expression/de-visualization - Visualize results
- differential-expression/batch-correction - Handle batch effects
- pathway-analysis/go-enrichment - Functional analysis of clusters
- temporal-genomics/circadian-rhythms - Circadian rhythm detection for time-course data
- temporal-genomics/temporal-clustering - Cluster genes by temporal expression profile
- temporal-genomics/trajectory-modeling - GAM trajectory fitting for temporal expression data
- temporal-genomics/temporal-grn - Dynamic GRN inference from bulk time-series data
