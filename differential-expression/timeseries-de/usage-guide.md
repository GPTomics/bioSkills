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

- Choose spline degrees based on number of timepoints (more timepoints = higher df possible)
- Include biological replicates at each timepoint for statistical power
- Test for group:time interaction to find condition-specific temporal dynamics
- Use ns() for natural splines or bs() for B-splines in design formulas
- maSigPro works well for experiments with multiple conditions and many timepoints
