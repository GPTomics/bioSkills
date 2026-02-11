# Biodiversity Metrics Usage Guide

## Overview

Quantifies biodiversity using the Hill number framework for alpha diversity (iNEXT coverage-based rarefaction/extrapolation), classic diversity indices (vegan), and beta diversity partitioning into turnover and nestedness components (betapart). Coverage-based standardization enables fair comparison of assemblages sampled with different effort, and the Hill number family unifies richness (q=0), Shannon (q=1), and Simpson (q=2) diversity into a single framework.

## Prerequisites

- R with iNEXT (`install.packages('iNEXT')`)
- iNEXT.3D for phylogenetic and functional diversity (`install.packages('iNEXT.3D')`)
- vegan for classic diversity metrics (`install.packages('vegan')`)
- betapart for beta diversity decomposition (`install.packages('betapart')`)
- ggplot2 for visualization (`install.packages('ggplot2')`)

## Quick Start

Tell your AI agent what you want to do:

- "Compare species diversity across my sampling sites using Hill number rarefaction"
- "Build rarefaction/extrapolation curves for my eDNA species abundance data"
- "Decompose beta diversity into turnover and nestedness components"
- "Estimate asymptotic species richness for each site using Chao1"
- "Calculate coverage-standardized diversity at 95% coverage across sites"

## Example Prompts

### Alpha Diversity

> "I have species abundance data from 12 sites. Calculate Hill numbers (q=0, 1, 2) using iNEXT with coverage-based rarefaction and plot the curves."

> "Estimate the true species richness at each site using Chao1 asymptotic estimation, and tell me how complete my sampling is."

### Beta Diversity

> "Decompose pairwise beta diversity into turnover and nestedness using betapart Sorensen indices. Which component dominates?"

> "Calculate multi-site beta diversity for my sampling region and determine whether community differences are driven by species replacement or richness gradients."

### Multidimensional Diversity

> "I have a species abundance matrix, a phylogenetic tree, and a trait distance matrix. Calculate taxonomic, phylogenetic, and functional diversity using iNEXT.3D."

### Comparison

> "Compare Shannon diversity between two habitat types and test whether the difference is statistically significant."

## What the Agent Will Do

1. Load the species-by-site abundance or incidence matrix
2. Compute Hill numbers at orders q = 0, 1, and 2 using iNEXT with coverage-based rarefaction/extrapolation
3. Generate rarefaction curves, sample completeness profiles, and coverage-based diversity comparisons
4. Calculate asymptotic diversity estimates and sample coverage
5. Decompose beta diversity into turnover and nestedness components using betapart
6. Produce publication-ready plots of diversity curves and beta diversity patterns
7. Report whether assemblage differences are driven by species replacement or nested subsets

## Tips

- Coverage-based rarefaction (iNEXT type 3) is preferred over sample-size rarefaction because it accounts for sampling completeness rather than arbitrary read counts
- The Hill number at q = 1 (Shannon exponential) is the most balanced measure, weighting species proportionally to their frequency
- For eDNA data, use `datatype = 'incidence_freq'` if working with detection/non-detection across replicates rather than read counts
- betapart requires presence/absence data; convert abundance matrices with `ifelse(x > 0, 1, 0)` before use
- When turnover dominates beta diversity, communities replace species along gradients; when nestedness dominates, depauperate sites are subsets of richer ones
- For phylogenetic diversity with iNEXT.3D, the tree must be ultrametric; use `ape::chronos()` if needed
- Bootstrap confidence intervals (nboot = 200) are essential for determining whether diversity differences are significant
- At least 1000 reads per sample are needed for reliable rarefaction; below this, extrapolation becomes unreliable

## Related Skills

- edna-metabarcoding - Generate species tables from eDNA data
- community-ecology - Constrained ordination of assemblages
- microbiome/diversity-analysis - 16S microbiome diversity metrics
- data-visualization/ggplot2-fundamentals - Customize diversity plots
