# Differential Networks - Usage Guide

## Overview

Compare co-expression and regulatory networks between conditions to find rewired gene-gene relationships. DiffCorr tests correlation differences with Fisher's z (gained/lost/reversed edges); DiffCoEx finds rewired modules; DINGO/iDINGO test direct (partial-correlation) rewiring; CoDiNA compares many networks. The key insight: differential connectivity is NOT differential expression -- a gene can rewire its partners with no change in mean expression and be the key signal (Hudson's myostatin). The dominant pitfall is power: rewiring needs more samples than DE, and pairwise testing has a p^2/2 multiple-testing explosion, so most small-cohort "rewired hub" claims are underpowered noise.

## Prerequisites

```r
# R packages
install.packages('DiffCorr')
install.packages('igraph')

# DGCA (archived from CRAN 2024-05, GitHub-only)
devtools::install_github('andymckenzie/DGCA')
```

```bash
# Python alternative
pip install pandas numpy scipy statsmodels networkx matplotlib
```

Minimum 20 samples per condition recommended (15 absolute floor).

## Quick Start

Tell your AI agent what you want to do:
- "Compare co-expression networks between disease and control"
- "Find genes with rewired regulatory connections between conditions"
- "Which gene pairs gain or lose co-expression in my treatment group?"
- "Identify hub genes with the most differential connections"

## Example Prompts

### DiffCorr Analysis
> "I have normalized RNA-seq data from 25 disease and 25 control samples. Find differentially correlated gene pairs."

> "Run DiffCorr to compare co-expression between treated and untreated groups."

### Network Comparison
> "Which genes have the most rewired connections between my tumor and normal samples?"

> "Show me gained, lost, and reversed edges in the differential network."

### Direct vs Module-Level Rewiring
> "I want direct (not indirect) rewiring between conditions. Run DINGO on partial correlations."

> "Find modules that rewire between conditions with DiffCoEx instead of testing every edge."

### Visualization
> "Visualize the differential co-expression network colored by edge type."

> "Create a network plot showing rewired connections for the top 50 most differentially connected genes."

## What the Agent Will Do

1. Separate expression data by condition
2. Filter to top variable genes (2000-5000)
3. Compute correlation matrices for each condition
4. Test differential correlation using Fisher's z-transform
5. Classify edges as gained, lost, reversed, or unchanged
6. Apply FDR correction and identify rewired hub genes
7. Visualize differential network

## Tips

- Connectivity is not expression - report differential connectivity and DE separately; a non-DE gene can be the key rewired hub (Hudson's myostatin).
- Power is the binding constraint - rewiring needs more samples than DE; below ~15-20 per group, results are mostly noise. Treat small-cohort findings as exploratory.
- Tame the p^2/2 explosion - pre-filter to variable genes and apply strict FDR, or use module-level DiffCoEx to avoid per-edge testing entirely.
- Marginal vs direct - DiffCorr gained edges are marginal (may reflect a shifted common driver); use DINGO partial correlations when directness matters.
- FDR method - statsmodels `multipletests` defaults to Holm-Sidak; pass `method='fdr_bh'` for BH.
- Effect-size filter - require absolute correlation change > 0.3 on top of significance to avoid trivial differences.
- DGCA archived - removed from CRAN May 2024; install from GitHub if needed.

## Related Skills

- coexpression-networks - build the per-condition networks being compared
- grn-inference - VIPER differential protein activity between conditions
- scenic-regulons - TF regulon activity differences as a complementary readout
- differential-expression/de-results - differential expression of means (the orthogonal question)
- pathway-analysis/go-enrichment - functional enrichment of rewired gene sets
- temporal-genomics/temporal-grn - time-resolved network change across stages
