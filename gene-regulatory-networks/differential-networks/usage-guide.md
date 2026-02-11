# Differential Networks - Usage Guide

## Overview

Compare co-expression networks between biological conditions to identify rewired gene-gene regulatory relationships. Uses Fisher's z-transform to test whether correlations differ significantly between groups, revealing gained, lost, and reversed co-expression edges. Complements differential expression analysis by showing how gene relationships change, not just expression levels.

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

- **Sample size** - Need at least 15 per group for stable correlations; 20+ is recommended
- **Gene filtering** - Restrict to top 2000-5000 variable genes to reduce multiple testing burden
- **Effect size** - Filter for meaningful differences (absolute correlation change > 0.3) in addition to statistical significance
- **DGCA archived** - DGCA was removed from CRAN in May 2024; install from GitHub via devtools if needed
- **Complement with DE** - Differential network analysis reveals regulatory rewiring that DE analysis misses; a gene can be rewired without changing expression level
- **Python alternative** - Use scipy and NetworkX for a pure Python approach if R is not available

## Related Skills

- coexpression-networks - Build WGCNA networks for individual conditions first
- scenic-regulons - TF regulon inference from scRNA-seq with pySCENIC
- differential-expression/deseq2-basics - DE analysis to complement network rewiring
- pathway-analysis/go-enrichment - Functional enrichment of rewired gene sets
- temporal-genomics/temporal-grn - Time-delayed regulatory inference from temporal data
