# Differential Analysis - Usage Guide

## Overview
Differential analysis identifies cell populations that differ in abundance or marker expression between experimental conditions using rigorous statistical methods.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('diffcyt', 'CATALYST', 'edgeR', 'limma'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Compare cell population frequencies between treatment groups"
- "Test for differential marker expression within clusters"
- "Run diffcyt analysis on my CyTOF experiment"

## Example Prompts
### Differential Abundance
> "Test which clusters are more abundant in treatment vs control"
> "Run edgeR-based differential abundance analysis"
> "Show a volcano plot of cluster abundance changes"

### Differential State
> "Test for differential marker expression within each cluster"
> "Compare Ki67 expression between conditions across all populations"
> "Run limma-based differential state analysis"

### Experimental Design
> "Set up a paired differential analysis for my pre/post treatment samples"
> "Account for batch effects in my differential abundance model"
> "Run differential analysis with a continuous covariate"

## What the Agent Will Do
1. Create experiment design matrix from sample metadata
2. Set up contrast for conditions of interest
3. Run differential abundance (DA) analysis with edgeR
4. Run differential state (DS) analysis with limma
5. Apply multiple testing correction and summarize results

## Tips
- Need biological replicates (n >= 3 per group) for valid statistics
- DA uses cluster counts/proportions; DS uses marker expression
- diffcyt wraps edgeR (DA) and limma (DS) in a cytometry workflow
- Correct for number of clusters (DA) or clusters x markers (DS)
- Consider paired designs when samples are matched

## Analysis Types

| Type | Question | Method | Output |
|------|----------|--------|--------|
| DA | Are proportions different? | edgeR | Fold change in frequency |
| DS | Is expression different? | limma | Marker changes per cluster |

## Interpretation

| Metric | Meaning |
|--------|---------|
| logFC > 0 | Higher in treatment group |
| logFC < 0 | Higher in control group |
| p_adj < 0.05 | Statistically significant |

## References
- diffcyt: doi:10.1038/s41467-017-00707-4
- CITRUS: doi:10.1073/pnas.1408792111
