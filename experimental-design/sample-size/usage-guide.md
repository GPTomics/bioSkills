# Sample Size Estimation Usage Guide

## Overview

This guide covers estimating required sample sizes for differential expression, ChIP-seq, methylation, and proteomics studies.

## Prerequisites

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('ssizeRNA', 'DESeq2', 'powsimR'))
```

## Quick Start

Tell your AI agent what you want to do:
- "How many samples do I need for my RNA-seq experiment?"
- "Estimate sample size for detecting 2-fold changes with 80% power"
- "What's the minimum replicates for a ChIP-seq experiment?"
- "Calculate sample size for my single-cell RNA-seq study"

## Example Prompts

### Bulk RNA-seq

> "I'm planning a bulk RNA-seq study comparing tumor vs normal. I expect to find about 500 DE genes out of 15,000. How many samples per group do I need?"

> "Estimate sample size for RNA-seq with expected fold change of 1.5 and dispersion around 0.2"

### Single-cell

> "How many cells and samples do I need for a scRNA-seq experiment comparing two conditions?"

### Multi-omic Studies

> "I'm planning a multi-omic study with RNA-seq and proteomics. What sample sizes do I need for each?"

## What the Agent Will Do

1. Identify assay type and study design
2. Gather required parameters (effect size, power, FDR)
3. Apply appropriate sample size formula
4. Account for expected dropouts or failed samples
5. Provide practical recommendations

## Tips

- Always add 10-20% extra samples for potential failures
- Biological replicates are more valuable than technical replicates
- For scRNA-seq, both cell count and sample count matter
- Consider pilot studies for parameter estimation
- Document assumptions for grant applications
