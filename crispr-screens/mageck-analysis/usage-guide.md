# MAGeCK Analysis - Usage Guide

## Overview
MAGeCK is the standard tool for analyzing pooled CRISPR screens. It handles count normalization, identifies significantly enriched/depleted genes, and performs pathway analysis.

## Prerequisites
```bash
pip install mageck
# or via conda
conda install -c bioconda mageck
```

## Quick Start
Tell your AI agent what you want to do:
- "Run MAGeCK analysis on my CRISPR screen count data"
- "Compare treatment vs control samples to find essential genes"
- "Identify drug resistance genes from my screen"

## Example Prompts

### Basic Analysis
> "I have count data from a CRISPR dropout screen. Run MAGeCK test to compare treatment vs control and find essential genes."

> "Analyze my CRISPR screen with MAGeCK using median normalization. My samples are T1, T2 as treatment and C1, C2 as control."

### Multi-Condition Analysis
> "I have a time-course CRISPR screen with day 0, day 7, and day 14 samples. Use MAGeCK MLE to analyze all timepoints together."

> "Run MAGeCK MLE on my screen data accounting for batch effects between screening rounds."

### Interpretation
> "I have MAGeCK results. Help me interpret the RRA scores and identify the top depleted genes."

> "Compare my MAGeCK results between positive and negative selection. Which genes show the strongest effect?"

## What the Agent Will Do
1. Load count data and verify format
2. Select appropriate MAGeCK mode (test vs MLE)
3. Run normalization and statistical analysis
4. Generate gene-level rankings
5. Identify significantly enriched/depleted genes
6. Create visualization of results

## Tips
- Use `mageck test` for simple two-group comparisons
- Use `mageck mle` for complex designs with multiple conditions or covariates
- Median normalization is most robust for typical screens
- Quality thresholds: sgRNA mapping >70%, Gini <0.2, replicate correlation >0.8
- FDR < 0.1 for discovery, FDR < 0.05 for high confidence hits

## Key Concepts

### RRA Score
Robust Rank Aggregation score combining sgRNA-level statistics to gene level.

### Beta Score (MLE)
Effect size from maximum likelihood estimation, comparable across genes.

## References
- MAGeCK paper: doi:10.1186/s13059-014-0554-4
- MAGeCK-MLE: doi:10.1186/s13059-015-0843-6
- Documentation: https://sourceforge.net/p/mageck/wiki/Home/
