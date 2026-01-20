# CRISPR Screen QC - Usage Guide

## Overview
Quality control is critical for CRISPR screens. Poor library representation or technical issues can lead to false positives/negatives and unreliable results.

## Prerequisites
```bash
pip install pandas numpy scipy matplotlib seaborn
# MAGeCK for count processing
pip install mageck
```

## Quick Start
Tell your AI agent what you want to do:
- "Run QC on my CRISPR screen count data"
- "Check library representation and replicate correlation"
- "Validate my screen by checking essential gene dropout"

## Example Prompts

### Library Quality
> "Calculate the Gini index and zero-count percentage for my CRISPR screen to check library representation."

> "How many sgRNAs have zero counts in my screen? Is my library representation acceptable?"

### Replicate Correlation
> "Calculate Pearson correlation between my screen replicates. Are they consistent enough to proceed?"

> "Generate a correlation heatmap for all my screen samples to identify outliers."

### Essential Gene Validation
> "Check if known essential genes are depleted in my dropout screen. Calculate the AUC against DepMap essentials."

> "Validate my screen worked by checking the separation between essential and non-essential genes."

### Comprehensive QC
> "Run comprehensive QC on my CRISPR screen: check representation, calculate Gini, correlation, and essential gene recovery."

> "I'm worried about my screen quality. Run all standard QC metrics and tell me if there are problems."

## What the Agent Will Do
1. Load count data and calculate basic statistics
2. Compute Gini index for distribution evenness
3. Count zero and low-count sgRNAs
4. Calculate replicate correlations
5. Test essential gene dropout (if references provided)
6. Generate QC plots and summary report

## Tips
- Target Gini index <0.2 (ideal) or <0.3 (acceptable)
- Zero-count sgRNAs should be <1% (ideal) or <5% (acceptable)
- Replicate Pearson r >0.8 is ideal, >0.6 minimum
- Essential gene AUC >0.85 validates screen performance
- Always use log-transformed counts for correlation

## Key QC Metrics

| Metric | Ideal | Acceptable | Action if Failing |
|--------|-------|------------|-------------------|
| Zero-count sgRNAs | <1% | <5% | Increase cell coverage |
| Gini index | <0.2 | <0.3 | Reduce PCR cycles |
| Replicate correlation | >0.8 | >0.6 | Process replicates together |
| Essential gene AUC | >0.85 | >0.75 | Extend selection time |

## Common Issues and Fixes
- **Low representation**: Increase cell number, optimize PCR
- **Skewed distribution**: Reduce PCR cycles, new library prep
- **Poor correlation**: Process replicates together, check batch effects
- **No essential dropout**: Extend selection, verify cell model

## References
- Hart et al. essentials: doi:10.1016/j.molcel.2017.06.014
- DepMap: https://depmap.org/
