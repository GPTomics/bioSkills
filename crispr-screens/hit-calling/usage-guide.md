# Hit Calling - Usage Guide

## Overview
Multiple methods exist for calling hits in CRISPR screens. The choice depends on screen design, reference data availability, and desired stringency.

## Prerequisites
```bash
# MAGeCK
pip install mageck
# BAGEL2
pip install bagel
# Or via conda
conda install -c bioconda mageck bagel2
```

## Quick Start
Tell your AI agent what you want to do:
- "Call significant hits from my CRISPR screen results"
- "Identify essential genes using multiple methods"
- "Find consensus hits across MAGeCK and BAGEL2"

## Example Prompts

### Single Method Analysis
> "Run MAGeCK RRA on my screen and call hits at FDR < 0.1."

> "Use BAGEL2 to identify essential genes. I have reference essential and non-essential gene lists."

> "Apply drugZ analysis to my drug screen data to find resistance genes."

### Threshold Selection
> "What's the appropriate FDR threshold for my exploratory CRISPR screen?"

> "I need high-confidence hits. What BAGEL2 Bayes Factor cutoff should I use?"

### Consensus Approach
> "Run MAGeCK, BAGEL2, and drugZ on my data and find genes called by at least 2 methods."

> "Compare hit lists between different analysis methods. Which genes are consistent?"

### Selection Direction
> "Find genes that drop out (negative selection) in my essentiality screen."

> "Identify resistance genes (positive selection) from my drug treatment screen."

## What the Agent Will Do
1. Select appropriate hit-calling method(s)
2. Run statistical analysis
3. Apply significance thresholds
4. Separate positive and negative selection hits
5. Generate ranked gene lists
6. Optionally find consensus across methods

## Tips
- Use multiple methods and take consensus for robust hits
- FDR < 0.1 for discovery, FDR < 0.05 for high confidence
- BAGEL2 Bayes Factor > 5 indicates strong evidence
- Z-score |Z| > 3 corresponds to approximately p < 0.003
- Consider both positive and negative selection directions

## Method Comparison

| Method | Approach | Best For |
|--------|----------|----------|
| MAGeCK RRA | Rank aggregation | General screens, no training data |
| MAGeCK MLE | Maximum likelihood | Complex designs, multiple conditions |
| BAGEL2 | Bayesian | Well-characterized systems with reference genes |
| drugZ | Z-score | Drug screens, simple interpretation |

## Significance Thresholds

| Method | Suggestive | Strong | Very Strong |
|--------|------------|--------|-------------|
| MAGeCK | FDR < 0.1 | FDR < 0.05 | FDR < 0.01 |
| BAGEL2 | BF > 3 | BF > 5 | BF > 10 |
| Z-score | |Z| > 2 | |Z| > 3 | |Z| > 4 |

## References
- BAGEL2: doi:10.1186/s13059-019-1749-z
- drugZ: doi:10.1186/s13073-019-0665-3
