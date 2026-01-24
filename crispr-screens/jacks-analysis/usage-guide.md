# JACKS CRISPR Screen Analysis

## Overview

JACKS (Joint Analysis of CRISPR/Cas9 Knockout Screens) models both gene essentiality and sgRNA efficacy simultaneously. It's particularly powerful for multi-experiment analyses where sgRNA performance can be learned across screens.

## Prerequisites

```bash
pip install jacks
# or from source:
git clone https://github.com/felicityallen/JACKS.git
cd JACKS && pip install -e .
```

## Quick Start

Tell your AI agent what you want to do:
- "Run JACKS on my CRISPR screen count data"
- "Analyze multiple screens jointly with JACKS"
- "Identify essential genes and estimate sgRNA efficacy"
- "Compare JACKS results with MAGeCK"

## Example Prompts

### Basic Analysis
> "Run JACKS analysis on my CRISPR screen counts"

> "Identify essential genes from my dropout screen using JACKS"

> "Estimate sgRNA efficacy across my experiments"

### Multi-Screen Analysis
> "Jointly analyze my three CRISPR screens with JACKS"

> "Learn shared sgRNA efficacy across experiments"

> "Compare gene essentiality between different cell lines"

### Interpretation
> "Which genes are significantly essential in my JACKS results?"

> "Identify poor-performing sgRNAs from the efficacy scores"

> "Create a volcano plot of JACKS gene scores"

### Comparison
> "Compare JACKS and MAGeCK results for my screen"

> "Correlate gene rankings between the two methods"

## What the Agent Will Do

1. Prepare input files (counts, replicate map, guide-gene map)
2. Run JACKS inference with appropriate parameters
3. Extract gene-level essentiality scores with significance
4. Analyze sgRNA efficacy distributions
5. Generate visualizations (volcano plot, efficacy histogram)
6. Identify significant hits and poor-performing guides

## Tips

- JACKS is best when analyzing multiple screens simultaneously
- Use 50,000+ iterations for publication-quality results; 10,000 for exploration
- sgRNA efficacy scores help identify guides to exclude from future libraries
- Negative JACKS scores indicate essential genes (dropout phenotype)
- FDR log10 < -1 corresponds to FDR < 0.1; use -2 for FDR < 0.01
- Compare with MAGeCK for validation; high correlation expected for robust hits
- JACKS is slower than MAGeCK but provides additional sgRNA-level insights
