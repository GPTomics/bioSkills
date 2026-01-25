# Power Analysis Usage Guide

## Overview

This guide covers calculating statistical power for RNA-seq, ATAC-seq, and other sequencing experiments.

## Prerequisites

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('RNASeqPower', 'ssizeRNA'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Calculate power for my RNA-seq experiment with 3 replicates per group"
- "How many samples do I need to detect 2-fold changes in RNA-seq?"
- "What's the minimum effect size I can detect with my current sample size?"
- "Run a power analysis for my ATAC-seq experiment"

## Example Prompts

### RNA-seq Power

> "I'm planning an RNA-seq experiment comparing drug-treated vs control cells. I expect about 0.3 CV and want to detect 1.5-fold changes. How many replicates do I need for 80% power?"

> "Calculate the power of my RNA-seq study with 4 samples per group, 20x coverage, assuming CV of 0.4"

### ATAC-seq Power

> "Run power analysis for my ATAC-seq comparing two conditions with 3 replicates each"

### Budget Optimization

> "I have budget for 10 samples total. Should I do 5 vs 5 or 4 vs 6 for my RNA-seq experiment?"

## What the Agent Will Do

1. Identify experiment type (RNA-seq, ATAC-seq, etc.)
2. Determine key parameters (CV, effect size, depth)
3. Use appropriate power calculation package
4. Report power or required sample size
5. Provide interpretation and recommendations

## Tips

- Estimate CV from pilot data or literature for your tissue type
- For human samples, assume higher CV (0.3-0.5) than cell lines (0.1-0.2)
- Power > 0.8 is standard; > 0.9 for critical studies
- More biological replicates beats deeper sequencing after ~20M reads
- Consider practical constraints (cost, sample availability) alongside statistics
