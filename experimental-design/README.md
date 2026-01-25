# experimental-design

## Overview

Skills for planning sequencing experiments, calculating power and sample sizes, and designing studies to minimize batch effects.

**Tool type:** r | **Primary tools:** RNASeqPower, ssizeRNA, qvalue, sva

## Skills

| Skill | Description |
|-------|-------------|
| power-analysis | Statistical power calculations for RNA-seq, ATAC-seq experiments |
| sample-size | Sample size estimation for differential expression studies |
| multiple-testing | FDR, Bonferroni, and q-value correction methods |
| batch-design | Experimental design to minimize and correct batch effects |

## Example Prompts

- "How many samples do I need for my RNA-seq experiment to detect 2-fold changes?"
- "Calculate power for my ATAC-seq study with 4 replicates per group"
- "Help me assign 24 samples to 3 sequencing batches without confounding"
- "Which multiple testing correction should I use for my differential expression results?"
- "What's the minimum effect size I can detect with 6 samples per group?"

## Requirements

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('RNASeqPower', 'ssizeRNA', 'qvalue', 'sva', 'limma'))

# Optional
install.packages('designit')
```

## Related Skills

- **differential-expression** - Run differential expression after proper design
- **single-cell** - scRNA-seq experimental design considerations
- **read-qc** - Quality control to validate experimental success
