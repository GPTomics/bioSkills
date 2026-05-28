# Power Analysis Usage Guide

## Overview

This guide covers calculating statistical power for high-dimensional genomics experiments: bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, and proteomics. The central idea is that power in a sequencing experiment is a per-gene quantity that depends on each gene's mean expression and dispersion, so the honest summary is the marginal (average) power across the transcriptome at a target false discovery rate. Closed-form formulas using a single coefficient of variation give a fast approximation; simulation from the empirical mean-dispersion trend is the defensible default for count data. The guide also covers the depth-versus-replicate tradeoff and why observed (post-hoc) power should never be used to interpret a null result.

## Prerequisites

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('RNASeqPower', 'PROPER', 'DESeq2', 'edgeR'))
install.packages('pwr')

# powsimR (GitHub; pin a commit for reproducibility)
# remotes::install_github('bvieth/powsimR')
```

## Quick Start

Tell your AI agent what you want to do:
- "How many replicates do I need to detect 1.5-fold changes in RNA-seq at 80% power?"
- "Estimate marginal power across the transcriptome at FDR 0.05 by simulation"
- "Should I sequence deeper or add samples on a fixed budget?"
- "Estimate power for my scRNA-seq disease-versus-control comparison"
- "Is reporting observed power a valid way to explain my non-significant result?"

## Example Prompts

### Bulk RNA-seq Power

> "I'm comparing drug-treated versus control cells, expect biological CV around 0.3, and want to detect 1.5-fold changes. How many replicates for 80% power, and can you confirm with simulation?"

> "Calculate the marginal power of a 5-versus-5 RNA-seq design at FDR 0.05 using PROPER with pilot dispersions."

### Single-cell

> "How does power scale with number of patients versus number of cells per patient for a scRNA-seq differential expression study?"

### Budget and Design Tradeoffs

> "I have budget for either 20M reads on 4 samples or 10M reads on 8 samples per group. Which gives more power to detect DE genes?"

> "A reviewer says my study is underpowered because the observed power was 0.3. How should I respond?"

## What the Agent Will Do

1. Identify the assay and whether pilot data exist.
2. Choose closed-form (RNASeqPower) for a quick number or simulation (PROPER/powsimR) for the reported figure.
3. Report power as a function of replicate number at the target FDR, not a single transcriptome-wide value.
4. Advise on the depth-versus-replicate allocation under budget constraints.
5. Flag misuse of observed power and recommend effect-size/CI reporting instead.

## Tips

- Power is per-gene; report marginal power at a target FDR, not one number for the whole transcriptome.
- Use simulation from the mean-dispersion trend for the reported figure; closed-form for a quick sanity check.
- Estimate the coefficient of variation or dispersion from pilot data; a literature value can be off by a factor of two.
- Past roughly 10-20 million mapped reads, add biological replicates rather than depth.
- Power to the minimum biologically meaningful effect, never the pilot-observed effect.
- Never use observed (post-hoc) power to interpret a non-significant result; report the confidence interval.

## Related Skills

- sample-size - The inverse problem: minimum replicates for a target power at a target FDR
- randomization-blocking - The experimental unit defines what is replicated
- batch-design - Account for batch/blocking factors in the power model
- differential-expression/deseq2-basics - Estimating dispersions from pilot data
- single-cell/preprocessing - Pseudobulk model underlying scRNA-seq power
- clinical-biostatistics/power-and-sample-size - Power for regulated clinical-trial endpoints
