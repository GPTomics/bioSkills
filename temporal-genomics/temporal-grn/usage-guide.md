# Temporal Gene Regulatory Network Inference - Usage Guide

## Overview

Infers directed, time-delayed regulatory relationships from bulk time-series expression data. Applies Granger causality, dynGENIE3 (ODE-based tree regression), and dynamic Bayesian networks to identify how transcription factor activity propagates through gene regulatory networks across time.

## Prerequisites

### Python
```bash
pip install statsmodels pandas numpy networkx matplotlib
```

### R
```r
install.packages('bnlearn')
devtools::install_github('vahuynh/dynGENIE3/dynGENIE3R')
```

### Data Requirements
- Time-series expression matrix with at least 6 timepoints (more improves lag estimation)
- Known or predicted transcription factor list (AnimalTFDB, PlantTFDB, or similar)
- Multiple biological replicates improve dynGENIE3 derivative estimation
- Evenly spaced timepoints preferred for Granger causality; uneven spacing acceptable for dynGENIE3

## Quick Start

Tell the AI agent what to infer:
- "Infer regulatory relationships between my TFs and targets from time-series expression data"
- "Run Granger causality to find time-delayed gene regulation"
- "Build a dynamic gene regulatory network from my temporal RNA-seq data"
- "Compare regulatory networks between conditions over time"

## Example Prompts

### Granger Causality
> "I have 12 timepoints of RNA-seq data. Test Granger causality between my transcription factors and target genes to find regulatory links."

> "Run pairwise Granger causality tests on my time-series expression data with lag 1 and lag 2."

### dynGENIE3
> "Use dynGENIE3 to infer a gene regulatory network from my time-series expression data with known TF regulators."

> "I have 3 biological replicates of a developmental time course. Run dynGENIE3 to identify key regulators."

### Dynamic Bayesian Networks
> "Learn a dynamic Bayesian network structure from my temporal expression data using hill-climbing with BIC."

> "Bootstrap a dynamic Bayesian network to find confident regulatory edges between my TFs and targets."

### Network Comparison
> "Compare the regulatory networks between treated and control time-course experiments. Which edges are gained or lost?"

> "Track network rewiring across my developmental stages using Jaccard similarity of edge sets."

## What the Agent Will Do

1. Load time-series expression data and TF list
2. Check stationarity (for Granger) and preprocess as needed
3. Run selected inference method (Granger, dynGENIE3, or DBN)
4. Filter significant edges by p-value, importance score, or bootstrap confidence
5. Build directed adjacency matrix of regulatory relationships
6. Optionally compare networks across conditions
7. Generate network visualizations with edge weights
8. Export ranked edge list and adjacency matrix

## Tips

- Granger causality requires stationarity; apply first differencing if the ADF test fails (p > 0.05)
- maxlag for Granger should satisfy n > 3 * maxlag; with 12 timepoints, maxlag=3 is the maximum
- dynGENIE3 benefits strongly from multiple biological replicates; 3+ replicates recommended
- Restrict regulators to known TFs for cleaner networks; genome-wide inference is noisy
- DBN bootstrap with R=200 and strength threshold 0.7 provides conservative but reliable edges
- Combine methods: edges detected by both Granger and dynGENIE3 are higher confidence
- Use Jaccard similarity < 0.3 as a rough threshold for substantial network rewiring

## Related Skills

- gene-regulatory-networks/coexpression-networks - Static co-expression networks
- gene-regulatory-networks/scenic-regulons - Single-cell regulon inference with pySCENIC
- gene-regulatory-networks/differential-networks - Condition-specific network comparison
- data-visualization/network-visualization - Network plotting with NetworkX and Cytoscape
