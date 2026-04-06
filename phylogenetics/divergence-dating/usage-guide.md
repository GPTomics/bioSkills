# Divergence Dating - Usage Guide

## Overview

Molecular clock models and divergence time estimation using BEAST2, MCMCTree, and TreePL. Covers clock model selection (strict, UCLN, autocorrelated), fossil calibration strategies (node calibration, tip-dating, Fossilized Birth-Death), and practical workflows for inferring absolute divergence times on phylogenies.

## Prerequisites

```bash
# BEAST2 (includes BEAUti, LogCombiner, TreeAnnotator)
# Download from https://www.beast2.org/
# Or via conda:
conda install -c bioconda beast2

# MCMCTree (part of PAML)
conda install -c bioconda paml

# TreePL
conda install -c bioconda treepl

# Python utilities
pip install biopython ete3
```

## Quick Start

Tell your AI agent what you want to do:
- "Estimate divergence times for my phylogeny using BEAST2"
- "Set up fossil calibrations for my molecular clock analysis"
- "Choose between strict and relaxed clock models for my dataset"
- "Prepare an MCMCTree control file for genome-scale divergence dating"
- "Date my large phylogeny quickly with TreePL"
- "Compare node calibration versus tip-dating approaches"

## Example Prompts

### Clock Model Selection
> "I have a multi-gene alignment of primates. Should I use a strict or relaxed molecular clock?"

> "How do I interpret the ucld.stdev parameter from my BEAST2 relaxed clock analysis?"

> "My dataset spans mammals and birds. Which clock model handles this level of rate variation?"

### Fossil Calibration
> "I have three fossil calibrations for my insect phylogeny. Help me set up appropriate prior distributions in BEAST2."

> "How do I verify that my specified calibration priors match the effective priors in my BEAST2 analysis?"

> "Set up MCMCTree calibrations using soft bounds for a vertebrate phylogeny with five fossil constraints."

> "I only have secondary calibrations from a previous study. What precautions should I take?"

### BEAST2 Workflow
> "Walk me through a complete BEAST2 divergence dating analysis from alignment to dated tree."

> "My BEAST2 run has low ESS values for the clock rate. How do I improve convergence?"

> "Set up a Fossilized Birth-Death analysis in BEAST2 with morphological and molecular data."

### MCMCTree
> "Generate an MCMCTree control file for dating a 500-gene dataset with approximate likelihood."

> "My MCMCTree acceptance proportions are outside 20-40%. How do I tune the finetune parameters?"

> "Compare results from two independent MCMCTree runs to check convergence."

### TreePL
> "Date my phylogeny of 5000 species using TreePL with cross-validation."

> "How do I set up fossil calibrations in TreePL format?"

### Interpretation
> "My divergence time estimates have very wide credible intervals. What could be wrong?"

> "The posterior on one node age looks identical to the prior. Is my data uninformative for that node?"

## What the Agent Will Do

1. Assess the dataset size and complexity to recommend the appropriate tool (BEAST2, MCMCTree, or TreePL)
2. Select an appropriate clock model based on dataset characteristics and rate variation
3. Help design fossil calibration strategies with justified prior distributions
4. Generate configuration files (BEAST2 XML parameters, MCMCTree control files, TreePL config)
5. Set up the two-step approximate likelihood workflow for MCMCTree when working with large datasets
6. Recommend running from the prior first to verify effective calibration priors
7. Diagnose convergence issues (low ESS, non-stationarity, multimodality) and suggest solutions
8. Interpret results including credible intervals, rate variation, and prior-posterior comparisons
9. Flag common pitfalls such as prior interaction, secondary calibration error compounding, and substitution saturation

## Tips

- Always run from the prior first (no data) in BEAST2 to verify effective priors match expectations; calibration prior interactions through the tree prior can silently distort node age priors
- Use at least 2 independent runs with different random seeds in both BEAST2 and MCMCTree; comparing runs is the most reliable convergence diagnostic
- For genome-scale data (>100 loci), MCMCTree with approximate likelihood is orders of magnitude faster than BEAST2
- TreePL is ideal for exploratory analyses or very large trees (>1000 taxa) where Bayesian methods are too slow, but it does not produce posterior distributions
- UCLN (uncorrelated lognormal) is the default relaxed clock choice for most datasets; switch to autocorrelated rates only with specific biological justification
- Fossils provide minimum age constraints; maximum ages require additional geological or biogeographic arguments
- When using secondary calibrations (estimates from previous studies), widen the uncertainty substantially to avoid error compounding
- Check `ucld.stdev` (BEAST2) or `sigma2` (MCMCTree) posteriors to assess whether rate variation is present; near-zero values suggest a strict clock may suffice
- For deep divergences (>500 Ma), watch for substitution saturation; use slow-evolving markers or amino acid sequences
- The birth-death tree prior parameters can influence uncalibrated node ages; perform sensitivity analyses on these settings

## Related Skills

- bayesian-inference - MCMC convergence diagnostics and Bayesian fundamentals
- modern-tree-inference - ML tree topology estimation before dating
- species-trees - Species tree estimation (date after resolving species tree)
- tree-manipulation - Rooting trees (required input for dating)
- tree-io - Reading and converting tree file formats
- alignment/multiple-alignment - Prepare alignments for dating analyses
