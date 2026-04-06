# Bayesian Phylogenetic Inference - Usage Guide

## Overview

Bayesian phylogenetic inference using MrBayes, BEAST2, RevBayes, and PhyloBayes. Covers MCMC setup and convergence diagnostics, posterior probability interpretation, model comparison via stepping-stone sampling, prior sensitivity analysis, and site-heterogeneous models (CAT-GTR) for deep phylogenies prone to long branch attraction.

## Prerequisites

```bash
# MrBayes
conda install -c bioconda mrbayes

# BEAST2
conda install -c bioconda beast2
# Or download from https://www.beast2.org/

# Tracer (convergence diagnostics GUI)
# Download from https://github.com/beast-dev/tracer/releases

# PhyloBayes MPI
conda install -c bioconda phylobayes-mpi

# RevBayes
conda install -c bioconda revbayes
# Or download from https://revbayes.github.io/

# RWTY (R package for topological convergence)
# install.packages('rwty')

# Python dependencies for convergence scripts
pip install biopython numpy pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Run a Bayesian phylogenetic analysis with MrBayes on my alignment"
- "Check MCMC convergence for my BEAST2 analysis"
- "Compare two substitution models using stepping-stone sampling in MrBayes"
- "Run PhyloBayes with CAT-GTR to check for long branch attraction"
- "Parse my MrBayes parameter files and check ESS values"

## Example Prompts

### Basic Bayesian Analysis
> "Set up a MrBayes analysis for my Nexus alignment with GTR+I+G and two runs of 1 million generations"

> "Run BEAST2 on my alignment with bModelTest for automatic model averaging"

### Convergence Diagnostics
> "Check if my MrBayes runs have converged by examining ESS, PSRF, and trace plots"

> "My BEAST2 analysis has low ESS for the tree likelihood. What should I do?"

> "Assess topological convergence of my MrBayes runs using RWTY in R"

> "Parse the .p files from my two MrBayes runs and calculate ESS for all parameters"

### Model Comparison
> "Compare GTR+G vs GTR+I+G using stepping-stone sampling in MrBayes"

> "Calculate Bayes factors between two BEAST2 models using path sampling"

> "Set up nested sampling in BEAST2 to compare clock models"

### Prior Sensitivity
> "Run my MrBayes analysis sampling from the prior only to check if the data are informative"

> "Test whether my divergence time estimates are sensitive to the calibration prior"

### Deep Phylogenies and LBA
> "Run PhyloBayes with CAT-GTR on my amino acid alignment to test for long branch attraction"

> "My ML and Bayesian trees disagree on the placement of a long-branched taxon. How do I resolve this?"

> "Check convergence of my PhyloBayes chains with bpcomp and tracecomp"

### RevBayes
> "Write a RevBayes script for a GTR+G analysis with two independent runs"

> "Set up a custom model in RevBayes that is not available in MrBayes or BEAST2"

## What the Agent Will Do

1. Assess whether Bayesian analysis is appropriate (vs ML) based on the scientific question
2. Prepare the alignment in the correct format (Nexus for MrBayes, XML via BEAUti for BEAST2)
3. Set up the substitution model, priors, and MCMC parameters
4. Configure multiple independent runs for convergence assessment
5. Run the analysis or generate the command/script
6. Check convergence: ESS >= 200 for all parameters, PSRF < 1.01, trace plots, topological ESS
7. Summarize posterior trees with appropriate burn-in
8. Perform model comparison via stepping-stone sampling if requested (never harmonic mean)
9. Flag potential issues: low ESS, non-convergent runs, inflated posterior probabilities
10. Recommend PhyloBayes CAT-GTR when LBA is suspected

## Tips

- Always run at least two independent MCMC analyses and compare them for convergence
- Check ESS for ALL parameters, not just the likelihood. Branch lengths and model parameters matter too
- If ESS is low, run the chain longer; do not just increase thinning (samplefreq)
- BEAST2 does not run multiple chains by default; manually run independent analyses with different seeds
- Never use the harmonic mean estimator for model comparison; use stepping-stone or path sampling
- Posterior probability = 1.00 does not mean certainty; compare with ML bootstrap for robustness
- For deep phylogenies with suspected LBA, try PhyloBayes CAT-GTR before concluding on topology
- PhyloBayes can take weeks to converge on large datasets; check `bpcomp` maxdiff < 0.3
- MrBayes `nst=mixed` averages over substitution models during MCMC, avoiding model selection bias
- Report both posterior probabilities and bootstrap support when possible

## Related Skills

- modern-tree-inference - ML tree inference as the default starting point
- phylogenetics/divergence-dating - Molecular dating with BEAST2
- phylogenetics/species-trees - Coalescent methods for multi-locus datasets
- phylogenetics/tree-io - Parse MrBayes and BEAST2 output tree files
- phylogenetics/tree-visualization - Visualize posterior consensus trees
- alignment/multiple-alignment - Alignment quality affects all downstream inference
