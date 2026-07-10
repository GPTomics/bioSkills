# Bayesian Phylogenetic Inference - Usage Guide

## Overview

A Bayesian phylogeny is not a tree; it is a posterior distribution over trees conditioned on the data and on the priors, approximated by an MCMC that may not have converged. The deliverable is a distribution WITH diagnostics. This skill covers the decisions that drive the result: which tool (MrBayes, BEAST2, RevBayes, PhyloBayes-MPI) fits the question, how to prove MCMC convergence (ESS, PSRF, ASDSF, and the crucial scalar-vs-topology distinction), why posterior probabilities run higher than bootstrap and are overconfident under model misspecification, why the default branch-length prior inflates tree length, why model comparison must use stepping-stone sampling and never the harmonic mean, and when site-heterogeneous CAT-GTR is required at depth. The honest one-sentence framing: PP near 1 means the data are decisive GIVEN the model, not that the clade is true.

## Prerequisites

```bash
# MrBayes
conda install -c bioconda mrbayes

# BEAST2 (and BEAUti / TreeAnnotator / LogCombiner)
conda install -c bioconda beast2
# Or download from https://www.beast2.org/

# Tracer (visual ESS / trace diagnostics GUI)
# Download from https://github.com/beast-dev/tracer/releases

# PhyloBayes MPI (CAT / CAT-GTR for deep phylogeny)
conda install -c bioconda phylobayes-mpi

# RevBayes (custom graphical models)
conda install -c bioconda revbayes

# RWTY (R package for topological convergence)
# install.packages('rwty')

# Python for the convergence-check script
pip install biopython numpy pandas
```

Conceptual prerequisites: a tree is a point estimate, not an observation, and support is not accuracy: a posterior probability of 1.00 can still be wrong under model misspecification. Bayesian model comparison requires PROPER priors. At least two independent runs are mandatory, because PSRF and ASDSF are between-run metrics.

## Quick Start

Tell your AI agent what you want to do:
- "Run a Bayesian phylogenetic analysis with MrBayes on my Nexus alignment"
- "Check whether my MrBayes runs converged using ESS, PSRF, and ASDSF"
- "Compare GTR+G vs GTR+I+G using stepping-stone sampling in MrBayes"
- "Run PhyloBayes CAT-GTR to test whether a long-branched taxon is an LBA artifact"
- "Parse my two MrBayes .p files and compute ESS and PSRF for every parameter"

## Example Prompts

### Basic Bayesian Analysis
> "Set up a MrBayes analysis for my Nexus alignment with GTR+I+G, the compound-Dirichlet branch-length prior, and two runs of 10 million generations"

> "Run BEAST2 on my alignment with bModelTest for substitution-model averaging, then summarize with TreeAnnotator"

### Convergence Diagnostics
> "Check whether my MrBayes runs converged: ESS for all parameters, PSRF, ASDSF, and the trace plots"

> "My continuous parameters all have high ESS but I am not sure the topology converged: assess tree-space convergence with RWTY"

> "My BEAST2 analysis has low ESS for a clock rate. What should I change?"

### Model Comparison
> "Compare two partition schemes by stepping-stone marginal likelihood in MrBayes and report the Bayes factor"

> "Calculate Bayes factors between two BEAST2 models using the MODEL_SELECTION package"

### Priors
> "Switch my MrBayes branch-length prior off exp(10) and check whether the tree length was inflated"

> "Sample from the prior only to see whether the data are informative for each parameter"

### Deep Phylogenies and LBA
> "Run PhyloBayes CAT-GTR on my amino-acid alignment and check convergence with bpcomp and tracecomp"

> "My ML and Bayesian trees disagree on a long-branched taxon at depth: how do I tell whether it is an LBA artifact?"

### RevBayes
> "Write a RevBayes script for a GTR+G analysis with two independent runs"

> "Set up a custom hierarchical model in RevBayes that is not available in MrBayes or BEAST2"

## What the Agent Will Do

1. Decide whether Bayesian inference is the right tool (vs ML in modern-tree-inference, or dating in divergence-dating) from the question.
2. Pick the tool: MrBayes for a standard alignment, BEAST2 for time trees, RevBayes for custom models, PhyloBayes-MPI for deep / compositionally-heterogeneous data.
3. Set a multiple-hit-correcting substitution model and the compound-Dirichlet branch-length prior (not exp(10)).
4. Configure >= 2 independent runs and MC3 heated chains, with live ASDSF and auto-stop where available.
5. Run the analysis or emit the command/script.
6. Prove convergence: ESS > 200 for EVERY parameter, PSRF ~ 1.00, ASDSF < 0.01 (or bpcomp maxdiff < 0.1), trace stationarity, and tree-space agreement across runs.
7. Summarize the posterior trees with a data-driven burn-in.
8. Do model comparison via stepping-stone / path sampling if requested (never the harmonic mean), and interpret with Kass-Raftery thresholds.
9. Flag overconfident posteriors: PP = 1.0 that disagrees with the bootstrap, or high PP on a near-zero internode (star-tree paradox).
10. Recommend CAT-GTR (PhyloBayes) or PMSF (IQ-TREE) when LBA from compositional heterogeneity is suspected at depth.

## Tips

- Always run at least two independent runs from different starting trees; a single run cannot demonstrate convergence (PSRF and ASDSF are between-run metrics).
- Check ESS for EVERY parameter, not just the likelihood: a nuisance parameter can have ESS=12 while the likelihood reads 2000.
- Scalar convergence does not prove topology convergence. Always add a tree diagnostic: ASDSF (MrBayes), bpcomp maxdiff (PhyloBayes), or the RWTY tree-space MDS plot.
- If ESS is low, run more generations; do not just raise samplefreq, which discards information without raising ESS.
- Do not trust a default exp(10) branch-length prior on datasets with many taxa or long branches; use the compound Dirichlet (gammadir) and compare tree length to an ML reference.
- Never select a model with the harmonic-mean estimator; use stepping-stone or path sampling, and make sure all priors are proper.
- PP = 1.0 means "decisive under this model", not truth. When PP and bootstrap disagree sharply, distrust the more confident number.
- A high PP on a clade subtended by a near-zero internal branch is where Bayesian support is least trustworthy (the star-tree paradox).
- For deep phylogenies with suspected LBA, run PhyloBayes CAT-GTR (or IQ-TREE PMSF) before concluding on the topology; budget days-to-weeks and gate on maxdiff < 0.1.
- BEAST2 does not run multiple chains by default; launch independent seeds and combine with LogCombiner before inspecting in Tracer.

## Related Skills

- modern-tree-inference - ML inference, bootstrap/SH-aLRT support, concordance factors, and IQ-TREE C60/PMSF
- divergence-dating - BEAST2 clock models, calibration, and time-scaled posteriors
- species-trees - coalescent species-tree estimation when ILS dominates
- tree-io - reading and summarizing MrBayes/BEAST2 output trees without dropping posteriors and HPDs
