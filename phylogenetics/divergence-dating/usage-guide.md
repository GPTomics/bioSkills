# Divergence Dating - Usage Guide

## Overview

Divergence dating converts a rooted, branch-length tree into absolute node ages under a molecular-clock model. The central idea: because a branch length is the product rate x time, sequences alone fix only the relative rate-time tree, and calibrations are the one thing that turns it into millions of years. A node-age posterior is therefore dominated by the calibration prior and how the tree prior and neighboring calibrations reshape it, far more than by the sequence data. This guide covers clock-model choice (strict, UCLN, autocorrelated, random local clocks), calibration strategy (node priors, soft bounds, tip-dating, fossilized birth-death), the dating engines (BEAST2, MCMCTree/PAML, TreePL, LSD2), and the two checks that separate a real date from a prior reflected back: the sample-from-prior run and the temporal-signal test.

## Prerequisites

```bash
# BEAST2 (includes BEAUti, LogCombiner, TreeAnnotator); add the SA package for FBD tip-dating
conda install -c bioconda beast2

# MCMCTree (part of PAML) for genome-scale approximate-likelihood dating
conda install -c bioconda paml

# TreePL (penalized likelihood, very large trees)
conda install -c bioconda treepl

# LSD2 ships inside IQ-TREE 2 as --date; TempEst is a separate GUI download (beast.community/tempest)
conda install -c bioconda iqtree

# Python helpers for control files and reading dated trees
pip install biopython dendropy
```

Conceptual prerequisites: a rooted tree with branch lengths (in substitutions), at least one justified calibration (a fossil minimum or sampling dates), and, for tip-dating, heterochronous samples with verified temporal signal.

## Quick Start

Tell your AI agent what you want to do:
- "Estimate divergence times for my phylogeny and report HPD intervals, not point ages"
- "Run my BEAST2 analysis sampling from the prior first and compare effective priors to the posterior"
- "Set up fossil calibrations as minima with soft bounds, not point constraints"
- "Prepare an MCMCTree control file for genome-scale dating with approximate likelihood"
- "Check whether my virus dataset has temporal signal before I tip-date it"
- "Date a 5000-taxon tree with TreePL, cross-validating the smoothing parameter"
- "Decide whether to use node calibration or the fossilized birth-death process"

## Example Prompts

### Clock Model Selection
> "I have a multi-gene primate alignment. Should I use a strict or relaxed clock, and how do I decide?"

> "How do I read the ucld.stdev posterior and the coefficient of variation to tell whether a strict clock suffices?"

> "My data span mammals and birds with very different rates. Which clock model fits, and how do I justify it?"

### Calibration Strategy
> "I have three insect fossil calibrations. Help me set lognormal-offset priors that treat each fossil as a minimum."

> "Set up MCMCTree soft-bound calibrations with 2.5% tails for a vertebrate phylogeny."

> "I have a rich fossil record and a morphological matrix. Should I use the fossilized birth-death process instead of node priors?"

> "I only have secondary calibrations from a previous study. How do I avoid laundering their uncertainty?"

### Effective Prior and Convergence
> "Run my BEAST2 model from the prior only and tell me whether the posterior on each calibrated node differs from the effective prior."

> "My clock-rate ESS is below 200 in Tracer. How do I fix convergence?"

> "Two independent MCMCTree runs disagree on a deep node age. How do I diagnose it?"

### Tip-Dating (Viruses, Ancient DNA)
> "Check temporal signal in my viral dataset with a root-to-tip regression and a date-randomization test before dating."

> "Set up tip-dating in BEAST2 using collection dates and a coalescent tree prior."

> "Run a fast LSD2 dating pass with confidence intervals on my large outbreak tree."

### Engine Choice and Interpretation
> "I have 500 loci. Is BEAST2 or MCMCTree the right engine, and why?"

> "My divergence-time credible intervals are very wide on both rates and the root. What is wrong?"

> "The posterior on one node looks identical to the prior. Does that mean the data are uninformative there?"

## What the Agent Will Do

1. Confirm the input is a rooted, branch-length tree with at least one justified calibration, and route topology/rooting/model selection to the upstream skills.
2. Recommend a dating engine from dataset size and goal: BEAST2 for posteriors and fossils-as-tips, MCMCTree for genome-scale data, TreePL/LSD2 for very large or fast tip-dated trees.
3. Choose a clock model by testing clocklikeness (coefficient of variation / `ucld.stdev`), not by defaulting, and justify strict-vs-relaxed via marginal-likelihood comparison.
4. Design calibrations that treat each fossil as a minimum with a backward tail or soft bounds, justified per the Parham et al. 2012 checklist, and prefer the fossilized birth-death process when several fossils exist.
5. For heterochronous data, run a TempEst root-to-tip regression and a date-randomization test, and refuse to date if there is no temporal signal.
6. Run the analysis sampling from the prior first (`sampleFromPrior="true"` / `usedata=0`) and compare specified-vs-effective-vs-posterior on every calibrated node.
7. Run at least two independent chains, check ESS > 200 and chain agreement in Tracer, and discard adequate burn-in.
8. Report median + 95% HPD per node, the effective-prior-vs-posterior comparison, and all model and calibration choices: never a bare point, and never a PL date without a bootstrap CI.

## Tips

- Always run from the prior first; topological truncation and neighboring calibrations make the effective prior differ from the density you typed, and a tight posterior is often just that truncated prior reflected back.
- Treat every fossil as a minimum, never a point; a near-delta calibration forces a guaranteed-too-young, falsely precise node.
- Use soft bounds (canonical 0.025 tail) rather than hard walls so one bad fossil or an over-tight maximum cannot dominate.
- For viruses and ancient DNA, verify temporal signal (positive root-to-tip slope plus a passing date-randomization test) before any dating run; short-time-span data often have none.
- Check the coefficient of variation / `ucld.stdev`: abutting zero means a strict clock suffices; clearly above zero means the relaxed clock is doing necessary work.
- Prefer the fossilized birth-death process when several fossils are available; it uses all of them coherently instead of multiplying ad hoc node densities.
- If credible intervals do not shrink as you add sequence data, you are at the calibration floor: invest in better fossils, not more sites (check the MCMCTree infinite-sites plot).
- TreePL and r8s give point estimates only; cross-validate the smoothing parameter and bootstrap for confidence intervals, and never report a bare PL date.
- Always run at least two independent chains; agreement between chains is the most reliable convergence diagnostic, with ESS > 200 on every reported parameter.
- For deep divergences, watch for substitution saturation; use slower markers, amino acids, or codon models and check that per-partition branch lengths are reliable before trusting MCMCTree.

## Related Skills

- bayesian-inference - MCMC convergence, ESS, marginal-likelihood model comparison, and site-heterogeneous models
- modern-tree-inference - the rooted, branch-length ML tree and model selection that dating consumes
- tree-manipulation - rooting as a separate inference and the input tree required before dating
- tree-io - reading and writing the dated MCC tree without dropping HPD intervals on node ages
- epidemiological-genomics/phylodynamics - effective population size and Re estimation downstream of tip-dated trees
