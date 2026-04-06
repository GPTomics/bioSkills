---
name: bio-phylo-bayesian-inference
description: Run Bayesian phylogenetic analysis with MrBayes, BEAST2, RevBayes, and PhyloBayes including MCMC convergence diagnostics and model comparison. Use when needing posterior probability support, Bayesian model averaging, site-heterogeneous models for deep phylogenies, or formal model comparison via stepping-stone sampling.
tool_type: mixed
primary_tool: MrBayes
---

## Version Compatibility

Reference examples tested with: MrBayes 3.2.7+, BEAST2 2.7+, Tracer 1.7+, RevBayes 1.2+, PhyloBayes MPI 1.9+, RWTY (R package)

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mb --version`, `beast -version`, `rb --version`, `pb --version`
- Python: `pip show biopython` then `help(module.function)` to check signatures
- R: `packageVersion('rwty')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Bayesian Phylogenetic Inference

**"Run a Bayesian phylogenetic analysis"** -> Infer posterior distribution of trees and parameters via MCMC sampling, producing posterior probability support and enabling formal model comparison.
- CLI: `mb` (MrBayes), `beast` (BEAST2), `rb` (RevBayes), `pb`/`bpcomp` (PhyloBayes)
- Python: BioPython `Bio.Phylo` for parsing output trees; `arviz`/`pandas` for trace diagnostics

## When to Use Bayesian vs ML

| Factor | ML (IQ-TREE/RAxML-NG) | MrBayes | BEAST2 | RevBayes | PhyloBayes |
|--------|------------------------|---------|--------|----------|------------|
| Primary use | Topology, hypothesis testing | Posterior probabilities on topology | Divergence times, demographics, tip-dating | Non-standard models, full graphical model control | Deep phylogenies with compositional heterogeneity |
| Speed | Fast | Moderate | Slow | Moderate | Very slow |
| Ease of use | Easy | Easy | Moderate (BEAUti GUI) | Steep learning curve | Moderate |
| Model flexibility | Good | Good | Extensive (packages) | Maximum (scripted) | CAT/CAT-GTR |
| Default choice | Yes, for topology | When posterior probabilities needed | When time calibration needed (see divergence-dating) | When no existing model fits | When LBA from model misspecification suspected |
| Multi-run | Automatic | 2 runs x 4 chains default | Must run independently | Manual | Must run independently |

**Default recommendation:** Start with ML (modern-tree-inference). Move to Bayesian when posterior probabilities, divergence times, or site-heterogeneous models are specifically needed.

## MrBayes

### Basic Analysis

**Goal:** Obtain a Bayesian phylogeny with posterior probability support values.

**Approach:** Load a Nexus alignment into MrBayes, set the substitution model and priors, run two independent MCMC analyses with Metropolis-coupled chains (MC3), then summarize parameters and trees after discarding burn-in.

```
execute alignment.nex
lset nst=6 rates=invgamma
prset brlenspr=unconstrained:exp(10.0)
mcmc ngen=1000000 samplefreq=100 nchains=4 nruns=2
sump
sumt
```

Key settings:
- `nst=6 rates=invgamma`: GTR+I+G model. Use `nst=2` for HKY, `nst=1` for JC/F81
- `nchains=4 nruns=2`: Two independent runs each with 4 Metropolis-coupled chains (1 cold + 3 heated). This is the default and should not be reduced
- `samplefreq=100`: Sample every 100 generations. Adjust so total samples = ngen/samplefreq is at least 10,000
- `sump`: Summarizes parameter traces, reports ESS and PSRF
- `sumt`: Summarizes trees, builds majority-rule consensus with posterior probabilities

### MC3 (Metropolis-Coupled MCMC)

MrBayes runs Metropolis coupling by default: one cold chain explores the posterior while heated chains explore broader landscape and swap states with the cold chain. This improves mixing for multimodal posteriors. The acceptance rate for swaps between adjacent chains should be 20-70%; adjust `temp=0.2` (default) if swaps are too rare or too frequent.

### Mixed Models and Partitions

```
charset gene1 = 1-300
charset gene2 = 301-600
partition by_gene = 2: gene1, gene2
set partition = by_gene
lset applyto=(1) nst=6 rates=invgamma
lset applyto=(2) nst=2 rates=gamma
unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all)
prset applyto=(all) ratepr=variable
```

### Reversible-Jump Model Selection

MrBayes can average over all 203 time-reversible substitution models during MCMC:

```
lset nst=mixed rates=invgamma
```

This integrates over model uncertainty rather than fixing a single model. Preferred when model choice is uncertain.

## BEAST2

### Workflow

1. Prepare alignment in BEAUti (generates XML configuration)
2. Set substitution model, clock model, and tree prior in BEAUti
3. Run BEAST2 from command line
4. Check convergence in Tracer
5. Summarize trees with TreeAnnotator
6. Visualize in FigTree or equivalent

```bash
beast -threads 4 -seed 12345 analysis.xml
treeannotator -burnin 10 -heights median analysis.trees consensus.tree
```

### bModelTest for Bayesian Model Averaging

The bModelTest package averages over all 203 time-reversible nucleotide substitution models during MCMC, eliminating the need for separate model selection:

1. Install via BEAUti Package Manager
2. In BEAUti: Site Model tab -> select "bModelTest"
3. The analysis samples across models weighted by their posterior probability

This is analogous to MrBayes `nst=mixed` but implemented as a BEAST2 package.

### Running Multiple Independent Analyses

BEAST2 does NOT run multiple chains by default (unlike MrBayes). Always run at least two independent analyses with different seeds:

```bash
beast -threads 4 -seed 12345 analysis.xml
beast -threads 4 -seed 67890 analysis.xml
```

Then compare traces in Tracer to confirm both runs converge to the same posterior.

## MCMC Convergence Diagnostics

This is the most critical aspect of any Bayesian phylogenetic analysis. An unconverged analysis produces meaningless results regardless of the model or data.

### Effective Sample Size (ESS)

ESS measures the number of effectively independent samples after accounting for autocorrelation.

| ESS | Interpretation |
|-----|---------------|
| < 100 | Insufficient; results unreliable |
| 100-199 | Marginal; increase run length |
| >= 200 | Adequate (Tracer default threshold) |
| >= 625 | Conservative target for precise credible intervals |

**Critical rule:** If ESS < 200 for any parameter, run the chain longer. Do NOT simply increase thinning (samplefreq). Thinning discards information and does not improve ESS; it is only justified to reduce file size.

### Trace Plot Interpretation

A well-converged trace plot resembles a "hairy caterpillar," with rapid oscillation around a stable mean and no visible trend. Warning signs:
- Upward/downward trends: chain has not reached stationarity
- Steps or plateaus: chain stuck in local optima
- Different levels between runs: runs exploring different regions of parameter space
- Periodic oscillations: poor mixing, possibly due to correlated parameters

### PSRF (Potential Scale Reduction Factor)

MrBayes reports PSRF (Gelman-Rubin diagnostic) comparing variance within and between runs:

| PSRF | Interpretation |
|------|---------------|
| < 1.01 | Good convergence |
| 1.01-1.05 | Acceptable but monitor |
| > 1.05 | Not converged; run longer or check model |

### Topological Convergence (RWTY)

Convergence of continuous parameters (likelihood, branch lengths) does NOT guarantee topological convergence. The tree topology may still be poorly sampled even when ESS for likelihood is high.

**Goal:** Assess whether MCMC chains have adequately sampled tree space.

**Approach:** Use the RWTY R package to compute topological ESS and visualize treespace using multidimensional scaling of Robinson-Foulds distances between sampled trees.

```r
library(rwty)
trees_run1 <- load.trees('run1.t', type='nexus')
trees_run2 <- load.trees('run2.t', type='nexus')
rwty_output <- analyze.rwty(list(run1=trees_run1, run2=trees_run2), burnin=25)
makeplot.topology(rwty_output)
makeplot.treespace(rwty_output)
```

What to look for:
- Topological ESS should be >= 200 (same threshold as continuous ESS)
- Treespace plot: samples from independent runs should overlap in a single cluster
- If runs occupy different regions of treespace, they have not converged on topology
- Split frequency standard deviation (MrBayes `sumt`): average < 0.01 is excellent, < 0.05 acceptable

### Convergence Checklist (Minimum Requirements)

1. Run at least 2 independent MCMC analyses
2. Check ESS >= 200 for ALL parameters (not just likelihood)
3. Check PSRF < 1.01 for all parameters (MrBayes)
4. Examine trace plots for stationarity
5. Check topological ESS with RWTY
6. Verify independent runs converge to same posterior (overlapping treespace)
7. If any check fails, run longer. Do not just increase thinning

## Burn-in Selection

Burn-in removes initial samples collected before the chain reaches stationarity. Typical range: 10-25% of total samples.

**Examine trace plots** to determine appropriate burn-in rather than using a fixed percentage. The burn-in period ends where the trace stabilizes at its stationary distribution.

- Too-short burn-in: biases posterior estimates toward starting values
- Too-long burn-in: wastes samples but does not introduce bias
- When uncertain, err on the side of discarding more (conservative)

MrBayes default burn-in for `sump`/`sumt` is 25%. BEAST2 TreeAnnotator requires specifying burn-in explicitly.

## Prior Sensitivity Testing

### Sampling from the Prior

**Goal:** Determine whether the data are informative for each parameter or whether the posterior is dominated by the prior.

**Approach:** Run MCMC sampling from the prior only (no likelihood calculation), then compare prior and posterior distributions.

In MrBayes:
```
mcmc ngen=1000000 data=no
```

In BEAST2: check "Sample from prior" in BEAUti MCMC tab.

### Interpreting Prior vs Posterior

- Prior and posterior overlap substantially: data is uninformative for that parameter; results are prior-dependent
- Prior and posterior differ markedly: data is driving the estimate; results are robust to prior choice
- For divergence time calibrations: always check the joint time prior (marginal priors on individual nodes can interact to produce unexpected effective priors on other nodes)

### Sensitivity Analysis Protocol

1. Run the analysis with default priors
2. Run with alternative priors (e.g., exponential vs lognormal for branch lengths)
3. Compare posterior distributions
4. If posteriors differ substantially, report results under multiple prior schemes or justify the chosen prior biologically

## Bayesian Model Comparison

### Methods (Ranked by Reliability)

| Method | Reliability | Cost | Available In |
|--------|------------|------|-------------|
| Harmonic mean estimator | NEVER use | Low | MrBayes (legacy) |
| Path sampling (PS) | Good | High | MrBayes, BEAST2 |
| Stepping-stone sampling (SS) | Good | High | MrBayes, BEAST2 |
| Generalized stepping-stone (GSS) | Better than SS | High | BEAST2 (MODEL_SELECTION package) |
| Nested sampling (NS) | Good | Moderate | BEAST2 (NS package) |

The harmonic mean estimator is the default in older MrBayes versions. It is unreliable, typically overestimates marginal likelihoods, and can favor overly complex models. Never use it for model comparison.

### Stepping-Stone Sampling in MrBayes

```
ss ngen=1000000 diagnfreq=1000 nsteps=50 alpha=0.4
```

This estimates the marginal likelihood via a series of power posteriors between the prior and the posterior. The `alpha` parameter controls the spacing of stepping stones (0.4 is the recommended default).

### Bayes Factor Interpretation

Bayes factor (BF) = ratio of marginal likelihoods. Compute as 2 * ln(BF) for comparison:

| 2 * ln(BF) | Evidence |
|-------------|----------|
| < 2 | Not worth mentioning |
| 2-6 | Positive |
| 6-10 | Strong |
| > 10 | Decisive |

Equivalently using natural log: ln(BF) > 1.0 = positive, > 3.2 = strong, > 4.6 = decisive.

### BEAST2 Path/Stepping-Stone Sampling

Install the MODEL_SELECTION package via BEAUti Package Manager. Then modify the XML to use path sampling:

```bash
beast -threads 4 model_selection.xml
```

The output reports log marginal likelihood estimates for each model.

## Posterior Probability Caveats

- PP values can be inflated by model misspecification: an inadequate model concentrates posterior mass on a wrong tree, producing high PP for an incorrect topology
- PP = 1.00 does not mean absolute certainty; it means the posterior probability rounds to 1.00 given the model and data
- Star-tree paradox: under certain conditions, even random (non-phylogenetic) data can produce high posterior probabilities for resolved topologies
- Compare PP with ML bootstrap support: if a clade has PP = 1.00 but bootstrap < 50%, suspect model misspecification inflating the PP
- For robust support assessment, report both PP and bootstrap values; agreement between methods strengthens confidence

## PhyloBayes for Deep Phylogenies

### When to Use PhyloBayes

PhyloBayes implements the CAT and CAT-GTR models, which are site-heterogeneous mixture models that assign each site to a profile category. These models are the gold standard for mitigating long branch attraction (LBA) caused by model misspecification in deep phylogenies (e.g., animal phyla, eukaryote root, bacterial deep branches).

CAT-GTR is preferred over fixed-matrix models (LG, WAG) or even empirical mixture models (C60) when:
- Sequences span deep evolutionary time (>500 Mya)
- Compositional heterogeneity across sites is expected
- LBA is suspected or known to be a problem
- Other methods produce conflicting topologies depending on the model

### Running PhyloBayes

```bash
# Run two independent chains (mandatory for convergence assessment)
mpirun -np 8 pb_mpi -d alignment.phy -cat -gtr -x 10 5000 chain1 &
mpirun -np 8 pb_mpi -d alignment.phy -cat -gtr -x 10 5000 chain2 &

# After chains complete, check convergence
bpcomp -x 1000 chain1 chain2
tracecomp -x 1000 chain1 chain2
```

### PhyloBayes Convergence

| Diagnostic | Tool | Good | Acceptable | Poor |
|-----------|------|------|------------|------|
| maxdiff (bipartition) | `bpcomp` | < 0.1 | 0.1-0.3 | > 0.3 |
| rel_diff (continuous params) | `tracecomp` | < 0.1 | 0.1-0.3 | > 0.3 |
| effsize | `tracecomp` | > 300 | 100-300 | < 100 |

PhyloBayes chains are typically very slow to converge for large datasets. Running for weeks on a cluster is not unusual. If `maxdiff` remains high, consider:
- Running longer (most common solution)
- Removing fast-evolving or saturated sites
- Using `-dgam 4` to add discrete gamma rate variation
- Reducing dataset size (fewer taxa or genes) as a diagnostic

### Recoding for Convergence

Amino acid recoding (Dayhoff-6 or SR4) reduces the state space and can improve convergence and reduce compositional bias:

```bash
# Dayhoff-6 recoding in PhyloBayes
pb_mpi -d alignment.phy -cat -gtr -recode dayhoff6 -x 10 5000 chain_d6
```

Recoding loses information but can reveal whether results are driven by compositional signal versus phylogenetic signal.

## RevBayes

RevBayes implements a probabilistic programming language for phylogenetics. It provides maximum flexibility for specifying custom models but requires writing Rev scripts.

```rev
# Basic GTR+G analysis in Rev language
data <- readDiscreteCharacterData("alignment.nex")
taxa <- data.taxa()
n_taxa <- taxa.size()

# Substitution model: GTR
er ~ dnDirichlet(v(1,1,1,1,1,1))
pi ~ dnDirichlet(v(1,1,1,1))
Q := fnGTR(er, pi)

# Among-site rate variation: Gamma
alpha ~ dnExponential(1.0)
site_rates := fnDiscretizeGamma(alpha, alpha, 4)

# Tree and branch lengths
topology ~ dnUniformTopology(taxa)
for (i in 1:2*n_taxa-3) {
    br_lens[i] ~ dnExponential(10.0)
}
psi := treeAssembly(topology, br_lens)

# Phylogenetic CTMC
seq ~ dnPhyloCTMC(tree=psi, Q=Q, siteRates=site_rates, type="DNA")
seq.clamp(data)

# MCMC
moves = VectorMoves()
moves.append(mvSimplexElementScale(er, weight=3))
moves.append(mvSimplexElementScale(pi, weight=2))
moves.append(mvScale(alpha, weight=1))
moves.append(mvNNI(topology, weight=n_taxa))
moves.append(mvSPR(topology, weight=n_taxa/5))
for (i in 1:2*n_taxa-3) {
    moves.append(mvScale(br_lens[i]))
}

monitors = VectorMonitors()
monitors.append(mnModel(filename="output.log", printgen=100))
monitors.append(mnFile(psi, filename="output.trees", printgen=100))
monitors.append(mnScreen(printgen=1000))

mymcmc = mcmc(mymodel, monitors, moves, nruns=2)
mymcmc.run(generations=1000000)
```

RevBayes is ideal when existing software does not support the desired model (e.g., non-standard clock models, state-dependent diversification, custom priors).

## Related Skills

- modern-tree-inference - ML alternative for topology estimation
- phylogenetics/divergence-dating - BEAST2 for molecular dating (different workflow)
- phylogenetics/species-trees - Coalescent methods as alternative to concatenated Bayesian
- phylogenetics/tree-io - Read MrBayes/BEAST2 output trees
- phylogenetics/tree-visualization - Visualize posterior trees
- alignment/multiple-alignment - Alignment quality directly affects tree inference
