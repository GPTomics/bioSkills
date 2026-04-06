---
name: bio-phylo-divergence-dating
description: Estimate divergence times using molecular clock models with BEAST2, MCMCTree, and TreePL. Use when dating speciation events, calibrating phylogenies with fossils, choosing between strict and relaxed clock models, or estimating evolutionary rates across lineages.
tool_type: mixed
primary_tool: BEAST2
---

## Version Compatibility

Reference examples tested with: BEAST2 2.7+, MCMCTree (PAML 4.10+), TreePL 1.0+, Python 3.9+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `beast -version` then `beast -help` to confirm flags
- CLI: `mcmctree --help` or check PAML documentation for parameter names
- CLI: `treePL` to confirm installation
- Python: `pip show biopython` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Divergence Time Estimation

**"Estimate when lineages diverged"** -> Combine phylogenetic trees with fossil calibrations and molecular clock models to infer absolute divergence times.
- CLI: BEAUti + BEAST2 for moderate data with complex models
- CLI: MCMCTree (PAML) for genome-scale data with approximate likelihood
- CLI: TreePL for very large trees with penalized likelihood

## Clock Model Selection

| Model | Assumption | When to Use |
|-------|------------|-------------|
| Strict clock | Constant rate across all lineages | Only when a molecular clock test does not reject clocklike behavior; rarely appropriate for inter-species data |
| UCLN (uncorrelated lognormal) | Branch rates drawn independently from lognormal | **Default choice**; rates vary across lineages without parent-child correlation |
| Uncorrelated exponential | Rates from exponential distribution | More extreme rate variation; when UCLN coefficient of variation is very high |
| Autocorrelated (ACLN) | Rates correlated with parent branch | When evolutionary rates change gradually (e.g., body size evolution driving rate changes) |
| Random local clock | Few distinct rate regimes | When a small number of rate shifts are suspected rather than continuous variation |

### Clock Model Diagnostics

Under UCLN in BEAST2, the posterior of `ucld.stdev` measures clocklikeness:

| ucld.stdev | Interpretation |
|------------|----------------|
| ~0 | Data is clocklike; strict clock may suffice |
| 0.1-0.5 | Moderate rate variation; UCLN appropriate |
| >1.0 | Extreme rate variation; consider model adequacy or data partitioning |

In MCMCTree, the equivalent diagnostic is `sigma2` under clock=2 (independent rates). Near-zero values suggest a strict clock; very large values indicate substantial rate heterogeneity.

## Fossil Calibration Strategies

### Node Calibration (Traditional)

Prior distributions constrain the age of internal nodes based on fossil evidence:

| Distribution | When to Use | Parameterization |
|--------------|-------------|------------------|
| Lognormal | Most common; fossil sets minimum, soft maximum | Offset = fossil minimum age, mean and stdev control tail |
| Exponential | Strong belief divergence is close to fossil age | Rate parameter controls decay from minimum |
| Uniform | Only minimum and maximum known | Dangerous with narrow bounds; truncation distorts posteriors |
| Normal | Independently known divergence age (e.g., biogeographic event) | Mean = best estimate, stdev = uncertainty |

### Critical Calibration Pitfalls

These are the most common sources of error in divergence dating, and each can silently produce incorrect results:

1. **Prior interaction**: Multiple calibration priors interact through the tree prior (birth-death or Yule). The effective prior on a node may differ substantially from the specified prior. ALWAYS run the analysis sampling from the prior only (no data) first and compare effective priors to specified priors. If they differ, adjust calibration distributions.

2. **Truncation effects**: When a calibration density extends beyond the root age, the density gets truncated and renormalized, distorting the effective prior. This is especially problematic for uniform and lognormal priors near the root.

3. **Secondary calibrations**: Using divergence time estimates from a previous study as calibrations compounds uncertainty. The original confidence intervals become point-like constraints in the new analysis. When secondary calibrations are unavoidable, use substantially wider priors than the original confidence intervals.

4. **Fossil placement**: The single largest source of calibration error is incorrect phylogenetic placement of fossils. A fossil assigned to the wrong node can bias all downstream estimates. Verify fossil placement against current morphological and phylogenetic evidence.

5. **Minimum vs maximum confusion**: Fossils provide minimum age constraints (the lineage is at least as old as the fossil). Maximum ages require additional geological or biogeographic arguments and are harder to justify.

## Tip-Dating and Fossilized Birth-Death (FBD)

### FBD vs Node Calibration

| Factor | Node Calibration | Fossilized Birth-Death |
|--------|-----------------|----------------------|
| Fossil usage | Few well-placed fossils | All available fossils as tips |
| Prior distributions | Ad hoc (lognormal, uniform) | Coherent birth-death framework |
| Morphological data | Not required | Beneficial but optional |
| Extinct taxa | Excluded | Included as tips |
| Setup complexity | Lower | Higher (morphological matrix, FBD parameters) |
| Software | BEAST2, MCMCTree | BEAST2 (SA package) |

**When to prefer FBD**: Rich fossil record available, morphological data matrix exists, wanting to include extinct taxa, or concerned about ad hoc prior distributions.

**When to prefer node calibration**: Few well-placed fossils, no morphological data, simpler setup needed, or working with MCMCTree.

## Software Decision Framework

| Scenario | Recommended Tool | Rationale |
|----------|-----------------|-----------|
| Genome-scale data (>100 loci) | MCMCTree | Approximate likelihood is 100-1000x faster than exact |
| Moderate data, complex models | BEAST2 | Richest model ecosystem, FBD support |
| Very large tree (>1000 taxa) | TreePL | Penalized likelihood scales to 10,000+ taxa |
| Fossils as tips needed | BEAST2 + SA package | Only tool with FBD tip-dating |
| Mixture models + dating | IQ2MC (IQ-TREE 3 + MCMCTree) | Combines mixture model fit with Bayesian dating |
| Quick exploratory dating | TreePL | Fast; correlates well with Bayesian estimates |

## BEAST2 Dating Workflow

**Goal:** Estimate divergence times with calibrated relaxed clock in BEAST2.

**Approach:** Prepare XML via BEAUti, run BEAST2, verify convergence, summarize results.

### Step-by-Step Workflow

1. **Prepare alignment** in BEAUti: import FASTA/NEXUS, set data type
2. **Set substitution model**: use bModelTest for model averaging (recommended) or fix a model
3. **Set clock model**: UCLN is the default choice; set initial clock rate based on prior knowledge
4. **Set tree prior**: Birth-Death for extant-only data, FBD when including fossils as tips
5. **Add fossil calibrations**: specify prior distributions with justified parameters for each calibrated node
6. **Run from prior first**: sample from the prior only (uncheck "use data" or set `sampleFromPrior="true"`) to verify effective priors match expectations
7. **Run full analysis**: at least 2 independent runs with different random seeds
8. **Check convergence**: Tracer: all ESS values >= 200; compare posterior traces between runs
9. **Check topological convergence**: use RWTY or similar
10. **Combine and summarize**: LogCombiner (remove burn-in from each run), then TreeAnnotator for MCC tree

### BEAST2 Clock Rate Prior

**Goal:** Set an informable prior on the clock rate to avoid poor mixing.

**Approach:** Estimate the expected substitution rate from prior information and set the clock rate prior accordingly.

```bash
# Rough clock rate estimation:
# rate = divergence / (2 * time)
# For mammals: ~0.01 substitutions/site/Myr for nuclear genes
# For mitochondrial: ~0.02 substitutions/site/Myr
# Set lognormal prior with M = ln(expected_rate), S = 0.5-1.0
```

### Key BEAST2 Packages

| Package | Purpose |
|---------|---------|
| SA (sampled ancestors) | FBD tip-dating |
| bModelTest | Bayesian model averaging |
| SNAPP | Species tree from SNPs |
| ORC | Optimized relaxed clock operators |
| CoupledMCMC | Parallel tempering for difficult posteriors |

## MCMCTree Workflow

**Goal:** Date a large phylogeny using approximate likelihood in MCMCTree.

**Approach:** Prepare control file, run approximate likelihood calculation in two steps, then MCMC sampling.

### Two-Step Approximate Likelihood

The approximate likelihood method makes MCMCTree practical for genome-scale data:

```bash
# Step 1: Generate gradient and Hessian (in.BV file)
# Set usedata = 3 in control file, then run:
mcmctree mcmctree_step1.ctl

# Step 2: Run MCMC with approximate likelihood
# Set usedata = 2 in control file, then run:
mcmctree mcmctree_step2.ctl
```

### MCMCTree Calibration Format

MCMCTree uses specific notation in the tree file for calibrations. Use the `B()`, `L()`, `U()` format (NOT the `>` `<` notation, which has a known parsing bug in some versions):

| Notation | Meaning | Example |
|----------|---------|---------|
| `B(tL, tU, pL, pU)` | Soft bounds (lower, upper) | `B(0.6, 0.8, 0.025, 0.025)`: 60-80 Ma, 2.5% tail probability |
| `L(tL, p, c, pL)` | Minimum bound (truncated Cauchy) | `L(0.5, 0.1, 1.0, 0.025)`: min 50 Ma |
| `U(tU, pU)` | Maximum bound only | `U(1.0, 0.025)`: max 100 Ma |
| `S2(tL, p, c, pL)` | Skew-normal minimum bound | `S2(0.5, 0.1, 1.0, 0.025)` |
| `ST(tL, p, c, pL, tU, pU)` | Skew-t with upper bound | Full bounded calibration |

Times are in units of 100 Myr by convention (so 0.6 = 60 Ma), though any consistent unit works if `rgene_gamma` is adjusted.

### MCMCTree Practical Tips

- **Root calibration is mandatory.** MCMCTree will not run without a calibration on the root node; set via `RootAge` in the control file or directly in the tree file
- **Acceptance proportions** should be 20-40% (ideal ~30%); adjust `finetune` parameters if outside this range
- **Independent runs**: always run at least 2 with different seeds and compare posteriors
- **Prior vs posterior**: run with `usedata = 0` to sample from the prior; compare with posterior to assess data informativeness
- **BDparas**: birth-death parameters (lambda, mu, rho) affect the time prior on uncalibrated nodes; sensitivity analysis recommended
- **rgene_gamma**: sets the prior on the overall substitution rate; parameterize based on prior knowledge of the clade

### MCMCTree Control File Parameters

**Goal:** Generate a properly configured MCMCTree control file.

**Approach:** Set key parameters with appropriate defaults and document rationale for non-default choices.

```python
# Key MCMCTree control file parameters
params = {
    'seed': -1,                  # random seed (-1 for clock-based)
    'seqfile': 'alignment.phy',  # phylip format alignment
    'treefile': 'tree.nwk',      # newick tree with calibrations
    'outfile': 'mcmctree.out',   # output file
    'ndata': 1,                  # number of data partitions
    'usedata': 2,                # 0=no data, 1=exact, 2=approx, 3=generate in.BV
    'clock': 2,                  # 1=strict, 2=independent rates, 3=correlated rates
    'model': 4,                  # 0=JC69, 4=HKY85, 7=REV (GTR)
    'alpha': 0.5,                # gamma shape for rate variation among sites
    'ncatG': 4,                  # number of gamma categories
    'cleandata': 0,              # 0=keep ambiguity data, 1=remove
    'BDparas': '1 1 0.1',        # birth rate, death rate, sampling fraction
    'rgene_gamma': '2 20 1',     # gamma prior on rates: alpha, beta, dirichlet
    'sigma2_gamma': '1 10 1',    # gamma prior on sigma2 (rate drift)
    'print': 1,                  # output detail level
    'burnin': 50000,             # burn-in iterations
    'sampfreq': 50,              # sampling frequency
    'nsample': 20000,            # number of posterior samples
}
```

## TreePL

**Goal:** Rapidly estimate divergence times on very large phylogenies.

**Approach:** Penalized likelihood with cross-validation for smoothing parameter selection.

### TreePL Configuration

```
treefile = input.tre
smooth = 100
numsites = 1000
mrca = cal1 taxonA taxonB
min = cal1 50
max = cal1 80
mrca = cal2 taxonC taxonD
min = cal2 100
max = cal2 120
outfile = dated.tre
thorough
opt = 1
optad = 1
optcvad = 1
moredetail
cv
```

### TreePL Key Points

- **Cross-validation** (`cv` flag): determines optimal smoothing parameter; always enable for final analyses
- **Smoothing parameter**: controls the penalty on rate variation; too low = overfitting, too high = strict clock
- **Scalability**: handles trees with 10,000+ taxa where Bayesian methods are impractical
- **Correlation with Bayesian**: TreePL estimates correlate well (R-squared > 0.94) with BEAST2 and MCMCTree on benchmark datasets
- **Limitations**: no posterior distributions on node ages, no model comparison, limited clock model flexibility

## Common Pitfalls and Diagnostics

### Convergence Issues

| Problem | Symptom | Solution |
|---------|---------|----------|
| Low ESS | ESS < 200 in Tracer | Run longer, tune operators, simplify model |
| Non-stationarity | Trace trends up/down | Increase burn-in or run longer |
| Multimodality | Multiple peaks in marginal density | May indicate genuine uncertainty; try different priors or starting trees |
| Poor mixing | Trace stuck at single value | Adjust operator weights/tuning; check for conflicting calibrations |

### Rate and Time Confounding

Molecular clock analyses estimate the product rate x time. Calibrations break this confounding by fixing time at certain nodes. With too few calibrations, rate and time estimates can be poorly identified. Signs of confounding:

- Very wide credible intervals on both rates and times
- Strong negative correlation between rate and root age in posterior
- Posterior nearly identical to prior for some node ages (data is uninformative)

### Substitution Saturation

Deep divergences (>500 Ma) with fast-evolving markers can hit substitution saturation, where multiple substitutions at the same site erase phylogenetic signal. This can bias both topology and dates. Mitigation: use slow-evolving markers, amino acids instead of nucleotides, or codon models.

## Related Skills

- bayesian-inference - MCMC convergence diagnostics and Bayesian fundamentals
- modern-tree-inference - ML tree topology estimation (needed before dating)
- species-trees - Species tree estimation (date after resolving species tree)
- tree-manipulation - Rooting trees (required before dating)
- tree-io - Reading and converting tree files for input/output
- alignment/multiple-alignment - Prepare alignments for dating analyses
