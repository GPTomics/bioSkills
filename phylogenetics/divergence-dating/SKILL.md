---
name: bio-phylo-divergence-dating
description: Estimate divergence times under molecular-clock models with BEAST2, MCMCTree/PAML, TreePL, and LSD2, framing a date as a product of the calibration prior and the clock model far more than of the sequence data. Covers why branch length = rate x time is nonidentifiable so only calibrations convert relative rate-time into absolute age; why the effective (marginal) prior on a calibrated node differs from the density specified, mandating a sample-from-prior run; the fossil-as-minimum rule, soft bounds, tip-dating, and the fossilized birth-death process; the temporal-signal check (TempEst root-to-tip regression + date-randomization) required before dating viruses or ancient DNA; and clock-model choice via the coefficient of variation. Use when dating nodes, calibrating with fossils or sampling dates, choosing a clock or dating engine, or routing topology to modern-tree-inference, posteriors to bayesian-inference, and rooting to tree-manipulation.
tool_type: mixed
primary_tool: BEAST2
---

## Version Compatibility

Reference examples tested with: BEAST2 2.7+, MCMCTree/PAML 4.10+, TreePL 1.0+, TempEst 1.5+, LSD2 (IQ-TREE 2.2+ `--date`).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `beast -version`, `mcmctree` (PAML), `treePL`, `iqtree2 --version` then the tool's `-help`/`--help` to confirm flags
- Python: `pip show biopython dendropy` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

BEAST2 prior-only is `sampleFromPrior="true"`; the MCMCTree equivalent is `usedata = 0`. MCMCTree calibration syntax is `B()/L()/U()` in the tree file, not `>`/`<`. FBD tip-dating needs the BEAST2 SA package.

# Divergence Time Estimation -- A Date Is Mostly the Calibration Prior, Not the Sequence

**"Estimate when these lineages diverged"** -> Convert a relative rate-time tree into absolute ages using external calibrations, then report the posterior, not a point.
- CLI: BEAST2 (BEAUti XML) for full posteriors, FBD, and tip-dating
- CLI: MCMCTree (PAML) for genome-scale data via approximate likelihood
- CLI: TreePL / LSD2 for very large trees (point estimates / fast phylodynamics)

Scope: converting a rooted, branch-length tree into absolute node ages -- clock models, calibrations, dating engines, and the temporal-signal and effective-prior checks. Topology, model selection, and support -> modern-tree-inference. Posterior distributions, MCMC convergence, and site-heterogeneous models -> bayesian-inference. The clock-based rooting concept and re-rooting -> tree-manipulation. Reading/writing the dated MCC tree without dropping HPD intervals -> tree-io. Phylodynamic population-size / Re estimation -> epidemiological-genomics/phylodynamics.

## The Single Most Important Modern Insight

A divergence date is a product of the calibration priors and the clock model far more than of the sequence data. A branch length is the product b = rate x time (expected substitutions per site); the likelihood depends on b alone, so the pair (rate, time) is nonidentifiable -- doubling every rate and halving every time leaves the likelihood unchanged. Sequences therefore carry information about the relative rate-time tree only, and calibrations are the one thing that converts it to absolute millions of years. The posterior on a node age is consequently dominated by the (often subjective) calibration prior and by how the tree prior and neighboring calibrations reshape it. Three load-bearing facts:

1. **More sequence data does not escape the calibration.** Adding genes sharpens the relative tree but the residual uncertainty in node ages is bounded below by calibration uncertainty; a precise posterior on a badly placed fossil is a precise estimate of the wrong age (dos Reis and Yang 2011; see the infinite-sites plot below).
2. **The honest practice is to run the MCMC sampling from the prior (no data) and report the effective prior.** Topological ordering and neighboring calibrations truncate each specified density, so the marginal (effective) prior the model actually uses is generically NOT the density typed in (Heled and Drummond 2012; Warnock et al. 2012). Report three numbers per calibrated node: specified prior, effective prior, posterior.
3. **A tight credible interval is not evidence the data nailed the age.** If the posterior coincides with the effective prior, the data added nothing and the prior is being reported back; tightness usually traces to a tight (or mutually truncating) prior, not to the molecular signal (Brown and Smith 2018).

## Clock Model Selection

The clock governs how substitution rate varies across branches; choosing it wrong biases dates and misstates their uncertainty. Choose by testing clocklikeness, not by defaulting either way.

| Model | Assumption | When | Diagnostic |
|-------|------------|------|------------|
| Strict | one rate for the whole tree | clocklike data: intraspecific, or short-timescale viral, when a clock test does not reject | most efficient; tightest justified CIs |
| UCLN (uncorrelated lognormal) | each branch rate drawn independently from a lognormal | field default for multi-species data with rate variation | the `ucld.stdev` / coefficient-of-variation diagnostic |
| UCED (uncorrelated exponential) | branch rates drawn from an exponential | larger, less-Gaussian rate swings | available; rarely the first choice |
| Autocorrelated (ACLN) | descendant rate centered on parent (Brownian log-rate; Thorne et al. 1998) | deep trees where rate is heritable (generation time, metabolism); MCMCTree default | rate-variance `sigma2` (MCMCTree) |
| Random local clocks | a few inferred, discrete rate-shift points | episodic / clade-specific rate shifts (Drummond and Suchard 2010) | estimates where and how many shifts |

Relaxed-clock work was established by Drummond et al. 2006 (uncorrelated relaxed clocks, "dating with confidence"). The single most useful BEAST2 relaxed-clock diagnostic is the coefficient of variation (CoV) of branch rates, derived from `ucld.stdev`: a CoV posterior abutting 0 (`ucld.stdev` near 0) means rates are effectively constant and a strict clock suffices (gain precision by simplifying); a CoV clearly above 0 with 0 excluded means the relaxed clock is doing necessary work and a strict clock would be falsely precise. If the `ucld.stdev` posterior just recovers its prior, the data cannot indicate how clocklike the lineages are -- report that. Use the CoV for a quick read, but decide strict-vs-relaxed formally by marginal-likelihood comparison (path sampling / stepping-stone).

## Calibration Strategy

Calibrations are the dominant input. The bedrock rule: a fossil is a MINIMUM, not a point -- a clade is at least as old as a fossil assigned to it, and the true divergence is older, so a near-delta prior on a fossil age forces a guaranteed-too-young, falsely precise node.

| Strategy | Encodes | When | Pitfall |
|----------|---------|------|---------|
| Lognormal (offset) node prior | offset = fossil minimum; true age = min + a modest+ gap | one well-justified fossil, modest gap expected | mean/SD chosen by feel sets the answer |
| Exponential (offset) node prior | firm minimum, agnostic about gap size | good minimum, weak idea of the maximum | long tail can pull the node very old |
| Uniform + soft bounds | min from fossil, max from absence/strat, both leaky (Yang and Rannala 2006) | a defensible minimum AND maximum | hard bounds (no tail) over-dictate |
| Total-evidence / tip dating | fossils as dated, morphologically scored tips (Ronquist et al. 2012) | morphology available; want data-driven fossil placement | the morphological clock is shaky |
| Fossilized birth-death (FBD) | all fossils as samples of one diversification process (Heath et al. 2014) | several fossils; want coherent calibration | needs lambda/mu/psi/rho; sampled-ancestor handling |
| Tip dates (sampling times) | calibration from collection dates | measurably-evolving populations: viruses, ancient DNA | requires a verified temporal-signal check first |

Soft bounds (Yang and Rannala 2006) make a bound a quantile, not a wall: a small canonical 0.025 tail of probability is allowed beyond each soft min/max so one bad fossil cannot dominate. Justify every fossil per Parham et al. 2012 (specimen identity, apomorphy-based placement, geochronologic basis, monophyly of the calibration clade, stated reasoning). The total-evidence approach (Ronquist et al. 2012) includes fossils as dated tips scored for morphology so the data, not the user, place each fossil. The fossilized birth-death process (Heath et al. 2014) is the modern coherent tree prior: it models speciation, extinction, and fossil sampling jointly, uses ALL fossils, allows sampled ancestors, and replaces the incoherent practice of multiplying ad hoc node densities. Prefer FBD when several fossils exist; use simple node densities only for one or two transparent constraints. Secondary calibrations (an age borrowed from another study) launder uncertainty -- never use a point, use the full distribution, and flag it.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| BEAST2 | Bouckaert et al. 2019 | full hierarchical Bayesian MCMC; FBD, tip-dating, total-evidence, phylodynamic priors | need a posterior, fossils-as-tips, or complex models; 10s-100s taxa |
| MCMCTree / PAML | dos Reis and Yang 2011 | approximate likelihood: a two-step BASEML gradient + Hessian, then MCMC over a Taylor approximation | genome-scale / many loci where full BEAST is infeasible |
| TreePL / r8s | Smith and O'Meara 2012; Sanderson 2002 | penalized-likelihood point estimate; roughness penalty `lambda` set by cross-validation | very large trees (1000s-10000s taxa); accept point estimates + bootstrap CIs |
| LSD2 / treedater | To et al. 2016 | least-squares dating; native tip-dating, very fast | huge tip-dated viral trees; fast rooting and sanity check before a Bayesian run |

MCMCTree's approximate likelihood does the expensive Felsenstein pruning once (computing branch-length MLEs, gradient, and Hessian per partition) rather than every MCMC step, which is what makes thousands of loci tractable; check that BASEML converged because unreliable per-partition branch lengths (saturated/short loci) corrupt the approximation. The infinite-sites plot (dos Reis and Yang 2011) plots posterior CI width against posterior mean age across nodes: under infinite data the relationship becomes linear through the origin, with the residual width set entirely by calibration uncertainty -- if real-data points already hug that line, more sequence data will NOT narrow the dates and the answer is better fossils, not more sites. TreePL and r8s return one number per node and produce NO uncertainty; a bare PL date is a failure mode, not a result -- bootstrap for CIs and cross-validate `lambda`.

## Temporal Signal Before Tip-Dating

When samples are collected at different times and the population evolves fast enough to accumulate measurable substitutions between dates (a measurably-evolving population: viruses, ancient DNA), the sampling-date differences themselves are the calibrations -- no fossil needed. Verifying temporal signal first is non-negotiable; a short-time-span dataset frequently has no signal yet a Bayesian run will return a tight, entirely prior-driven date.

- **Root-to-tip regression (TempEst, Rambaut et al. 2016):** build a rough ML tree, regress each tip's root-to-tip genetic distance against its sampling date. A genuine signal gives a positive slope (the slope estimates the rate); a negative slope means no usable signal or a wrong root. R^2 is exploratory only (tips are non-independent): treat near-zero / much-below ~0.2 as a red flag, not a formal pass. The x-intercept estimates the TMRCA as a sanity check; large-residual tips are date/contamination/recombination outliers to investigate before dating.
- **Date-randomization test (Duchene et al. 2015):** the formal test. Re-run the dating with tip dates shuffled many times; the real-data clock-rate estimate must fall OUTSIDE the distribution of randomized replicates (no CI overlap). If the real estimate sits inside that cloud, there is no temporal signal and any date is an artifact of the prior.

## Run Prior-Only to Expose the Effective Prior

**Goal:** Determine whether the molecular data inform each calibrated node, or whether the reported posterior is just a truncated prior reflected back.

**Approach:** Run the same model first with no sequence data to obtain the effective (marginal) prior, then with data; compare specified-vs-effective-vs-posterior on every calibrated node.

```bash
# BEAST2: edit the XML so the MCMC samples from the prior only (no likelihood)
# set <run ... sampleFromPrior="true"> (BEAUti: MCMC panel, "Sample From Prior")
beast -seed 1 -prefix prioronly prioronly.xml      # effective prior on every node
beast -seed 1 -prefix withdata  withdata.xml       # full posterior

# MCMCTree: usedata=0 gives the effective prior; usedata=2 the approx-likelihood posterior
mcmctree mcmctree_prior.ctl    # control file has usedata = 0
mcmctree mcmctree_post.ctl     # control file has usedata = 2
```

```python
from Bio import Phylo

prior = Phylo.read('prioronly.mcc.tree', 'nexus')   # effective prior summary
post = Phylo.read('withdata.mcc.tree', 'nexus')      # posterior summary
for c_prior, c_post in zip(prior.get_nonterminals(), post.get_nonterminals()):
    # if the posterior median and HPD ~ the effective prior, the data did not inform this node
    print(c_prior.confidence, c_post.confidence)     # compare per-node summaries side by side
```

## Check Temporal Signal Before Tip-Dating

**Goal:** Confirm a heterochronous (virus / ancient-DNA) dataset actually contains clock signal before committing to a Bayesian tip-dated run.

**Approach:** Regress root-to-tip distance on sampling date (positive slope, sane intercept, outliers flagged), then run a date-randomization test; only date if the real estimate sits outside the randomized cloud.

```bash
# Build a quick ML tree to feed TempEst (modern-tree-inference)
iqtree2 -s seqs.fa -m GTR+G -T AUTO --prefix rttree
# TempEst (GUI): load rttree.treefile + a tab file of tip sampling dates;
# read the root-to-tip regression -- require a POSITIVE slope; inspect R^2 and residual outliers.

# Fast non-Bayesian tip-dating + CI as a cross-check (LSD2 via IQ-TREE)
iqtree2 -s seqs.fa -m GTR+G --date dates.tsv --date-ci 100 --prefix lsd2   # dates.tsv: tip <tab> date
```

## Per-Method Failure Modes

### Effective Prior Is Not the Specified Prior, Unchecked
**Trigger:** Several calibration densities plus a tree prior (Yule / birth-death / FBD), run straight to the posterior.
**Mechanism:** Every node must be older than its descendants, so a parent and child density truncate each other, and the tree prior is multiplied in; the marginal prior can look nothing like either typed density (Heled and Drummond 2012; Warnock et al. 2012).
**Symptom:** A tight posterior credible interval is read as "the data nailed it," when it is really the (truncated) prior.
**Fix:** Always run prior-only (`sampleFromPrior="true"` / `usedata=0`); report specified-vs-effective-vs-posterior per node; if posterior ~ effective prior, the data did not inform it.

### Fossil Treated as a Point, Not a Minimum
**Trigger:** A near-delta calibration density centered on a fossil age.
**Mechanism:** A fossil only bounds a clade from below; the true origin is older by an unknown gap, so a point prior forces a guaranteed-too-young age.
**Symptom:** Falsely precise, systematically too-young dates that propagate across the tree.
**Fix:** Use the fossil as the offset/minimum with a backward tail (lognormal/exponential or soft bounds); never a point.

### No Temporal-Signal Check Before Tip-Dating
**Trigger:** Tip-dating a short-time-span virus or ancient-DNA dataset without TempEst + a date-randomization test.
**Mechanism:** With too little accumulated substitution between sampling dates, the data carry no rate information and the prior drives the date.
**Symptom:** Plausible-looking but entirely prior-driven dates; tight HPDs on data that cannot support them.
**Fix:** Root-to-tip regression (positive slope) AND a date-randomization test (real estimate outside the randomized cloud) before any dating run.

### Penalized-Likelihood Point Estimate Reported With No Uncertainty
**Trigger:** A TreePL / r8s date reported as a single number.
**Mechanism:** PL maximizes a penalized likelihood and returns a point; it produces no posterior or CI, and `lambda` controls how clocklike the tree is forced to be.
**Symptom:** "Clade X is 45 Ma" with no interval, and a `lambda` chosen by default rather than cross-validation.
**Fix:** Cross-validate `lambda` (TreePL `prime` + `cv`); bootstrap sites/input trees and re-run to get CIs; never report a bare PL date.

### Over-Tight Calibrations Drive the Posterior
**Trigger:** One or two narrow calibration densities dominating the timescale.
**Mechanism:** A narrow density propagates through the clock and tree prior to set ages everywhere; the data barely move them.
**Symptom:** The posterior barely differs from the prior, and conclusions reverse when a single calibration is tweaked.
**Fix:** Widen / soften bounds, sensitivity-analyze each calibration one at a time, and prefer FBD for coherence across many fossils.

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| ESS (posterior, prior, likelihood, every reported parameter) | > 200 | Rambaut et al. 2018 (Tracer); below ~100 unusable |
| Independent MCMC chains | >= 2, posteriors must overlap | convergence cannot be judged from one chain |
| Burn-in discarded | >= 10% (confirm by trace, not rote) | standard practice; verify stationarity |
| Soft-bound tail probability | 0.025 per bound | Yang and Rannala 2006; MCMCTree `pL = pU = 0.025` |
| Root-to-tip R^2 (TempEst) | exploratory; near-zero / << ~0.2 = weak signal; positive slope mandatory | Rambaut et al. 2016 (tips non-independent, not a formal test) |
| Date-randomization test | real-data rate estimate outside the randomized distribution (no CI overlap) | Duchene et al. 2015 |
| Coefficient of variation of branch rates | abutting 0 -> strict adequate; clearly > 0 (0 excluded) -> relaxed needed | Drummond et al. 2006 |
| Infinite-sites plot | points on the linear CI-width-vs-age line -> more sites will not help | dos Reis and Yang 2011 |
| MCMCTree acceptance proportion | ~20-40% (target ~30%); tune `finetune` if outside | PAML practice |
| Smoothing `lambda` (TreePL/r8s) | set by cross-validation, never default | Sanderson 2002; Smith and O'Meara 2012 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| MCMCTree refuses to run | no root calibration | set `RootAge` in the control file or a calibration on the root node |
| MCMCTree calibration silently ignored | used `>`/`<` notation (parsing bug) | use `B()`/`L()`/`U()` in the tree file |
| Posterior ~ prior for a node age | data uninformative for that node | report it honestly; do not claim the data estimated the age |
| Wide CIs on both rate and root age | rate-time confounding from too few calibrations | add a well-justified calibration; check the rate-vs-root-age correlation |
| Times off by 100x in MCMCTree | unit confusion (no fixed time unit; the user picks one, 100 Myr conventional so ages are O(1)) | keep calibrations and `rgene_gamma` in that same unit; then 0.6 = 60 Ma |
| Dates conflict wildly with independent evidence | unjustified fossil placement | apply the Parham et al. 2012 checklist; recheck the assigned clade |
| Deep dates biased | substitution saturation at fast sites | use slower markers, amino acids, or codon models; remove saturated partitions |

## References

Bouckaert R, Vaughan TG, Barido-Sottani J, Duchene S, Fourment M, et al. 2019. BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. *PLoS Computational Biology* 15(4):e1006650.
Drummond AJ, Ho SYW, Phillips MJ, Rambaut A. 2006. Relaxed phylogenetics and dating with confidence. *PLoS Biology* 4(5):e88.
Thorne JL, Kishino H, Painter IS. 1998. Estimating the rate of evolution of the rate of molecular evolution. *Molecular Biology and Evolution* 15(12):1647-1657.
Drummond AJ, Suchard MA. 2010. Bayesian random local clocks, or one rate to rule them all. *BMC Biology* 8:114.
Heath TA, Huelsenbeck JP, Stadler T. 2014. The fossilized birth-death process for coherent calibration of divergence-time estimates. *PNAS* 111(29):E2957-E2966.
Ronquist F, Klopfstein S, Vilhelmsen L, Schulmeister S, Murray DL, Rasnitsyn AP. 2012. A total-evidence approach to dating with fossils, applied to the early radiation of the Hymenoptera. *Systematic Biology* 61(6):973-999.
Dos Reis M, Yang Z. 2011. Approximate likelihood calculation on a phylogeny for Bayesian estimation of divergence times. *Molecular Biology and Evolution* 28(7):2161-2172.
Yang Z, Rannala B. 2006. Bayesian estimation of species divergence times under a molecular clock using multiple fossil calibrations with soft bounds. *Molecular Biology and Evolution* 23(1):212-226.
Sanderson MJ. 2002. Estimating absolute rates of molecular evolution and divergence times: a penalized likelihood approach. *Molecular Biology and Evolution* 19(1):101-109.
Smith SA, O'Meara BC. 2012. treePL: divergence time estimation using penalized likelihood for large phylogenies. *Bioinformatics* 28(20):2689-2690.
To T-H, Jung M, Lycett S, Gascuel O. 2016. Fast dating using least-squares criteria and algorithms. *Systematic Biology* 65(1):82-97.
Rambaut A, Lam TT, Carvalho LM, Pybus OG. 2016. Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). *Virus Evolution* 2(1):vew007.
Duchene S, Duchene D, Holmes EC, Ho SYW. 2015. The performance of the date-randomization test in phylogenetic analyses of time-structured virus data. *Molecular Biology and Evolution* 32(7):1895-1906.
Parham JF, Donoghue PCJ, Bell CJ, Calway TD, Head JJ, et al. 2012. Best practices for justifying fossil calibrations. *Systematic Biology* 61(2):346-359.
Heled J, Drummond AJ. 2012. Calibrated tree priors for relaxed phylogenetics and divergence time estimation. *Systematic Biology* 61(1):138-149.
Warnock RCM, Yang Z, Donoghue PCJ. 2012. Exploring uncertainty in the calibration of the molecular clock. *Biology Letters* 8(1):156-159.
Brown JW, Smith SA. 2018. The past sure is tense: on interpreting phylogenetic divergence time estimates. *Systematic Biology* 67(2):340-353.
Dos Reis M, Donoghue PCJ, Yang Z. 2016. Bayesian molecular clock dating of species divergences in the genomics era. *Nature Reviews Genetics* 17(2):71-80.
Rambaut A, Drummond AJ, Xie D, Baele G, Suchard MA. 2018. Posterior summarization in Bayesian phylogenetics using Tracer 1.7. *Systematic Biology* 67(5):901-904.

## Related Skills

- bayesian-inference - MCMC convergence, ESS, marginal-likelihood model comparison, and site-heterogeneous models
- modern-tree-inference - the rooted, branch-length ML tree and model selection that dating consumes
- tree-manipulation - rooting as a separate inference and the input tree required before dating
- tree-io - reading and writing the dated MCC tree without dropping HPD intervals on node ages
- epidemiological-genomics/phylodynamics - effective population size and Re estimation downstream of tip-dated trees
