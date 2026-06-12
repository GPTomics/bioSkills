---
name: bio-phylo-bayesian-inference
description: Frames Bayesian phylogenetics as approximating a posterior distribution over trees conditioned on data AND priors via an MCMC that must be proven to have converged, using MrBayes, BEAST2, RevBayes, and PhyloBayes-MPI. Covers why convergence (ESS, PSRF, ASDSF, topology vs scalar) is the load-bearing claim, why posterior probabilities are systematically higher than bootstrap and overconfident under model misspecification, why the default branch-length prior inflates tree length, why the harmonic-mean estimator must never select models (use stepping-stone), and when site-heterogeneous CAT-GTR is required at depth. Use when needing posterior clade support, model averaging, marginal-likelihood model comparison, or CAT models for deep phylogeny. Routes topology-only ML to modern-tree-inference, divergence times to divergence-dating, and tree summarization to tree-io.
tool_type: mixed
primary_tool: MrBayes
---

## Version Compatibility

Reference examples tested with: MrBayes 3.2.7+, BEAST2 2.7+, RevBayes 1.2+, PhyloBayes MPI 1.9+, Tracer 1.7+. R diagnostics: RWTY, coda. Python: BioPython 1.83+, NumPy, pandas.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mb --version`, `beast -version`, `rb --version`, `pb_mpi` then check the help banner to confirm flags
- R: `packageVersion('rwty')` then `?analyze.rwty` to verify parameters
- Python: `pip show biopython` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The MrBayes default branch-length prior changed at 3.2.3 from `unconstrained:exp(10)` to the compound `unconstrained:gammadir`; do not assume the old default on any version. `stoprule`/`stopval` auto-halt on ASDSF.

# Bayesian Phylogenetic Inference -- A Posterior Over Trees Conditioned on Priors; an Unconverged Run Is Confident Nonsense

**"Run a Bayesian phylogenetic analysis"** -> Approximate the posterior distribution over trees and parameters by MCMC, then prove convergence before trusting any clade posterior.
- CLI: `mb` (MrBayes), `beast`/BEAUti (BEAST2), `rb` (RevBayes), `pb_mpi`/`bpcomp`/`tracecomp` (PhyloBayes)
- Diagnostics: Tracer (visual ESS/trace), RWTY/coda (R), MrBayes `sump`/`sumt`, BioPython for parsing output trees

Scope: the posterior-over-trees + MCMC-convergence + prior + marginal-likelihood-model-comparison spine, plus site-heterogeneous CAT models. Topology-only ML, bootstrap, IC model selection -> modern-tree-inference. Divergence times, clock models, node/tip calibration, fossilized birth-death -> divergence-dating. Multispecies-coalescent species trees under ILS -> species-trees. Reading/converting/summarizing the annotated output trees without dropping posteriors and HPDs -> tree-io.

## The Single Most Important Modern Insight

A Bayesian phylogeny is not a tree. It is a posterior distribution over trees (topologies AND branch lengths AND model parameters) conditioned on the data AND on the priors, approximated by a finite MCMC sample that may not have converged. The deliverable is a distribution WITH diagnostics, not a point estimate; the consensus tree with clade posteriors is a lossy summary of that distribution. Three load-bearing facts:

1. The posterior is conditional on priors that may not have been chosen deliberately. Default priors are not neutral: the MrBayes default branch-length prior before 3.2.3 (`exp(10)` per branch) systematically inflates tree length, and that inflation can feed back into topology and node support (Brown 2010; Zhang 2012). A run that converged beautifully can be confidently wrong because it converged to the posterior implied by a bad prior. "Converged" answers "did the MCMC sample the target distribution?" -- not "was the target distribution the right one?".
2. The approximation quality is unknown until diagnostics are run. Most Bayesian-phylogenetics results in the wild are non-convergence reported as a result. Convergence is not a checkbox; it is the load-bearing claim. Without ESS, PSRF/ASDSF, and topological convergence, the output is a candidate, not a result.
3. Posterior probabilities are systematically higher (closer to 0 or 1) than bootstrap support and are overconfident under model misspecification (Suzuki 2002; Erixon 2003). PP=1.0 is routine and is not proof. Under the correct model PP is roughly calibrated but more decisive than the bootstrap; under an under-parameterized model PP becomes anti-conservative, assigning high support to wrong clades. This is structural: a Bayesian analysis conditions on the model being true, so a false model has no mechanism to express doubt. Support is not accuracy.

## When to Use Bayesian vs ML

| Factor | ML (modern-tree-inference) | MrBayes | BEAST2 | RevBayes | PhyloBayes-MPI |
|--------|----------------------------|---------|--------|----------|----------------|
| Reach for it when | topology, fast, IC model selection, bootstrap | posterior clade PP on a standard alignment; easiest entry; model averaging | TIME trees, tip-dating, demographics, phylodynamics | the model is not a built-in option; full graphical-model control | DEEP phylogeny with compositional heterogeneity / LBA |
| Output | point tree + bootstrap | posterior + clade PP | dated posterior + HPDs | posterior, custom | posterior under CAT/CAT-GTR |
| Speed | fast | moderate | slow | moderate | brutal (days-weeks, MPI) |
| Multi-run | automatic resampling | `nruns=2` default | run seeds independently | manual | run chains independently |
| Citation | Minh 2020 | Ronquist 2012 | Bouckaert 2019 | Hohna 2016 | Lartillot 2013 |

Default recommendation: start with ML (modern-tree-inference) for the topology, then move to Bayesian when posterior probabilities, model averaging, divergence times, or site-heterogeneous models are specifically needed. Need TIMES -> divergence-dating. Need a coalescent species tree -> species-trees.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| MrBayes 3.2 | Ronquist 2012 | NEXUS `lset`/`prset`; native MC3 + 2 runs; live ASDSF | general-purpose first choice for a concatenated alignment; partitioned data; `nst=mixed` model averaging |
| BEAST2 | Bouckaert 2019 | BEAUti -> XML; operator tuning; no MC3 by default | TIME trees, tip-dating, coalescent/birth-death demographics, phylodynamics (-> divergence-dating) |
| RevBayes | Hohna 2016 | Rev probabilistic-graphical-model scripting | a model that is not a checkbox: novel mixtures, custom hierarchical priors, biogeography |
| PhyloBayes-MPI | Lartillot 2013 (CAT: Lartillot 2004) | data-augmentation MCMC, MPI-parallel; CAT/CAT-GTR | deep phylogeny, compositional heterogeneity, LBA; `bpcomp`/`tracecomp` convergence |

Decision tree: standard alignment wanting clade PP, easiest path -> MrBayes (gammadir brlenspr, nruns>=2, gate on ASDSF<0.01). Need divergence TIMES -> BEAST2 (-> divergence-dating). Need a model that is not built-in -> RevBayes. Deep phylogeny with suspected compositional heterogeneity / LBA -> PhyloBayes CAT-GTR (or IQ-TREE C60/PMSF in modern-tree-inference if too big for CAT or fast bootstrap is wanted).

## MCMC Convergence -- the Heart

Almost all Bayesian-phylo failures live here, and the non-obvious core is that scalar-parameter convergence and tree-topology convergence are DIFFERENT problems. Continuous parameters (tree length, gamma shape, rates, frequencies) live in a smooth low-dimensional space that mixes well; topology lives in the vast discrete space of (2n-5)!! unrooted trees where a chain can get stuck on one island, mix perfectly within it, post a beautiful trace and ESS for every scalar, and never visit a competing topological peak. Diagnose BOTH; topology is the harder, more-often-skipped one.

- Question A (scalars converged?): ESS, autocorrelation, PSRF on the `.p`/`.log` traces. Tools: Tracer, coda, MrBayes `sump`.
- Question B (the TREE distribution converged?): between-run split-frequency agreement (ASDSF, max split-freq diff), topological ESS. Tools: MrBayes `sumt` (live ASDSF), PhyloBayes `bpcomp`, RWTY (Warren 2017). Passing A and stopping is the single most common silent error.

What each diagnostic measures:
- ESS (effective sample size) = N / integrated autocorrelation time. Check PER parameter (the likelihood can have ESS=2000 while a nuisance rate has ESS=12). Tracer flags <200 yellow, <100 red.
- PSRF (Gelman-Rubin Rhat) compares between-run to within-run variance across >=2 runs; PSRF -> 1.00 at convergence. Requires multiple runs. Scalar-only (Question A).
- ASDSF (average standard deviation of split frequencies) = THE between-run topology metric in MrBayes, computed automatically when `nruns>=2`: for each bipartition above a minimum frequency, take the SD of its frequency across runs and average. Converged runs give the same frequency in both -> ASDSF -> 0. Printed live during `mcmc`; `stoprule=yes stopval=0.01` auto-halts. Companion: the MAX split-freq diff (worst single split) catches one unsettled clade that a low average hides.
- Multiple independent runs are mandatory: PSRF and ASDSF are between-run metrics, undefined with one run, and a single run can look perfectly stationary on a local peak. A single long run is never sufficient evidence of convergence.
- MC3 (Metropolis-coupled MCMC) is MrBayes' defense against tree-space trapping: 1 cold chain (sampled) + heated chains whose flattened posterior crosses valleys, with periodic state swaps. Tuning signal is the chain-swap acceptance rate (target ~20-70%); if near 0, chains too far apart -> lower `temp`; near 100%, too similar -> raise `temp`. BEAST2 and PhyloBayes do not use MC3 by default.

Convergence checklist (all must pass before shipping):
1. Run >= 2 independent runs from different starting trees.
2. ESS > 200 for EVERY parameter (posterior, prior, likelihood, and all model params), not just the likelihood.
3. PSRF ~ 1.00 (<= 1.01) for all scalars (MrBayes `sump`).
4. Trace plots stationary and well-mixed ("fuzzy caterpillar"), independent runs overlay (Tracer).
5. Topology converged: ASDSF < 0.01 (MrBayes) / bpcomp maxdiff < 0.1 (PhyloBayes), and check RWTY tree-space MDS that runs overlap (not separate clouds).
6. Burn-in set from the trace / cumulative split-freq plots, not a dogmatic percentage. If any check fails, run longer; do NOT just increase thinning.

## Priors That Bite

Defaults are not neutral. The branch-length prior is the one that changes conclusions. The classic MrBayes default before 3.2.3 put an i.i.d. Exponential(10) prior on each branch with no control over the SUM, so as taxon number grows the implied prior on total tree length grows and the posterior is pulled toward implausibly long trees with degraded mixing -- and the inflation can feed back into TOPOLOGY and node support, so it is not cosmetic (Brown 2010). The fix is the compound (gamma-)Dirichlet prior (Zhang 2012): a diffuse Gamma on the whole tree length, partitioned among branches by a Dirichlet, decoupling "how long is the tree" from "how is length distributed". It yields posterior tree lengths close to the ML estimate and is robust to its hyperparameters; it became the MrBayes default at 3.2.3+. Rule: never trust a default `exp(10)` brlenspr on datasets with many taxa or long branches -- use `unconstrained:gammadir(...)`; inflated tree length vs an ML reference points first at the branch-length prior. Two further notes: a uniform topology prior is NOT uniform on clades; and Bayes-factor model comparison is only valid under PROPER priors (an unbounded "uninformative" prior leaves the marginal likelihood undefined).

The star-tree paradox and short-internode overconfidence: when the true internal branch is near zero (an effective polytomy), PP does NOT settle toward the uninformative 1/3 among the three resolutions as data accumulate -- it behaves erratically and can drive toward HIGH support for an arbitrary resolution (Lewis 2005; the theory in Yang 2007). So a high PP on a clade subtended by a very short internal branch is exactly where Bayesian support is least trustworthy. Treat near-zero internodes as soft, cross-check against bootstrap, and consider a polytomy-allowing reversible-jump prior (Lewis 2005).

## Model Comparison Done Right

Bayesian model comparison compares MARGINAL likelihoods (the data probability integrated over all parameters under the model's priors) via Bayes factors. The marginal likelihood is a hard high-dimensional integral. The correct estimators are stepping-stone (Xie 2011) and path sampling, which sample a series of power posteriors interpolating prior (beta=0) to posterior (beta=1); stepping-stone is more accurate per step and is the default recommendation. Baele 2012 showed by simulation and empirically that both substantially outperform the harmonic-mean estimator.

The harmonic-mean estimator (HME) is discredited and must never select a model. It is dominated by the smallest likelihoods in the posterior sample (the prior-favored tail the posterior rarely visits), giving it effectively infinite variance -- it does not converge as samples are added, is unstable run-to-run, systematically overstates the marginal likelihood, and favors over-parameterized models. Interpret Bayes factors on the 2*ln(BF) scale (Kass & Raftery 1995): 2-6 positive, 6-10 strong, >10 very strong/decisive, where 2 lnBF = 2*(lnML_1 - lnML_2). BFs require proper priors. The reversible-jump alternative `lset nst=mixed` sidesteps explicit BFs by sampling the substitution model itself (the 203 GTR rate-class groupings) and reporting each model's posterior probability -- the principled way to account for substitution-model uncertainty by averaging.

### Run MrBayes and Verify Convergence Before Trusting the Tree

**Goal:** Produce a Bayesian phylogeny whose clade posteriors are trustworthy, gated on proven convergence rather than a single stationary-looking run.

**Approach:** Set a multiple-hit-correcting model and the compound-Dirichlet branch-length prior, run two MC3 runs, watch live ASDSF, then summarize only after the topology and scalar diagnostics pass.

```
begin mrbayes;
    lset nst=6 rates=invgamma;                          [ GTR+I+G; nst=mixed = rjMCMC model averaging ]
    prset brlenspr=unconstrained:gammadir(1,0.1,1,1);   [ compound Dirichlet, NOT exp(10): avoids tree-length inflation ]
    mcmc ngen=10000000 nruns=2 nchains=4 temp=0.1       [ 2 runs x (1 cold + 3 heated MC3 chains) ]
         samplefreq=1000 printfreq=1000 diagnfreq=5000
         stoprule=yes stopval=0.01;                     [ auto-halt when ASDSF < 0.01 (topology converged) ]
    sump burninfrac=0.25 relburnin=yes;                 [ scalar PSRF + ESS, discarding first 25% ]
    sumt burninfrac=0.25 relburnin=yes;                 [ consensus tree + clade PP + ASDSF ]
end;
```

After the run: confirm `sump` PSRF ~ 1.00 and ESS > 200 for every parameter, `sumt` ASDSF < 0.01 with a small max split-freq diff, and check tree-space convergence with RWTY (`analyze.rwty(list(run1=..., run2=...), burnin=25)`, then `makeplot.treespace`). Distrust any clade where PP is high but the bootstrap or concordance factor (modern-tree-inference) is low.

### Compare Two Models by Stepping-Stone

**Goal:** Select between two models (e.g. GTR+G vs GTR+I+G, or partition schemes) by marginal likelihood, never by the harmonic mean.

**Approach:** Estimate each model's marginal log-likelihood by stepping-stone sampling, then convert the difference to a Kass-Raftery Bayes factor.

```
[ run once per model on the same proper priors; nsteps = power-posterior stones ]
ss ngen=1000000 nsteps=50 diagnfreq=1000;
[ MrBayes prints the stepping-stone marginal log-likelihood; record it per model ]
[ 2 lnBF = 2 * (lnML_model1 - lnML_model2); interpret on the Kass-Raftery scale ]
```

BEAST2: install the `MODEL_SELECTION` package and run `PathSampler` (path sampling / stepping-stone); RevBayes: `powerPosterior()` + `steppingStoneSampler()`. All require proper priors, or the marginal likelihood is undefined.

## Site-Heterogeneous CAT Models for Deep Phylogeny

Standard models (GTR+G, LG/WAG+G) assume all sites share one set of equilibrium frequencies -- they are site-HOMOGENEOUS. Real proteins are not: a buried hydrophobic site and a surface charged site have different profiles. At DEEP timescales this misleads, because saturated sites convergently acquire similar compositions in unrelated lineages and a homogeneous model misreads convergent composition as shared ancestry -> long-branch attraction with HIGH PP (the overconfidence-under-misspecification mechanism made concrete). The answer is CAT / CAT-GTR (Lartillot 2004): an infinite-mixture (Dirichlet-process) model assigning sites to an unknown number of categories, each with its own amino-acid frequency profile, all inferred from the data; CAT-GTR adds one shared GTR exchangeability matrix (the PhyloBayes default with 4-category gamma). It is far more robust to LBA at depth than homogeneous models, at the cost of slow convergence (millions of cycles, days-weeks), which is why PhyloBayes-MPI parallelizes it (Lartillot 2013).

CAT chains converge slowly, so diagnostics are mandatory (run >= 2 chains with different names): `bpcomp -x <burnin> <every> chain1 chain2` reports maxdiff (the largest bipartition-frequency discrepancy across chains, the topology metric) and writes a pooled consensus; `tracecomp -x <burnin> chain1 chain2` reports per-statistic effsize and rel_diff (a PSRF analog). Thresholds: maxdiff < 0.1 AND effsize > 300 = converged; maxdiff < 0.3 AND effsize > 50 = acceptable; otherwise keep running. The ML treatment of the same artifact is IQ-TREE C60/PMSF (Wang 2018): fit the C60 profile mixture once, fix each site's posterior-mean profile, then do fast tree search and bootstrap. CAT-GTR gives a posterior (slow, MPI); PMSF gives bootstrap (fast, scales to huge matrices) -- different output, same disease. Routes to modern-tree-inference.

```bash
mpirun -np 8 pb_mpi -d alignment.phy -cat -gtr -dgam 4 chain1   # CAT-GTR + discrete gamma
mpirun -np 8 pb_mpi -d alignment.phy -cat -gtr -dgam 4 chain2   # second independent chain
bpcomp -x 1000 10 chain1 chain2     # maxdiff (topology); writes bpcomp.con.tre
tracecomp -x 1000 chain1 chain2     # effsize + rel_diff (scalars)
```

## Per-Method Failure Modes

### Non-Convergence Reported as a Result
**Trigger:** A single short chain, a stationary-looking trace, a point estimate shipped with no ESS / PSRF / ASDSF and no second run.
**Mechanism:** The "result" is an arbitrary draw from a chain that never demonstrably sampled the posterior; the clade PPs are an artifact of where the chain happened to sit.
**Symptom:** No diagnostics reported; or low ESS / PSRF > 1.05 / ASDSF > 0.01 quietly ignored.
**Fix:** Always >= 2 runs; gate on ESS > 200 (all params) AND ASDSF < 0.01 / maxdiff < 0.1; inspect in Tracer/RWTY; never ship without the diagnostics.

### Scalar Convergence Mistaken for Topology Convergence
**Trigger:** Every scalar has ESS > 1000 and PSRF = 1.00, so convergence is declared.
**Mechanism:** The chain is trapped on one tree-space island, mixing perfectly within it while never sampling a competing topology; scalars are necessary but not sufficient.
**Symptom:** Perfect Tracer traces but ASDSF stuck high, or RWTY tree-space MDS shows two separate clouds for the two runs.
**Fix:** Always run a TOPOLOGY diagnostic (ASDSF / bpcomp maxdiff / RWTY MDS), not just Tracer.

### Branch-Length Prior Inflation (Prior Domination)
**Trigger:** Default `exp(10)` brlen prior (old MrBayes) on many taxa or long branches, or a too-strong tree/clock prior on weak data.
**Mechanism:** The i.i.d.-per-branch prior places mass on long total tree length; the posterior converges to the WRONG distribution, and inflated branch lengths can distort topology and support.
**Symptom:** Posterior tree length far exceeds the ML estimate; topology shifts when the prior is changed.
**Fix:** Use `unconstrained:gammadir(...)`; do a prior-sensitivity check (sample the prior alone, vary hyperparameters); compare tree length to an ML reference.

### Harmonic-Mean Model Selection
**Trigger:** Selecting a substitution/clock/partition model by HME (or AICM) Bayes factors.
**Mechanism:** The HME has effectively infinite variance, is unstable, and favors over-parameterized models -> the wrong model is chosen.
**Symptom:** Bayes factors that change on rerun, or a consistent preference for the most complex model.
**Fix:** Stepping-stone or path sampling for marginal likelihoods; interpret with Kass-Raftery 2lnBF thresholds; ensure proper priors.

### PP Overconfidence Under Misspecification
**Trigger:** Reporting PP = 1.0 as proof, especially site-heterogeneous data analyzed under a homogeneous model.
**Mechanism:** A Bayesian analysis conditions on the model; a misspecified model has no way to hedge, so the posterior concentrates with false confidence on LBA-driven clades.
**Symptom:** PP = 1.0 where the bootstrap is 55%, or two conflicting topologies each with high PP across analyses.
**Fix:** Report the model; cross-check PP against bootstrap (sharp disagreement = warning); use CAT-GTR / PMSF for deep data; treat PP = 1.0 as "decisive under THIS model", not truth.

### CAT Non-Convergence at Scale
**Trigger:** Stopping a PhyloBayes CAT run too early because wall-clock budget ran out.
**Mechanism:** CAT chains converge slowly (millions of cycles); maxdiff is still ~0.3 and effsize < 50 when the tree is read out.
**Symptom:** bpcomp maxdiff well above 0.1, tracecomp effsize below 300.
**Fix:** Run long, >= 2 chains, gate on maxdiff < 0.1 AND effsize > 300; budget days-weeks of MPI time; consider removing the fastest sites or recoding.

## Quantitative Thresholds

| Quantity | Threshold | Source |
|----------|-----------|--------|
| ESS per parameter | > 200 (floor 100; 500-1000 for reported point estimates) | Rambaut 2018 |
| PSRF / Rhat | ~ 1.00, accept <= 1.01; > 1.05 = not converged | Ronquist 2012 |
| ASDSF (MrBayes between-run topology) | < 0.01 for a confident result; < 0.05 adequate for hard data | Ronquist 2012 |
| Max split-freq diff | small, consistent with ASDSF < 0.01 | Ronquist 2012 |
| bpcomp maxdiff (PhyloBayes topology) | < 0.1 good; < 0.3 acceptable | Lartillot 2013 |
| tracecomp effsize (PhyloBayes scalars) | > 300 good; > 50 acceptable; rel_diff small (<~0.1) | Lartillot 2013 |
| Burn-in | discard 10-25% (data-driven; up to 50% for slow chains) | Ronquist 2012 |
| MC3 chain-swap acceptance | ~20-70% between adjacent chains; tune `temp` if outside | Ronquist 2012 |
| Bayes factor (2 lnBF) | 2-6 positive, 6-10 strong, > 10 very strong/decisive | Kass & Raftery 1995 |
| Independent runs | minimum 2 (4 for difficult data) | Ronquist 2012 |
| Posterior probability | >= 0.95 nominal but overconfident under misspecification; PP = 1.0 != proof | Suzuki 2002 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Convergence declared on one run | PSRF/ASDSF undefined with a single run | run nruns >= 2; gate on between-run metrics |
| Perfect ESS but wrong PPs | tree-space trapping not caught by scalars | add a topology diagnostic (ASDSF / bpcomp / RWTY MDS) |
| Posterior tree length >> ML | default `exp(10)` branch-length prior | switch to `unconstrained:gammadir(...)` |
| Model choice flips on rerun | harmonic-mean estimator used | use stepping-stone / path sampling |
| Bayes factor is undefined / arbitrary | an improper "uninformative" prior | make all priors proper before SS/PS |
| PP = 1.0 on an LBA clade at depth | site-heterogeneous data under a homogeneous model | use PhyloBayes CAT-GTR or IQ-TREE PMSF |
| High PP on a near-zero internode | star-tree paradox fabricates resolution | treat the internode as soft; consider a polytomy prior |
| Low ESS after thinning more | thinning discards information, does not raise ESS | run more generations, not higher samplefreq |

## References

Ronquist F, Teslenko M, van der Mark P, Ayres DL, Darling A, Hohna S, Larget B, Liu L, Suchard MA, Huelsenbeck JP. 2012. MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. *Systematic Biology* 61(3):539-542.
Bouckaert R, Vaughan TG, Barido-Sottani J, Duchene S, Fourment M, et al. 2019. BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. *PLoS Computational Biology* 15(4):e1006650.
Hohna S, Landis MJ, Heath TA, Boussau B, Lartillot N, Moore BR, Huelsenbeck JP, Ronquist F. 2016. RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language. *Systematic Biology* 65(4):726-736.
Lartillot N, Philippe H. 2004. A Bayesian mixture model for across-site heterogeneities in the amino-acid replacement process. *Molecular Biology and Evolution* 21(6):1095-1109.
Lartillot N, Rodrigue N, Stubbs D, Richer J. 2013. PhyloBayes MPI: phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. *Systematic Biology* 62(4):611-615.
Wang HC, Minh BQ, Susko E, Roger AJ. 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. *Systematic Biology* 67(2):216-235.
Warren DL, Geneva AJ, Lanfear R. 2017. RWTY (R We There Yet): an R package for examining convergence of Bayesian phylogenetic analyses. *Molecular Biology and Evolution* 34(4):1016-1020.
Rambaut A, Drummond AJ, Xie D, Baele G, Suchard MA. 2018. Posterior summarization in Bayesian phylogenetics using Tracer 1.7. *Systematic Biology* 67(5):901-904.
Xie W, Lewis PO, Fan Y, Kuo L, Chen MH. 2011. Improving marginal likelihood estimation for Bayesian phylogenetic model selection. *Systematic Biology* 60(2):150-160.
Baele G, Lemey P, Bedford T, Rambaut A, Suchard MA, Alekseyenko AV. 2012. Improving the accuracy of demographic and molecular clock model comparison while accommodating phylogenetic uncertainty. *Molecular Biology and Evolution* 29(9):2157-2167.
Kass RE, Raftery AE. 1995. Bayes factors. *Journal of the American Statistical Association* 90(430):773-795.
Brown JM, Hedtke SM, Lemmon AR, Lemmon EM. 2010. When trees grow too long: investigating the causes of highly inaccurate Bayesian branch-length estimates. *Systematic Biology* 59(2):145-161.
Zhang C, Rannala B, Yang Z. 2012. Robustness of compound Dirichlet priors for Bayesian inference of branch lengths. *Systematic Biology* 61(5):779-784.
Suzuki Y, Glazko GV, Nei M. 2002. Overcredibility of molecular phylogenies obtained by Bayesian phylogenetics. *PNAS* 99(25):16138-16143.
Erixon P, Svennblad B, Britton T, Oxelman B. 2003. Reliability of Bayesian posterior probabilities and bootstrap frequencies in phylogenetics. *Systematic Biology* 52(5):665-673.
Lewis PO, Holder MT, Holsinger KE. 2005. Polytomies and Bayesian phylogenetic inference. *Systematic Biology* 54(2):241-253.
Yang Z. 2007. Fair-balance paradox, star-tree paradox, and Bayesian phylogenetics. *Molecular Biology and Evolution* 24(8):1639-1655.

## Related Skills

- modern-tree-inference - ML inference, bootstrap/SH-aLRT support, concordance factors, and IQ-TREE C60/PMSF
- divergence-dating - BEAST2 clock models, calibration, and time-scaled posteriors
- species-trees - coalescent species-tree estimation when ILS dominates
- tree-io - reading and summarizing MrBayes/BEAST2 trees without dropping posteriors and HPDs
