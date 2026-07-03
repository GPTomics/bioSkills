---
name: bio-epidemiological-genomics-phylodynamics
description: Estimates time-scaled phylogenies, molecular-clock rates, effective reproduction number R_e, and population dynamics from dated pathogen genomes using TreeTime (maximum-likelihood) and BEAST2 (Bayesian; strict/relaxed clocks; coalescent, Bayesian-Skyline, Skygrid, Birth-Death-Skyline, and sampled-ancestor priors; structured coalescent via MASCOT). Covers root-to-tip clock QC via TempEst, date-randomisation tests, recombination masking via Gubbins/ClonalFrameML before clock inference for recombining bacteria, BDSKY origin-vs-rootHeight pitfalls, sampling-bias correction, multi-chain convergence diagnostics, and reconciling phylodynamic R_e with case-based R_t. Use when dating outbreak origins, estimating substitution rates, inferring R_e through time, building time-calibrated Nextstrain Augur trees, choosing between strict and relaxed clocks, fitting Birth-Death-Skyline models, diagnosing temporal-signal failure, running MASCOT for structured-population analyses, or using UShER for pandemic-scale placement.
tool_type: mixed
primary_tool: BEAST2
---

## Version Compatibility

Reference examples tested with: BEAST 2.7.6+, BDSKY 1.5+, BEASTLabs 2.0+, feast 9.5+, ORC 1.1.2+, MASCOT 3.0+, BEAGLE 4.0+, TreeTime 0.11+, IQ-TREE 2.3.6+, Gubbins 3.3+, ClonalFrameML 1.13+, UShER 0.6+, matUtils 0.6+, BactDating 1.1+ (R), bdskytools 1.1+ (R), coda 0.19+ (R), ape 5.8+ (R), ggplot2 3.5+, BioPython 1.84+, dendropy 4.6+, baltic 0.2+.

Before using code patterns, verify installed versions match. If versions differ:
- BEAST: `beast -version`; `packagemanager -list` for BEAST2 packages and versions
- Python: `pip show treetime`; `help(treetime.TreeTime)`
- R: `packageVersion('bdskytools')`; `?bdskytools::bdskytools_plot`
- CLI: `gubbins --version`; `iqtree --version`; `clonalframeml --version`

If BEAST2 throws `IllegalArgumentException` on XML load, the BDSKY / feast / BEASTLabs minor version probably moved; check the XML against the installed package's example XML in `~/beast/examples/`. BEAST2 XML is NOT robust across minor releases; pin the BEAST package version in any published analysis.

# Phylodynamics

**"How fast is this outbreak growing, and when did it start?"** -> Combine a dated set of pathogen genomes with a molecular-clock model to time-scale the phylogeny, then fit a population-dynamic model (constant / exponential / Bayesian Skyline / BICEPS / Birth-Death-Skyline) to read off R_e, growth rate, and origin date. Choice of clock (strict vs UCLN vs ORC) and tree prior (BSP coalescent vs BICEPS vs BDSKY birth-death) is load-bearing; the same data can yield different R_e under different priors. For bacterial pathogens, recombination MUST be masked first; running BEAST on a *Streptococcus pneumoniae* or *E. coli* core-genome alignment without Gubbins / ClonalFrameML inflates the clock rate 2-5x and the date-randomisation test is NOT a sufficient guard.

- CLI: `treetime --tree raw.nwk --aln aln.fasta --dates dates.tsv --coalescent skyline --clock-filter 4` -- fast ML phylodynamics
- Java/CLI: BEAST2 with BDSKY XML (BEAUti-generated) -- full Bayesian birth-death-skyline with R_e per epoch
- CLI: `gubbins --prefix gubbins core.full.aln` -- recombination masking before bacterial clock inference
- R: `bdskytools::bdskytools_plot` for BDSKY post-processing; `BactDating` for fast Bayesian dating after Gubbins

## The Single Most Important Modern Insight -- Recombination passed unmasked into clock inference inflates the clock rate 2-5x and breaks every downstream estimate

The date-randomisation test is NOT a guard against unmasked recombination. For any recombining bacterium (*S. pneumoniae*, *N. gonorrhoeae*, *E. coli*, *Klebsiella pneumoniae*, *Campylobacter*, *Helicobacter pylori*), run Gubbins (Croucher 2015 *NAR* 43:e15) or ClonalFrameML (Didelot & Wilson 2015 *PLoS Comput Biol* 11:e1004041) on the core-SNP alignment FIRST, rebuild the tree on the recombination-masked alignment, THEN run TreeTime / BEAST. Only *M. tuberculosis* and a handful of clonal pathogens are exempt -- and even those benefit from a recombination check on cross-lineage analyses. The Mostowy 2017 *Mol Biol Evol* 34:1167 fastGEAR paper documents the limits of mask-based approaches for highly recombinogenic species; for those (*N. gonorrhoeae*, *S. pneumoniae* lineages with strong recombination), residual signal persists post-masking and biases downstream R_e estimates downward. Second-order insight: Volz & Frost 2014 *J R Soc Interface* 11:20140945 showed that BEAST coalescent priors are biased under realistic preferential sampling; BDSKY models the sampling proportion explicitly and is the correct tool when sampling rate varies; MASCOT-Skyline / MASCOT-GLM (Müller 2018 *Bioinformatics* 34:3843) further correct for sampling-deme covariation. Third-order insight: Featherstone & Duchêne 2023 *Mol Biol Evol* 40:msad132 quantified that for shallow trees with many samples, sampling times dominate over sequence information for R_e inference -- biased sampling drives biased R_e estimates regardless of how much sequence data is added.

## Algorithmic Taxonomy

| Tool / model | Mechanism | Outputs | Strength | Fails when |
|--------------|-----------|---------|----------|------------|
| TreeTime ML (Sagulenko 2018 *Virus Evol* 4:vex042) | ML joint optimisation of clock + dates with optional coalescent skyline prior | Time-scaled tree + clock rate + root-to-tip regression | 100-1000x faster than BEAST; ideal for outbreak-scale data | Strict-clock assumption; no posterior; no R_e directly |
| BEAST2 + BICEPS (Bouckaert 2022 *Syst Biol* 71:1549) | Bayesian skyline with analytic Ne integration per epoch and new tree-flexing operators | Ne(t) posterior | Modern default skyline; weeks-of-BSP becomes hours | Replaces BSP for many use cases; check current BEAST2 tutorials |
| BEAST2 + BDSKY (Stadler 2013 *PNAS* 110:228) | Birth-death-skyline with explicit sampling | R_e(t), become-uninfectious rate, sampling proportion | Direct R_e estimation | `origin` vs `rootHeight` confusion; rho-and-turnover unidentifiability with flat priors (Legried & Terhorst 2022 *PNAS* 119:e2119513119) |
| BEAST2 + MASCOT (Müller 2018 *Bioinformatics* 34:3843) | Marginal-approximation structured coalescent | Per-deme Ne + migration | Correct for structured sampling; replaces biased DTA (De Maio 2015 *PLoS Genet* 11:e1005421) | Migration unidentifiable with <20 sequences per deme |
| BEAST2 + MASCOT-Skyline / MASCOT-GLM | Time-varying migration with covariates | Migration-rate trajectories tied to predictors | Sampling-aware; addresses Volz & Frost 2014 sampling bias | ~10x slower than DTA; many users still default to DTA |
| BEAST2 + Sampled-Ancestor BDSKY | BDSKY with internal-node sampling | R_e + ancestral / longitudinal samples | Right for ancient-DNA, within-host longitudinal sampling | Specialist parameterisation |
| BactDating (Didelot 2018 *NAR* 46:e134) | Bayesian Poisson / mixedgamma / relaxedgamma clock on a fixed tree | Time-scaled tree + clock rate posterior | Fast Bayesian dating after Gubbins; right for large bacterial trees | `mixedgamma` mixes poorly; `poisson` is the cleaner default |
| Gubbins (Croucher 2015 *NAR* 43:e15) | Sliding-window elevated SNP density detection | Recombination-masked alignment + recombination GFF | Standard for clonal bacterial alignments | Cannot detect ancient recombination; mis-masks mutation hotspots |
| ClonalFrameML (Didelot & Wilson 2015 *PLoS Comput Biol* 11:e1004041) | Coalescent-with-recombination model on a fixed tree | Recombination-masked alignment + r/m | Model-based alternative to Gubbins | Slow on large trees |
| UShER + matUtils (Turakhia 2021 *Nat Genet* 53:809) | Parsimony placement on a daily-updated mutation-annotated tree | Subtrees, lineage assignments, RIPPLES recombination calls | Pandemic-scale (millions of genomes) | Parsimony branch lengths systematically shorter than ML; re-estimate branch lengths before downstream R_e |
| TempEst (Rambaut 2016 *Virus Evol* 2:vew007) | Root-to-tip linear regression | Clock signal R^2 | First-line temporal-signal diagnostic | Slope can be artificially good with biased sampling |
| Date randomisation (Ramsden 2009 *Mol Biol Evol* 26:143; Duchêne 2015 *Mol Biol Evol* 32:1895) | Shuffle dates, compare clock-rate estimate | Pass / fail | Detects spurious clock signal | Can "pass" with narrow sampling windows (false negative) |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| "Estimate Ne(t) for this virus" | BEAST2 + BICEPS (or Skygrid if BICEPS unavailable) | Constant coalescent without checking flatness; BSP if BICEPS available (BSP mixes poorly) |
| "Estimate R_e from sequences" | BEAST2 + BDSKY with sampling-process explicit; document sampling proportion per epoch | BSP-style Ne -> R_e conversion via Wallinga-Lipsitch loses sampling information |
| "Multi-deme analysis with migration" | BEAST2 + MASCOT for <=10 demes with >=20 sequences each | Exact structured coalescent (intractable >5 demes); BEAST DTA (Lemey 2009) is sampling-biased |
| "Pandemic-scale (>10k sequences)" | UShER + matUtils for placement; TreeTime for dates; BDSKY on lineage-specific subsets | Full BEAST on full dataset is intractable; using UShER branch lengths directly biases R_e |
| "Date a bacterial tree" | Snippy -> Gubbins -> IQ-TREE -> BactDating OR BEAST2; recombination mask FIRST | Skipping recombination masking inflates the clock rate 2-5x |
| "Date a fast-evolving virus" | TreeTime (Nextstrain pipeline) for routine; BEAST2 + UCLN for headline analyses | Strict clock by default underestimates rate variation on shallow trees |
| "Test for temporal signal" | TempEst root-to-tip first (R^2 >= 0.3 minimum); date-randomisation as secondary check | Skipping the diagnostic; trusting date-randomisation alone (can pass with narrow window) |
| "Reconcile phylodynamic R_e with case R_t" | Report both with explicit assumptions; investigate disagreement (sampling bias, lineage-specific signal) | Reporting one as "the" R_e |
| "Bacterial outbreak phylogenetics start-to-finish" | snippy + snippy-core -> Gubbins (on `core.full.aln`) -> IQ-TREE -> BactDating or BEAST2 + BDSKY | Skipping Gubbins; running Gubbins on `core.aln` instead of `core.full.aln` |
| "Migration / phylogeography source-attribution" | MASCOT or MASCOT-GLM (sampling-aware); never BEAST DTA for attribution claims | BEAST DTA inherits Lemey 2009 sampling bias; published source-attribution remains biased toward heavily-sampled locations |

Methodology evolves; before any high-stakes phylodynamic analysis, web-search "BEAST2 BDSKY tutorial 2025" and "MASCOT-Skyline benchmark" for current best practice.

## Time-Scaling With TreeTime

**Goal:** Produce a time-scaled phylogeny with per-node date estimates and a global clock rate, ready for downstream R_e estimation or visualisation -- in minutes rather than hours.

**Approach:** Build the raw topology with IQ-TREE 2 (Minh 2020 *Mol Biol Evol* 37:1530) for outbreak-scale data (RAxML-NG for larger trees); pass to TreeTime jointly optimising the molecular clock and date assignments with a coalescent skyline prior; inspect `root_to_tip_regression.pdf` BEFORE trusting any downstream output; apply `--clock-filter 4` to drop tips with root-to-tip residuals exceeding 4 SDs.

```bash
iqtree -s aln.fasta -m GTR+G -B 1000 -T AUTO -pre raw_tree

treetime \
    --tree raw_tree.treefile \
    --aln aln.fasta \
    --dates dates.tsv \
    --coalescent skyline \
    --clock-filter 4 \
    --confidence \
    --reroot best \
    --outdir timetree

treetime clock \
    --tree raw_tree.treefile \
    --dates dates.tsv \
    --reassign-dates \
    --outdir date_randomisation
```

Outputs: `timetree.nexus` is the time-scaled tree; `dates.tsv` records per-tip filter status; `root_to_tip_regression.pdf` is the temporal-signal diagnostic. If R^2 < 0.3, the data do not support time-scaled inference; report uncertainty and consider extending the sampling window. For published clock rates, run the date-randomisation analysis (Duchêne 2015 *Mol Biol Evol* 32:1895) and report the clock-rate distribution under shuffled dates.

## BDSKY in BEAST2 With Recombination Masking First

**Goal:** Estimate R_e through time from a bacterial outbreak alignment with explicit handling of recombination, sampling proportion, and convergence diagnostics.

**Approach:** Snippy + snippy-core to build the core-genome alignment; Gubbins on `core.full.aln` (full positions including invariant) to mask recombinant tracts; IQ-TREE on the masked alignment; BEAUti to set up BDSKY XML in BEAST 2 (Bouckaert 2019 *PLoS Comput Biol* 15:e1006650 for BEAST 2.5+) with origin -- NOT rootHeight; sampling proportion per epoch; become-uninfectious rate fixed from epi knowledge; run with 3-4 independent chains from different seeds; combine chains only after marginal posteriors overlap.

```bash
snippy-core --ref reference.fa --prefix core snippy_out/*

run_gubbins.py --prefix gubbins core.full.aln

iqtree -s gubbins.filtered_polymorphic_sites.fasta -m GTR+G+ASC -B 1000 -T AUTO -pre masked_tree

beast -threads 4 -beagle -seed 42 bdsky_analysis.xml
beast -threads 4 -beagle -seed 17 bdsky_analysis.xml
beast -threads 4 -beagle -seed 99 bdsky_analysis.xml
beast -threads 4 -beagle -seed 7 bdsky_analysis.xml

logcombiner -log run_42.log -log run_17.log -log run_99.log -log run_7.log -burnin 10 -o combined.log
logcombiner -log run_42.trees -log run_17.trees -log run_99.trees -log run_7.trees -burnin 10 -o combined.trees -decimalPlaces 6
loganalyser -burnin 10 combined.log
```

`run_gubbins.py` input MUST be `core.full.aln` (full alignment with reference). Passing `core.aln` (variable-only) gives wrong recombination calls because Gubbins cannot estimate background SNP density without invariant positions. IQ-TREE `+ASC` is the ascertainment-bias correction required for SNP-only input.

## MASCOT Structured Coalescent

**Goal:** Infer per-deme effective population size and migration rates from sequences sampled in multiple subpopulations (countries, hospitals, ward types) without inheriting the Lemey 2009 BEAST DTA sampling bias.

**Approach:** BEAUti -> MASCOT template; one trait per deme; require >=20 sequences per deme for migration identifiability; consider MASCOT-GLM if migration rates plausibly depend on observable covariates (travel volume, geographic distance); pre-2024 default is MASCOT, but MASCOT-Skyline / MASCOT-GLM is preferred when sampling intensity varies over time.

```bash
beast -threads 4 -beagle -seed 42 mascot_analysis.xml
loganalyser -burnin 10 mascot.log
```

## Per-Method Failure Modes

### Recombination passed unmasked into clock inference

**Trigger:** BEAST2 or TreeTime run on a core-genome alignment of a recombining bacterium (S. pneumoniae, N. gonorrhoeae, E. coli, Klebsiella, Campylobacter, Helicobacter pylori).

**Mechanism:** Recombination imports SNPs from a divergent donor lineage in a single event. The clock model interprets these as accumulated point mutations across the branch, inflating the apparent clock rate and distorting node-date estimates. Recombination is non-clocklike, so date-randomisation tests may still pass -- the problem is silent.

**Symptom:** Estimated clock rate is 2-5x the literature consensus for the species (S. pneumoniae core clock ~1.5e-6 subs/site/year per literature; rates >5e-6 indicate unmasked recombination). Per-branch dN/dS profile is wildly heterogeneous. Some branches show implausibly recent divergence dates.

**Fix:** Mask recombinant regions before clock inference. Build initial tree with IQ-TREE on the core-SNP alignment; run Gubbins or ClonalFrameML to detect recombinant tracts; rebuild tree on the recombination-masked alignment; THEN run TreeTime or BEAST. For *M. tuberculosis* (rare recombination), masking is optional but defensible for cross-lineage analyses.

### Date-randomisation test passes despite no real temporal signal

**Trigger:** Outbreak with narrow sampling window (e.g. all isolates collected within 3 months of a year-long outbreak).

**Mechanism:** Date randomisation (Ramsden 2009 *Mol Biol Evol* 26:143; Duchêne 2015 *Mol Biol Evol* 32:1895) tests whether shuffling dates degrades the clock-rate estimate. With insufficient temporal sampling, the true clock estimate is also poorly informed, so randomised and true estimates overlap by chance.

**Symptom:** TempEst root-to-tip regression R^2 < 0.1; date-randomisation test still "passes" (HPD overlap); 95% HPD on clock rate spans an order of magnitude.

**Fix:** Inspect root-to-tip regression FIRST (TempEst). If R^2 < 0.3 and there is no strong a-priori clock-rate prior from the literature, the data do not support time-scaled inference. Options: (1) use a strong informative prior on clock rate from outside data; (2) extend the sampling window before re-running; (3) report a tree without time-scaling and discuss uncertainty.

### BDSKY origin specified as the root height

**Trigger:** BEAST2 BDSKY XML where `origin` is set to the same value as `rootHeight` (often because the user inferred from the tutorial that "origin = tree depth").

**Mechanism:** Stadler 2013 *PNAS* 110:228 defines `origin` as the time from the start of the epidemic to the most recent sample -- strictly larger than the tree root height (tMRCA). Setting `origin = tMRCA` causes the MCMC to start in an inconsistent state and systematically biases R_e estimates upward (because turnover is forced into a shorter time window).

**Symptom:** R_e estimates implausibly high in early epochs; chains mix poorly; origin-date estimate clusters at the lower bound of the prior.

**Fix:** Initialise `origin` to (tMRCA + 0.1*tMRCA) or use prior knowledge (e.g., for SARS-CoV-2 within a country, the documented import date). BEAUti default is sensible; hand-edited XML often gets this wrong.

### BSP / BDSKY ESS < 200 reported as a result

**Trigger:** Single BEAST chain reaching the planned MCMC length; some parameters with ESS in the 50-150 range; user reports the posterior anyway.

**Mechanism:** Tracer's ESS calculation is a single-chain effective sample size. For phylogenetic posteriors, parameters may be "mixed within chain" but "unmixed across chains" -- multiple independent chains can converge to different parts of the posterior. ESS > 200 is necessary but not sufficient.

**Symptom:** Reported HPD intervals from a single chain; reviewers from Stadler / Bouckaert / Suchard schools reject the analysis.

**Fix:** Run >=3-4 chains from different starting trees / seeds; examine marginal posteriors per chain; only after they overlap can chains be combined (`logcombiner`). Report Gelman-Rubin diagnostic via R `coda::gelman.diag` on the per-chain log files.

### MASCOT migration with too few sequences per deme

**Trigger:** MASCOT analysis with <20 sequences per deme; user reports migration-rate posterior.

**Mechanism:** MASCOT migration rates are jointly identifiable only with sufficient sequences per deme to inform within-deme coalescent. With few sequences, migration rate and Ne become confounded; the posterior reflects the prior more than the data.

**Symptom:** Migration-rate 95% HPD spans 2+ orders of magnitude; estimates implausibly extreme.

**Fix:** Pool nearby demes; accept wide HPD intervals; report MASCOT-GLM if covariates are available; cite Müller 2018 for the identifiability requirement.

### DTA used for phylogeography source-attribution

**Trigger:** BEAST DTA (Lemey 2009 *PLoS Comput Biol* 5:e1000520) used to claim a geographic source for an outbreak.

**Mechanism:** De Maio 2015 *PLoS Genet* 11:e1005421 demonstrated DTA is biased toward heavily-sampled locations; the example was Ebola DTA implausibly concluding humans seeded outbreaks (truth: sylvatic reservoir spillover). DTA treats sampling as random; in reality, source-attribution-relevant locations are often the most under-sampled.

**Symptom:** Inferred "source" is the heavily-sampled location; conclusion is sensitive to sub-sampling; reviewers familiar with MASCOT push back.

**Fix:** Use MASCOT or MASCOT-Skyline / MASCOT-GLM for source-attribution claims; if infeasible, frame results as "consistent with" rather than "demonstrates" and quantify sampling bias.

### UShER branch lengths used directly for R_e estimation

**Trigger:** Downstream BDSKY analysis on a UShER-placed subtree using UShER's parsimony branch lengths.

**Mechanism:** Turakhia 2021 *Nat Genet* 53:809 UShER places sequences on the MAT via parsimony; branch lengths under parsimony are systematically shorter than maximum-likelihood or Bayesian branch lengths (because parsimony minimises substitutions). R_e estimates that use UShER branch lengths directly are biased.

**Symptom:** Implausibly high R_e from a UShER-derived subtree; estimates inconsistent with case-based R_t.

**Fix:** Use UShER for placement only; re-estimate branch lengths with TreeTime or BEAST before downstream R_e estimation. This caveat is in the UShER documentation but routinely ignored.

### R_e reported as R_0

**Trigger:** Phylodynamic R_e estimate (current immunity / intervention context) described in text as R_0 (fully susceptible population).

**Mechanism:** Phylodynamic methods estimate the *effective* reproduction number R_e (or R_t), not the basic reproduction number R_0. The two diverge substantially: R_0 for ancestral SARS-CoV-2 was ~5-8; R_e during Omicron waves was 1.0-1.4.

**Symptom:** Figure legends say R_e; discussion text says R_0; readers conflate the two; comparisons to non-phylodynamic R_0 estimates inappropriate.

**Fix:** Use R_e (or R_t) consistently. Explicit footnote: "R_e estimated in the current epidemic context; the basic reproduction number R_0 was higher and is not what BDSKY estimates."

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Phylodynamic R_e and case-based R_t disagree | Sampling bias on one side; lineage-specific signal in phylodynamic; under-ascertainment in case data | Report both with explicit assumptions; investigate sampling profile |
| TreeTime R^2 = 0.5 but BEAST clock rate confidence wide | TreeTime point estimate vs Bayesian posterior with prior; informative prior may collapse | Compare BEAST clock to literature; check for prior dominance |
| BDSKY R_e implausibly high in early epochs | `origin` mis-specified; sampling proportion too high for early epoch | Re-check `origin` definition; allow sampling proportion to vary per epoch |
| Gubbins and ClonalFrameML give different recombination masks | Different model assumptions; both approximate | Either is defensible; run sensitivity by re-doing clock inference with the alternate mask |
| MASCOT and BEAST DTA disagree on migration | DTA is sampling-biased; MASCOT is sampling-aware | Trust MASCOT; report DTA only with caveat |
| UShER MAT subtree vs full BEAST disagree on TMRCA | Parsimony branch lengths shorter than ML | Re-estimate branch lengths on the UShER subtree with TreeTime before downstream |
| BICEPS Ne(t) trajectory differs from BSP on the same data | BICEPS analytic Ne integration vs BSP segment uniform | Trust BICEPS; BSP suffered from edge artifacts and slow mixing |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| TempEst root-to-tip R^2 minimum | >=0.3 | Rambaut 2016 *Virus Evol* 2:vew007 convention |
| BEAST ESS per parameter (single chain) | >=200 | Standard convention; necessary but not sufficient |
| BEAST burn-in | 10% of chain length | Convention; visually verify trace |
| MASCOT minimum sequences per deme | >=20 | Müller 2018 *Bioinformatics* 34:3843 identifiability requirement |
| BDSKY MCMC length (>=100 tips) | 10^7-10^8 states | Stadler 2013 *PNAS* 110:228; BDSKY mixes slowly |
| TreeTime `--clock-filter` default | 4 (SD multiplier on root-to-tip residual) | TreeTime convention; tighter for outlier-sensitive analyses |
| Gubbins input | `core.full.aln` (full positions, NOT `core.aln`) | Cannot estimate background SNP density from variable positions only |
| IQ-TREE ascertainment bias correction | `+ASC` for SNP-only alignment | IQ-TREE convention |
| S. pneumoniae core clock (recombination-masked) | ~1.5e-6 subs/site/year | Croucher et al *Science* 2013 literature; >5e-6 indicates unmasked recombination |
| SARS-CoV-2 clock rate (ancestral) | ~8e-4 subs/site/year | SARS-CoV-2 substitution-rate literature (varies by lineage) |
| M. tuberculosis clock rate | ~0.3-0.5 SNPs/genome/year | Walker 2013 *Lancet Infect Dis* 13:137 literature (varies by lineage) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Clock rate 3-5x literature | Unmasked recombination | Run Gubbins / ClonalFrameML first |
| BEAST chain stuck at single tree topology | Operator-weight imbalance or extreme prior | Check operator schedule; relax prior |
| `treetime --clock-filter False` rejected | `--clock-filter` takes a numeric SD multiplier | Pass a number (typically 4) or omit |
| `beast --threads 4` rejected | BEAST uses single-dash flags | `-threads 4` |
| BDSKY origin same as rootHeight | Tutorial confusion | Set `origin > rootHeight` per Stadler 2013 |
| MASCOT migration HPD spans 2+ orders | <20 seqs per deme | Pool demes; accept uncertainty |
| Single-chain ESS reported as proof of convergence | ESS necessary but not sufficient | Run multi-chain; combine post-overlap |
| `run_gubbins.py core.aln` (not `core.full.aln`) | Variable-only alignment | Use `core.full.aln`; cite Croucher 2015 |
| UShER branch lengths fed directly to BDSKY | Parsimony branch lengths biased low | Re-estimate via TreeTime or BEAST first |
| BEAST2 XML breaks on minor-release upgrade | XML fragile across versions | Pin BEAST + every BEAST2 package version |
| `--reroot best` produces implausible root | Sampling-biased outliers | Cross-check with TempEst; consider `--reroot oldest` or manual root |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Was temporal signal checked?" | TempEst R^2 reported; date-randomisation test run; both interpreted (date-randomisation can pass with narrow windows) |
| "Was recombination masked?" | Gubbins on `core.full.aln`; rebuilt tree on masked alignment; cite Croucher 2015 |
| "Why BDSKY and not BSP?" | BICEPS / BDSKY model sampling explicitly; BSP / Skygrid assume uniform sampling, biased under preferential surveillance |
| "How many BEAST chains?" | 3-4 independent chains; marginal posterior overlap checked; logcombiner only after agreement |
| "Was MASCOT or DTA used for migration?" | MASCOT (or MASCOT-GLM with covariates); DTA inherits Lemey 2009 sampling bias for source attribution |
| "BDSKY origin vs rootHeight?" | `origin > rootHeight` by Stadler 2013 convention; documented |
| "Is the clock rate consistent with the literature?" | Compared to species-specific reference; recombination-masked S. pneumoniae expected ~1.5e-6 |
| "Was the case-based R_t reconciled?" | Reported both; disagreement attributed to sampling bias (phylodynamic is lineage-specific; case data is population mean) |

## References

- Stadler T, Kühnert D, Bonhoeffer S, Drummond AJ (2013) Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). *Proc Natl Acad Sci USA* 110(1):228-233. doi:10.1073/pnas.1207965110
- Sagulenko P, Puller V, Neher RA (2018) TreeTime: maximum-likelihood phylodynamic analysis. *Virus Evol* 4(1):vex042. doi:10.1093/ve/vex042
- Rambaut A, Lam TT, Carvalho LM, Pybus OG (2016) Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). *Virus Evol* 2(1):vew007. doi:10.1093/ve/vew007
- Duchêne S, Duchêne D, Holmes EC, Ho SY (2015) The performance of the date-randomization test in phylogenetic analyses of time-structured virus data. *Mol Biol Evol* 32(7):1895-1906. doi:10.1093/molbev/msv056
- Ramsden C, Holmes EC, Charleston MA (2009) Hantavirus evolution in relation to its rodent and insectivore hosts: no evidence for codivergence. *Mol Biol Evol* 26(1):143-153. doi:10.1093/molbev/msn234
- Bouckaert R, Vaughan TG, Barido-Sottani J et al (2019) BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. *PLoS Comput Biol* 15(4):e1006650. doi:10.1371/journal.pcbi.1006650
- Bouckaert RR (2022) An Efficient Coalescent Epoch Model for Bayesian Phylogenetic Inference. *Syst Biol* 71(6):1549-1560. doi:10.1093/sysbio/syac015
- Müller NF, Rasmussen DA, Stadler T (2018) MASCOT: parameter and state inference under the marginal structured coalescent approximation. *Bioinformatics* 34(22):3843-3848. doi:10.1093/bioinformatics/bty406
- Lemey P, Rambaut A, Drummond AJ, Suchard MA (2009) Bayesian phylogeography finds its roots. *PLoS Comput Biol* 5(9):e1000520. doi:10.1371/journal.pcbi.1000520
- De Maio N, Wu CH, O'Reilly KM, Wilson D (2015) New routes to phylogeography: A Bayesian structured coalescent approximation. *PLoS Genet* 11(8):e1005421. doi:10.1371/journal.pgen.1005421
- Volz EM, Frost SDW (2014) Sampling through time and phylodynamic inference with coalescent and birth-death models. *J R Soc Interface* 11(101):20140945. doi:10.1098/rsif.2014.0945
- Croucher NJ, Page AJ, Connor TR et al (2015) Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Res* 43(3):e15. doi:10.1093/nar/gku1196
- Didelot X, Wilson DJ (2015) ClonalFrameML: Efficient inference of recombination in whole bacterial genomes. *PLoS Comput Biol* 11(2):e1004041. doi:10.1371/journal.pcbi.1004041
- Mostowy R, Croucher NJ, Andam CP et al (2017) Efficient inference of recent and ancestral recombination within bacterial populations. *Mol Biol Evol* 34(5):1167-1182. doi:10.1093/molbev/msx066
- Didelot X, Croucher NJ, Bentley SD, Harris SR, Wilson DJ (2018) Bayesian inference of ancestral dates on bacterial phylogenetic trees. *Nucleic Acids Res* 46(22):e134. doi:10.1093/nar/gky783
- Turakhia Y, Thornlow B, Hinrichs AS et al (2021) Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. *Nat Genet* 53(6):809-816. doi:10.1038/s41588-021-00862-7
- Minh BQ, Schmidt HA, Chernomor O et al (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol* 37(5):1530-1534. doi:10.1093/molbev/msaa015
- Legried B, Terhorst J (2022) A class of identifiable phylogenetic birth-death models. *Proc Natl Acad Sci USA* 119(35):e2119513119. doi:10.1073/pnas.2119513119
- Featherstone LA, Duchêne S (2023) Decoding the fundamental drivers of phylodynamic inference. *Mol Biol Evol* 40(6):msad132. doi:10.1093/molbev/msad132
- Walker TM, Ip CLC, Harrell RH et al (2013) Whole-genome sequencing to delineate Mycobacterium tuberculosis outbreaks: a retrospective observational study. *Lancet Infect Dis* 13(2):137-146. doi:10.1016/S1473-3099(12)70277-3

## Related Skills

- pathogen-typing - Defines isolates and clusters that feed phylodynamic inference
- transmission-inference - Phylodynamic R_e + transmission tree are complementary
- variant-surveillance - Lineage assignment runs upstream of skygrid / BDSKY by lineage
- phylogenetics/divergence-dating - Calibrated trees for non-pathogen contexts; deeper review of clock models
- phylogenetics/bayesian-inference - BEAST / RevBayes / MrBayes details beyond BDSKY
- phylogenetics/modern-tree-inference - IQ-TREE / RAxML topology before time-scaling
- phylogenetics/tree-io - Tree parsing and format conversion
- comparative-genomics/whole-genome-alignment - Core-genome alignment input for bacterial phylodynamics
- variant-calling/vcf-basics - Per-isolate variant calls feeding core-SNP alignment
- read-alignment/bwa-alignment - Read mapping upstream of variant calling
- data-visualization/multipanel-figures - Skyline / R_e trajectory plotting
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
