---
name: bio-experimental-design-randomization-blocking
description: Structures biological experiments so inference is valid by construction, covering Fisher's principles (randomization, replication, local control), the experimental-vs-observational unit distinction and pseudoreplication (Hurlbert 1984; Lazic 2018), randomization mechanics (complete, restricted, stratified, rerandomization, run-order), blocking layouts (randomized complete block, Latin square, incomplete block), factorial designs and interactions, and the split-plot/nested error strata hidden inside multi-batch genomics. Use when deciding the experimental unit and what counts as a replicate, planning randomization and run order, choosing a blocked/factorial/split-plot/nested layout, avoiding pseudoreplication in cell-culture or animal studies, or specifying the random-effects structure of the analysis model. For assigning samples to sequencing batches/lanes/plates and batch-effect correction see experimental-design/batch-design; for regulated clinical-trial randomization see clinical-biostatistics.
tool_type: r
primary_tool: designit
---

## Version Compatibility

Reference examples tested with: designit 0.5+, lme4 1.1-35+, lmerTest 3.1+, pwr 1.3+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt the example to the actual API rather than retrying. designit is an R6 package whose `BatchContainer$new()`, `optimize_design()`, and `*_score_generator()` signatures evolve between releases; confirm against the installed vignette (`vignette(package = 'designit')`) before relying on argument names.

# Randomization and Blocking

**"Design the experiment so the statistics will be valid"** -> Decide what the experimental unit is, randomize treatments to those units, replicate the unit (not the measurement), and remove known nuisance variation by blocking — so that the analysis model mirrors how the experiment was actually run.
- R: `designit::optimize_design()` for constrained randomization; `lme4::lmer()` / `lmerTest` for the matching mixed model
- The design and the analysis are one decision: "analyze as randomized"

## The Single Most Important Modern Insight -- The Experimental Unit, Not the Measurement, Is the n

The most consequential and most violated idea in biological design: the **experimental unit (EU)** is the smallest entity *independently assigned to a treatment* — and it, not the number of measurements, is the sample size for inference. Lazic 2018 *PLoS Biol* 16:e2005282 separates three entities: the **biological unit** (what conclusions are about), the **experimental unit** (what randomization acts on, = the n), and the **observational unit** (what is measured). When observational units are counted as independent replicates, the standard error shrinks illegitimately and p-values become meaningless — **pseudoreplication** (Hurlbert 1984 *Ecol Monogr* 54:187). Ten thousand cells from three mice are n = 3, not n = 10,000, for a between-mouse question; mice co-housed in a cage dosed through the chow make the *cage* the EU, not the mouse. Lazic et al. found ~46% of surveyed animal studies pseudoreplicated. The fix is structural: model the design's hierarchy (random effects) or aggregate to the EU before testing — pseudoreplication is, formally, an omitted random effect.

A second, deeper point (Fisher): **randomization is what licenses the p-value.** It supplies the physical basis for the error term and converts systematic lurking-variable bias into random error balanced *in expectation*. Model-based tests are approximations to the randomization distribution. Skip randomization and the causal claim rests entirely on assumptions.

## Algorithmic Taxonomy

| Design | Controls / estimates | When to use | Fails / costs when |
|--------|----------------------|-------------|--------------------|
| Completely randomized (CRD) | error variance only | units homogeneous; no known nuisance | inefficient if real nuisance structure exists |
| Randomized complete block (RCBD) | one known nuisance (day, litter, chip, donor) | nuisance factor identifiable and blockable | costs error df; harmful if block variance is ~0 ("blocking on noise") |
| Latin square | two orthogonal nuisances (day x technician) | n² runs affordable for n treatments | assumes no interaction among row/col/treatment |
| (Balanced) incomplete block | one nuisance, block smaller than #treatments | plate/chip holds fewer samples than treatments | analysis more complex; needs balance for efficiency |
| Factorial | main effects + interactions, "hidden replication" | >1 factor; interaction is of interest | #runs grows multiplicatively |
| Fractional factorial / screening | main effects under sparsity-of-effects | many factors, few runs (Plackett-Burman) | aliases effects; cannot resolve all interactions |
| Split-plot | two EU sizes, two error strata | one factor hard to randomize finely (lane, incubator, batch) | wrong error term if analyzed as a flat factorial -> anti-conservative |
| Nested / hierarchical | variance components across levels | sub-sampling within units (cells in mice in cages) | pseudoreplication if the nesting is ignored |
| Repeated measures | within-unit change over time | longitudinal sampling of the same EU | a split-plot in time; needs the within-unit error term |

## Decision Tree by Scenario

| Scenario | Recommended structure | Why |
|----------|----------------------|-----|
| Treatment given per animal, one tissue measured each | CRD or RCBD; n = animals | EU = animal |
| Many cells measured per animal, between-animal question | nested; aggregate to per-animal (pseudobulk) before testing | EU = animal, cells are observational units |
| Treatment delivered per cage (chow/water), several mice/cage | EU = cage; block or model cage as random | randomization acted on the cage |
| Two factors of interest (genotype x drug) | factorial; estimate the interaction | main effects uninterpretable if interaction is large |
| One factor fixed per run (incubator temp, sequencing lane) | split-plot; whole-plot = run, sub-plot = sample | two error strata; test whole-plot against whole-plot error |
| Known batch/day nuisance, all conditions fit per block | RCBD; include block in the model | removes nuisance from error; "analyze as randomized" |
| Plate holds fewer samples than conditions | incomplete block + include block term | balance preserves estimability |
| Assigning samples to sequencing batches/lanes | -> experimental-design/batch-design | constrained sample-to-batch allocation lives there |
| Regulated clinical trial randomization | -> clinical-biostatistics | confirmatory/regulated regime out of scope |

## Choosing and Counting the Experimental Unit

**Goal:** Identify the EU and therefore the true n before any power or analysis decision.

**Approach:** Trace the randomization: the EU is the smallest entity to which a treatment level was independently assigned. Anything measured below that level is an observational unit and is summarized (mean/sum) up to the EU, or modeled as a nested random effect — never counted as an independent replicate.

```r
# Between-condition question with multiple cells per donor:
# the donor is the experimental unit, NOT the cell.
# Correct: aggregate observational units to the EU, then test on EU-level values.
library(dplyr)
eu_level <- cells |>
  group_by(donor, condition) |>
  summarise(value = mean(measurement), .groups = 'drop')   # one row per experimental unit
# n for inference = number of donors per condition, not number of cells
```

## Randomization Mechanics

**Goal:** Assign treatments to units with a documented random mechanism, optionally restricted to guarantee balance on known factors.

**Approach:** Use a seeded pseudo-random generator (never "haphazard" order, which aliases treatment with processing position/time). For known prognostic factors, restrict the randomization (block/stratify) and then *include those factors in the model*. When finite-sample imbalance matters, rerandomize against a pre-specified balance criterion (Morgan & Rubin 2012 *Ann Stat* 40:1263) or use minimization for sequential enrollment (Pocock & Simon 1975 *Biometrics* 31:103).

```r
set.seed(20260528)                      # record the seed for reproducibility
units <- data.frame(id = sprintf('S%02d', 1:24),
                    block = rep(c('day1','day2','day3'), each = 8))

# Restricted (block) randomization: randomize treatment WITHIN each block
units$treatment <- ave(units$id, units$block,
                       FUN = function(ids) sample(rep(c('ctrl','treat'),
                                                      length.out = length(ids))))
# Also randomize RUN ORDER so processing position is not confounded with treatment
units$run_order <- sample(nrow(units))
```

## Blocking and Local Control

**Goal:** Remove a known nuisance source from the error term to sharpen the treatment comparison.

**Approach:** Group units into homogeneous blocks (day, litter, chip, donor), randomize treatments within block, and add the block as a term in the model. The paired t-test is the special case of an RCBD with block size 2. Block only on factors with real between-block variation; blocking on a noise factor spends error df for nothing.

```r
library(designit)                        # constrained assignment; verify API vs installed vignette
bc <- BatchContainer$new(dimensions = list(block = 3, position = 8))
bc <- assign_in_order(bc, samples = units)
bc <- optimize_design(
  bc,
  scoring = osat_score_generator(batch_vars = 'block',
                                 feature_vars = c('treatment')))  # balance treatment across blocks
```

## Split-Plot and Nested Designs -- the Genomics Trap

A **split-plot** has two experimental-unit sizes and therefore two error strata: a *whole-plot* factor that is hard to randomize finely (incubator temperature, the sequencing run/lane, the 10x chip, the staining batch) and a *sub-plot* factor randomized within each whole plot (the individual sample, the genotype). Analyzing a split-plot as a flat factorial uses the wrong, too-small error term for the whole-plot factor and gives **anti-conservative** tests for exactly the factor that was hardest to replicate. In genomics the lane/run/chip is almost always a whole plot; "batch effects" are frequently a split-plot structure to be modeled, not a nuisance to scrub.

**Goal:** Match the model's random-effects structure to the design's randomization structure.

**Approach:** Encode each randomization level as a random effect; fixed effects carry the questions. Crossed vs nested structure determines the denominator for each fixed effect; with few EUs use Satterthwaite or Kenward-Roger degrees of freedom (lmerTest / pbkrtest).

```r
library(lme4); library(lmerTest)
# Whole plot = run (random); sub-plot factor = condition (fixed); cells nested in sample
fit <- lmer(expression ~ condition + (1 | run/sample), data = df)   # run, and sample within run
anova(fit)                               # Satterthwaite df via lmerTest
```

## Factorial Designs and Interactions

A factorial design crosses factors so every observation informs every main effect ("hidden replication") and, uniquely, estimates **interactions** — the joint action one-factor-at-a-time (OFAT) cannot see. When an interaction is large, main effects are not interpretable alone; reporting a main effect while ignoring a strong interaction is the most common misreading of a 2x2 design. OFAT is less efficient and silently assumes additivity.

## Per-Method Failure Modes

### Pseudoreplication (observational units counted as n)
- **Trigger:** treating cells/wells/sections/technical aliquots as independent replicates.
- **Mechanism:** units within an EU are correlated; the SE is computed as if they were independent (Hurlbert 1984; Lazic 2018).
- **Symptom:** implausibly small p-values that fail to replicate; reviewers ask "what is n?".
- **Fix:** aggregate to the EU (pseudobulk) or add the EU as a random effect; the n is the number of EUs.

### Split-plot analyzed as a flat factorial
- **Trigger:** lane/run/incubator factor crossed with a within-run factor, fit with one error term.
- **Mechanism:** whole-plot factor tested against sub-plot error (too small).
- **Symptom:** the hard-to-randomize factor looks significant on thin evidence.
- **Fix:** two-stratum model; whole-plot factor uses whole-plot error (`(1 | run)`).

### Haphazard assignment mistaken for randomization
- **Trigger:** processing units "in the order they arrived".
- **Mechanism:** order aliases treatment with time/position/temperature gradients.
- **Symptom:** apparent treatment effect tracks run order.
- **Fix:** seeded PRNG assignment; randomize run order too; record the seed.

### Blocking on a noise factor
- **Trigger:** adding a block term with negligible between-block variance.
- **Mechanism:** spends error df without removing variance.
- **Symptom:** power lower than the unblocked design.
- **Fix:** block only on factors with documented between-block variation.

### Over-/under-specified random effects
- **Trigger:** maximal random structure that will not converge, or a structure missing a randomization level.
- **Mechanism:** maximal protects Type-I but may be singular (Barr 2013); too-lean inflates Type-I (pseudoreplication).
- **Symptom:** singular-fit warnings, or anti-conservative tests.
- **Fix:** keep the structure justified by the design; for small EU counts prune by a selection criterion (Matuschek 2017) and report the choice.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| EU = level of independent treatment assignment | Hurlbert 1984; Lazic 2018 | defines the n for inference |
| Paired design = RCBD with block size 2 | Fisher; standard | pairing is blocking |
| Latin square needs n² runs for n treatments | standard design theory | controls two nuisances orthogonally |
| Use Kenward-Roger/Satterthwaite df when EU count is small (roughly < ~10/group) | Kenward & Roger 1997 *Biometrics* 53:983 | naive F df are anti-conservative with few units |
| Maximal random effects for confirmatory; prune for small samples | Barr 2013; Matuschek 2017 | Type-I protection vs convergence/power tradeoff |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Significant result that will not replicate | pseudoreplication (cells as n) | aggregate to EU or add EU random effect |
| Whole-plot factor over-significant | split-plot analyzed flat | two-stratum mixed model |
| Treatment effect tracks processing order | no run-order randomization | seeded randomization of run order |
| Blocked design analyzed without block term | "design but don't analyze" | include block in the model; analyze as randomized |
| Main effect reported despite strong interaction | factorial misread | interpret simple effects within the interaction |
| Mixed model singular fit | over-specified random effects | prune to the design-justified structure (Matuschek 2017) |

## References

- Hurlbert SH. 1984. Pseudoreplication and the design of ecological field experiments. *Ecol Monogr* 54:187-211.
- Lazic SE, Clarke-Williams CJ, Munafò MR. 2018. What exactly is 'N' in cell culture and animal experiments? *PLoS Biol* 16:e2005282.
- Blainey P, Krzywinski M, Altman N. 2014. Points of significance: replication. *Nat Methods* 11:879-880.
- Krzywinski M, Altman N. 2014. Points of significance: designing comparative experiments. *Nat Methods* 11:597-598.
- Krzywinski M, Altman N. 2014. Points of significance: analysis of variance and blocking. *Nat Methods* 11:699-700.
- Morgan KL, Rubin DB. 2012. Rerandomization to improve covariate balance in experiments. *Ann Stat* 40:1263-1282.
- Pocock SJ, Simon R. 1975. Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial. *Biometrics* 31:103-115.
- Barr DJ, Levy R, Scheepers C, Tily HJ. 2013. Random effects structure for confirmatory hypothesis testing: keep it maximal. *J Mem Lang* 68:255-278.
- Matuschek H, Kliegl R, Vasishth S, Baayen H, Bates D. 2017. Balancing Type I error and power in linear mixed models. *J Mem Lang* 94:305-315.
- Kenward MG, Roger JH. 1997. Small sample inference for fixed effects from restricted maximum likelihood. *Biometrics* 53:983-997.
- Auer PL, Doerge RW. 2010. Statistical design and analysis of RNA sequencing data. *Genetics* 185:405-416.

## Related Skills

- batch-design - Assigning samples to sequencing batches/lanes and batch-effect correction
- sample-size - The experimental unit defines what is replicated and counted
- power-analysis - Blocking and nesting change the effective error variance
- multiple-testing - The design fixes what counts as a family of tests
- single-cell/preprocessing - Pseudobulk aggregation to the donor (experimental unit) for scRNA-seq
- differential-expression/deseq2-basics - The DE model that consumes the design's structure
- clinical-biostatistics/power-and-sample-size - Randomization and design in the regulated-trial regime
