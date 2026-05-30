---
name: bio-experimental-design-power-analysis
description: Calculates statistical power for high-dimensional genomics experiments (bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, proteomics) under negative-binomial count models using RNASeqPower, PROPER, and simulation via powsimR, distinguishing per-gene from marginal (transcriptome-wide) power, the role of mean expression and dispersion, and the sequencing-depth-versus-replicate tradeoff. Covers simulation as the honest default for overdispersed counts, FDR-aware average power versus single-test power, observed/post-hoc power as an anti-pattern, and the winner's-curse / Type-S / Type-M consequences of underpowering. Use when planning replicate number for a sequencing experiment, deciding whether to add depth or samples, choosing closed-form versus simulation power, estimating power from pilot dispersions, or justifying replication in a grant. For clinical-trial power see clinical-biostatistics/power-and-sample-size; for the inverse sample-size question see experimental-design/sample-size.
tool_type: r
primary_tool: RNASeqPower
---

## Version Compatibility

Reference examples tested with: RNASeqPower 1.42+, PROPER 1.34+, powsimR 1.2+ (GitHub), DESeq2 1.42+, edgeR 4.0+, pwr 1.3+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt to the actual API. Notes: `RNASeqPower::rnapower()` solves for whichever of `n` or `power` is omitted; PROPER is a multi-step pipeline (`RNAseq.SimOptions.2grp` -> `simRNAseq` -> `runSims` -> `comparePower`); powsimR is GitHub-only and its `estimateParam`/`Setup`/`simulateDE` signatures drift — pin a commit SHA for reproducible work. Verify each against the installed help before relying on argument names.

# Power Analysis for Genomics Experiments

**"How many replicates does my sequencing experiment need?"** -> Compute the probability of detecting a biologically meaningful effect given replicate number, sequencing depth, and biological variability — modeling counts as negative-binomial and recognizing that power is a per-gene quantity, not one number for the whole transcriptome.
- R: `RNASeqPower::rnapower()` — closed-form NB power/sample size; `PROPER`, `powsimR` — simulation from the mean-dispersion trend

## The Single Most Important Modern Insight -- Genomics Power Is Per-Gene; Simulate, and Never Report Observed Power

Power in a sequencing experiment is not a single number. It is a per-gene quantity that depends on that gene's mean expression and dispersion, so the honest summary is the **marginal (average) power** across the expression distribution at a target FDR — the expected discovery rate. A single coefficient of variation plugged into a closed-form formula mis-states power for low- and high-expressed genes alike, because dispersion varies systematically with the mean; the defensible default for count data is **simulation from the empirical mean-dispersion trend** (PROPER, Wu 2015 *Bioinformatics* 31:233; powsimR, Vieth 2017 *Bioinformatics* 33:3486). The second rule is negative: **observed (post-hoc) power is information-free.** Computed from the effect a study actually estimated, it is a one-to-one function of the p-value and cannot explain a null result (Hoenig & Heisey 2001 *Am Stat* 55:19). Power is a design-stage quantity, computed for hypothesized effects before data exist. Underpowering does not merely miss true effects — it makes the significant ones overstate magnitude (Type-M) and sometimes reverse sign (Type-S), lowering the chance a significant call is real (Button 2013 *Nat Rev Neurosci* 14:365; Gelman & Carlin 2014 *Perspect Psychol Sci* 9:641).

## Algorithmic Taxonomy

| Approach | Model | Tool | Strength | Fails / costs when |
|----------|-------|------|----------|--------------------|
| NB closed-form | negative-binomial, single CV/dispersion | `RNASeqPower::rnapower` | fast; transparent; grant-ready | one CV cannot represent the mean-dispersion trend |
| Simulation, parametric | NB with mean-dispersion relationship | `PROPER` | honest marginal power + EDR at target FDR | needs a dispersion model / pilot |
| Simulation, empirical | resampled from pilot (incl. dropout) | `powsimR` | bulk AND scRNA-seq; realistic | GitHub-only; heavier; version drift |
| Gaussian closed-form | t-test / Cohen's d | `pwr::pwr.t.test` | per-feature ATAC/proteomics after transform | wrong for raw counts; ignores overdispersion |
| Effect-inflation design analysis | retrodesign for Type-S/Type-M | `retrodesign` (Gelman) | exposes exaggeration in noisy small-n | needs a plausible true effect |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Bulk RNA-seq, pilot data available | PROPER/powsimR simulation from pilot dispersions | matches the real mean-dispersion trend |
| Bulk RNA-seq, no pilot, quick grant number | `rnapower()` with a literature CV, stated as approximate | transparent; flag as conservative-to-rough |
| scRNA-seq cross-condition DE | powsimR on a pseudobulk model; power scales with samples | population power is set by donors, not cells |
| ATAC/ChIP/methylation per-region | NB simulation (PROPER-style) or pwr after variance-stabilizing | overdispersed counts; per-region power |
| Proteomics (continuous, log-abundance) | `pwr::pwr.t.test` per protein with missingness caveat | Gaussian after transform; MNAR matters |
| Justifying a null result post-hoc | report CI / effect size, NOT observed power | post-hoc power is uninformative (Hoenig-Heisey) |
| Fixed budget: depth vs replicates | favor replicates past ~10-20M mapped reads | biological variance dominates (Liu 2014) |
| Clinical-trial endpoint | -> clinical-biostatistics/power-and-sample-size | regulated regime, different machinery |

## Closed-Form NB Power -- RNASeqPower

**Goal:** Get a fast, transparent power or replicate number for bulk RNA-seq from depth, biological CV, and fold change.

**Approach:** Supply per-gene depth, biological coefficient of variation, the fold change to detect, and alpha; supply `n` to get power, or `power` to get the required `n`. Treat the result as a single-gene approximation and sanity-check against simulation.

```r
library(RNASeqPower)
# depth = reads/gene; cv = biological coefficient of variation; effect = fold change
rnapower(depth = 20, n = 5, cv = 0.4, effect = 2, alpha = 0.05)          # solves for POWER
rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.80)   # solves for n per group
```

## Simulation-Based Power -- the Honest Default for Counts

**Goal:** Estimate marginal power and the true realized FDR across the whole expression distribution, accounting for the mean-dispersion trend.

**Approach:** Build (or fit from pilot) a simulation model of counts with a realistic dispersion-mean relationship and DE-effect distribution, simulate many datasets at each candidate sample size, run the intended DE test, and read the average power at the target FDR.

```r
library(PROPER)
sim_opts <- RNAseq.SimOptions.2grp(ngenes = 20000, p.DE = 0.05,
                                   lOD = 'cheung', lBaselineExpr = 'cheung')  # empirical dispersion/expr priors
sims <- runSims(Nreps = c(3, 5, 8, 12), sim.opts = sim_opts, nsims = 50,
                DEmethod = 'edgeR')
powr <- comparePower(sims, alpha.type = 'fdr', alpha.nominal = 0.05,
                     stratify.by = 'expr', delta = log(1.5))          # delta is NATURAL-log lfc in PROPER; marginal power by expression stratum
summaryPower(powr)
```

## Depth vs Replicates -- the Budget Question

For bulk RNA-seq differential expression, sequencing depth shows diminishing returns once it is adequate — Liu, Zhou & White 2014 (*Bioinformatics* 30:301) found the inflection near **~10 million mapped reads** in MCF7 (commonly generalized to a 10-20M band) — whereas adding biological replicates improves power across the whole range. Under a fixed budget, allocate to more biological units before more depth. ATAC/ChIP have their own depth floors (library complexity, peak detection), but the principle holds: biological variance, not read count, limits discovery once depth is adequate.

## CV / Dispersion Guidelines (estimate from pilot when possible)

| Material | Typical biological CV | Source / note |
|----------|----------------------|---------------|
| Cell lines (technical replicates) | 0.1-0.2 | low biological variability |
| Inbred mice | 0.2-0.3 | moderate |
| Primary cells / donor-derived | 0.3-0.4 | donor-dependent |
| Human population samples | 0.3-0.5 | high; Hart 2013 *J Comput Biol* 20:970 default examples |

These are starting points, not substitutes for a pilot estimate; real dispersion is study-specific and a literature CV can be off by a factor of two (estimate via DESeq2/edgeR `estimateDispersions` — see experimental-design/sample-size).

## Per-Method Failure Modes

### Single CV for the whole transcriptome
- **Trigger:** one `cv` plugged into `rnapower()` for all genes.
- **Mechanism:** dispersion varies with mean expression; a single CV mis-states low/high-expressed genes.
- **Symptom:** simulation gives materially different power than the closed form.
- **Fix:** simulation-based power (PROPER/powsimR) from the mean-dispersion trend.

### Observed (post-hoc) power
- **Trigger:** "non-significant, but observed power was 0.3, so add samples."
- **Mechanism:** observed power is a monotone function of the p-value (Hoenig-Heisey 2001).
- **Symptom:** circular reasoning that adds nothing to the CI.
- **Fix:** report effect size + CI; do prospective power for the next study.

### Powering to the expected (or pilot-observed) effect
- **Trigger:** setting the effect to the hoped-for or pilot point estimate.
- **Mechanism:** the pilot estimate is itself noisy; building it in bakes in the winner's curse.
- **Symptom:** chronic underpowering; inflated significant effects (Type-M).
- **Fix:** power to the minimum biologically meaningful effect; propagate pilot variance, not its mean.

### Depth instead of replicates
- **Trigger:** "we will sequence deeper rather than add samples."
- **Mechanism:** past ~10-20M reads, biological variance dominates technical (Liu 2014).
- **Symptom:** deep libraries, still underpowered.
- **Fix:** add biological replicates.

### scRNA-seq power computed on cells
- **Trigger:** "100k cells from 2 patients gives huge power."
- **Mechanism:** population DE power is set by the number of biological samples; cells are pseudoreplicates.
- **Symptom:** power estimate wildly optimistic; results do not replicate.
- **Fix:** power on a pseudobulk model over donors (powsimR); see randomization-blocking.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Power >= 0.80 standard; >= 0.90 for pivotal | convention | tolerable Type-II risk |
| Depth saturates ~10-20M mapped reads for DE | Liu 2014 *Bioinformatics* 30:301 | biological variance then dominates |
| >=6 biological replicates recover most true DE | Schurch 2016 *RNA* 22:839 | n=3 misses many true DE at realistic effects |
| Observed power is a function of the p-value | Hoenig-Heisey 2001 *Am Stat* 55:19 | never use it to interpret a null |
| Type-M exaggeration large in noisy small-n | Gelman-Carlin 2014 *Perspect Psychol Sci* 9:641 | significant effects overstated |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Closed-form and simulation power disagree | single CV vs mean-dispersion trend | use simulation for the reported number |
| "Underpowered (observed power 0.3)" to excuse a null | post-hoc power fallacy | report CI; prospective power only |
| Deep libraries still underpowered | depth over replicates | add biological replicates |
| scRNA-seq power absurdly high | power computed on cells | pseudobulk power over donors |
| Significant effect far larger than literature | winner's curse from underpowering | design analysis (Type-S/Type-M); replicate |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Where did the CV come from?" | estimated from pilot dispersions (DESeq2); literature value used only as a conservative cross-check |
| "Why simulation rather than a formula?" | count power is per-gene; simulation captures the mean-dispersion trend and reports marginal power at the target FDR |
| "Is the study powered?" | marginal power >= 0.8 at FDR 0.05 for the minimum meaningful fold change; power curve provided |
| "Why not just sequence deeper?" | depth saturates ~10-20M reads (Liu 2014); replicates added instead |
| "Observed power of the null?" | observed power is uninformative (Hoenig-Heisey); CI on the effect reported instead |

## References

- Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher JP. 2013. Calculating sample size estimates for RNA sequencing data. *J Comput Biol* 20:970-978.
- Wu H, Wang C, Wu Z. 2015. PROPER: comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics* 31:233-241.
- Vieth B, Ziegenhain C, Parekh S, Enard W, Hellmann I. 2017. powsimR: power analysis for bulk and single cell RNA-seq experiments. *Bioinformatics* 33:3486-3488.
- Liu Y, Zhou J, White KP. 2014. RNA-seq differential expression studies: more sequence or more replication? *Bioinformatics* 30:301-304.
- Schurch NJ, Schofield P, Gierliński M, et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22:839-851.
- Hoenig JM, Heisey DM. 2001. The abuse of power: the pervasive fallacy of power calculations for data analysis. *Am Stat* 55:19-24.
- Button KS, Ioannidis JPA, Mokrysz C, Nosek BA, Flint J, Robinson ESJ, Munafò MR. 2013. Power failure: why small sample size undermines the reliability of neuroscience. *Nat Rev Neurosci* 14:365-376.
- Gelman A, Carlin J. 2014. Beyond power calculations: assessing Type S (sign) and Type M (magnitude) errors. *Perspect Psychol Sci* 9:641-651.
- Ioannidis JPA. 2005. Why most published research findings are false. *PLoS Med* 2:e124.

## Related Skills

- sample-size - The inverse problem: minimum replicates for a target power at a target FDR
- randomization-blocking - The experimental unit defines what is replicated; blocking changes error variance
- batch-design - Account for batch/blocking factors in the power model
- differential-expression/deseq2-basics - Estimating dispersions from pilot data for the power model
- single-cell/preprocessing - Pseudobulk model underlying scRNA-seq power
- clinical-biostatistics/power-and-sample-size - Power for regulated clinical-trial endpoints
