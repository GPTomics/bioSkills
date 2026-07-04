---
name: bio-temporal-genomics-trajectory-modeling
description: Models continuous temporal trajectories from BULK or time-resolved omics where the x-axis is measured experimental time: penalized GAMs (mgcv) for smooth trends and changepoint detection (segmented, ruptures) for abrupt regime shifts. Use when deciding between a smooth GAM and a changepoint model; choosing the GAM distribution (nb() plus a library-size offset for raw counts vs Gaussian on vst/log-CPM); setting the basis-dimension ceiling k below the number of timepoints and letting REML pick wiggliness; handling residual autocorrelation across timepoints with corAR1/bam(rho=); testing whether two conditions' trajectories diverge with an ordered-factor difference smooth; and choosing a changepoint search/cost/penalty (Pelt/Binseg, l2/rbf). Not for single-cell pseudotime (see single-cell/trajectory-inference).
tool_type: mixed
primary_tool: mgcv
---

## Version Compatibility

Reference examples tested with: mgcv 1.9+, tradeSeq 1.16+, segmented 2.0+, ruptures 1.1+, numpy 1.26+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: the modeled quantity must be on a scale the family assumes. A Gaussian GAM is valid only on variance-stabilized/log-transformed expression (vst, rlog, log-CPM); raw RNA-seq counts require `family=nb()` with a `offset(log(library_size))`. Successive timepoints are correlated, so a plain `gam()` (which assumes independent residuals) inflates smooth-term significance unless the AR structure is modeled or the residual ACF is checked.

# Temporal Trajectory Modeling

**"Fit smooth curves to my gene expression over real time, compare trajectories, and find abrupt shifts"** -> model a continuous function of MEASURED time f(time), test whether it changes / differs between conditions, and locate discrete regime changes.
- R: `mgcv::gam()`/`gamm()`/`bam()` for penalized-spline GAMs; `segmented::segmented()` for slope breaks
- Python: `ruptures` for level/distribution changepoints

## The governing principle: measured time is not pseudotime, and the model class must match the mechanism

This skill models trajectories where the x-axis is ACTUAL experimental time (hours, days, developmental stage) or a pseudobulk value aggregated over real time. Time is measured and shared across every sample at a timepoint, so replicates are exchangeable, timepoints are few and fixed, and residual autocorrelation across ordered timepoints is real. This is categorically different from single-cell pseudotime, which is a latent per-cell ordering estimated with error and belongs to single-cell/trajectory-inference (and to tradeSeq's native use case). Conflating the two is the deepest error in the area.

Three decisions dominate correctness before any p-value is read:
1. Distribution. A Gaussian GAM on raw counts gives wrong SEs, wrong p-values, and can predict negatives. Model counts with `nb()` and a library-size offset, or fit Gaussian on a variance-stabilized scale.
2. Autocorrelation. Positive residual correlation across timepoints shrinks the effective sample size, so a plain `gam()` under-estimates SEs and over-calls temporal trends. Name the assumption and model it (`corAR1`/`bam(rho=)`) or at least inspect the residual ACF.
3. Mechanism. A smooth GAM assumes a gradually curving process; a changepoint model assumes a genuinely abrupt regime shift. Imposing changepoints on smooth data invents regime shifts; smoothing a true step Gibbs-rings over it. Match the model to the biology.

## Smooth GAM vs changepoint: choosing the model class

| Question | Model | Use when | Do NOT use when |
|----------|-------|----------|-----------------|
| Does expression change / curve over time? | GAM `s(time)` (mgcv) | process is gradual (induction/decay kinetics, developmental ramps) | the process is a discrete switch -> a smooth smears the discontinuity |
| Do two conditions' trajectories diverge? | ordered-factor difference smooth (mgcv) | testing whether treated shape departs from control | groups have no shared reference / only a constant offset differs (use a parametric term) |
| When does the regime shift (slope break)? | `segmented` (broken-line) | continuous piecewise-LINEAR change in slope | the shift is a level jump, not a slope change (use ruptures `l2`) |
| When does the regime shift (level/distribution)? | `ruptures` Pelt/Binseg | step change in mean (`l2`) or distribution (`rbf`) | the curve is smooth -> any liberal penalty fabricates breaks |
| Is a curve even warranted? | AIC/edf of `s(time)` vs linear | deciding non-linear vs linear is enough | over-interpreting edf near 1 as a real curve |

Methodology evolves; before committing verify current defaults and recommendations against the latest mgcv/ruptures/segmented documentation.

## mgcv GAM (R)

**Goal:** Fit a smooth non-linear curve to expression over measured time and test whether it changes, on the correct distributional scale.

**Approach:** Use a penalized regression spline `s(time)`; set the basis-dimension ceiling k generously but below the number of unique timepoints and let REML choose the realized wiggliness; use `nb()` + a library-size offset for raw counts, or Gaussian on a variance-stabilized scale.

```r
library(mgcv)

# k is a CEILING (max basis dimension), NOT the number of knots the curve will use.
# The REML-chosen penalty picks the realized wiggliness; edf (below) reports it.
# k must be < number of unique timepoints; realistically k <= (#timepoints - 1).
# method='REML': better-behaved objective than GCV, resists under/over-smoothing (Wood 2011).
fit <- gam(expression ~ s(time, k = 6, bs = 'tp'), data = gene_df, method = 'REML')

summary(fit)
# s.table columns: edf, Ref.df, F (Gaussian/unknown scale), p-value.
# edf ~ 1 => penalty shrank the smooth to linear; edf near k-1 => nearly full flexibility.
# The smooth-term p-value is APPROXIMATE (conditional on estimated lambda; Wood 2013):
# treat it as categorical significant/not, do not compare tiny magnitudes.
```

### Raw counts: NB family with a library-size offset

**Goal:** Model overdispersed RNA-seq counts on the count scale without violating the Gaussian assumption.

**Approach:** Fit `family=nb()` (mgcv estimates theta by REML) with `offset(log(library_size))` so the smooth describes rate, not depth.

```r
# Counts are mean-variance coupled and overdispersed: Gaussian is wrong on raw counts.
# nb() estimates theta jointly with the smoothing parameters under REML.
# offset(log(libsize)) absorbs sequencing depth so s(time) models expression rate.
fit_nb <- gam(counts ~ s(time, k = 6) + offset(log(library_size)),
              data = gene_df, family = nb(), method = 'REML')
# With a known-scale family the test column becomes Chi.sq, not F.
summary(fit_nb)$s.table
```

### Residual autocorrelation across timepoints

**Goal:** Prevent inflated smooth-term significance caused by correlation between successive timepoints.

**Approach:** Model a lag-1 AR structure grouped by the replication unit with `gamm(correlation=corAR1())`, or fix `rho` in `bam()` for genome-wide fits after reading the lag-1 residual ACF.

```r
# Plain gam() assumes independent residuals; positive AR(1) shrinks effective n,
# under-estimates SEs, and lets the smooth get too wiggly -> false temporal trends.
# corAR1 models e_t = rho * e_{t-1} + noise; group by subject/animal/plate.
fit_ar <- gamm(expression ~ s(time, k = 6),
               correlation = corAR1(form = ~ time | subject),
               data = gene_df, method = 'REML')
# fit_ar$gam holds the smooth; fit_ar$lme holds the correlation estimate.

# Genome-wide alternative: fit once without AR, read lag-1 residual ACF, set rho, refit.
# AR.start flags the first observation of each independent series.
fit_bam <- bam(expression ~ s(time, k = 6), data = gene_df,
               rho = 0.4, AR.start = series_start, method = 'fREML')
```

With independent biological replicates AT EACH timepoint the correlation is often weak or unidentifiable and plain `gam()` is defensible; a single series sampled repeatedly over many timepoints is where AR bites hardest. Always inspect the residual ACF before trusting the smooth p-value.

### Comparing conditions with an ordered-factor difference smooth

**Goal:** Directly test whether the treated trajectory's shape diverges from control, with its own p-value.

**Approach:** Make the grouping an ORDERED factor so `s(time, by=grp)` becomes a difference smooth (level minus reference); keep the reference global smooth AND the parametric main effect.

```r
# Unordered by= gives each group's curve vs zero -- NOT a divergence test.
# Ordered factor: s(time, by=grp) is the difference smooth; its single p-value
# directly tests whether the trajectories diverge. Extends cleanly to >2 groups.
gene_df$condition <- as.ordered(gene_df$condition)
fit_diff <- gam(expression ~ condition + s(time, k = 6) + s(time, k = 6, by = condition),
                data = gene_df, method = 'REML')
# The parametric 'condition' term is REQUIRED: centered smooths cannot carry the
# group's overall level, so without it a constant offset is misattributed to the smooth.
summary(fit_diff)
```

A numeric 0/1 `by=is_treated` indicator is a valid shortcut for a single 2-level contrast (the second smooth is the treatment deviation), but the ordered-factor form is the general, canonical idiom.

### Diagnostics: gam.check, k-index, concurvity

**Goal:** Decide whether the basis is adequate and whether smooth terms are mutually identifiable.

**Approach:** Read `gam.check()`/`k.check()`; respond to a low k-index by doubling k and refitting, not by reflexively cranking k; use `concurvity()` only for multi-smooth models.

```r
gam.check(fit)   # 4 residual plots + k.check(): reports k', edf, k-index, p-value

# k-index < 1 with a small p means residual pattern the basis is too rigid to capture.
# CORRECT diagnostic = double k and refit: if edf rises substantially, k was too low;
# if edf barely moves, k was fine and the low k-index reflects autocorrelation or
# a distributional problem -- do not just raise k.

# Concurvity = the smooth analog of collinearity; only meaningful for multi-term models.
# Near 1 => partial attribution between smooths is unstable; > 0.8 is a worry, not a cutoff.
concurvity(fit_diff, full = TRUE)
```

### Prediction and pointwise intervals

**Goal:** Visualize the fitted trajectory with an uncertainty band, within the sampled range only.

**Approach:** Predict on a fine grid with `se.fit=TRUE`; band = fit +/- 1.96*SE (pointwise, not simultaneous); never extrapolate.

```r
grid <- data.frame(time = seq(min(gene_df$time), max(gene_df$time), length.out = 200))
pred <- predict(fit, newdata = grid, se.fit = TRUE)
grid$fitted <- pred$fit
# 1.96*SE is a POINTWISE 95% band; whole-curve (simultaneous) coverage is < 95%,
# so overlapping condition bands are NOT a formal test -- use the difference smooth p-value.
grid$lower <- pred$fit - 1.96 * pred$se.fit
grid$upper <- pred$fit + 1.96 * pred$se.fit
# Beyond [min(time), max(time)] the spline and its SE diverge -- do not predict outside range.
```

### Genome-wide GAM fitting + FDR

**Goal:** Rank genes by temporal significance across the transcriptome.

**Approach:** Fit `s(time)` per gene, collect the smooth p-value, apply BH across genes (the per-gene p-values are approximate, so the FDR is approximate; permutation calibration is the gold standard for strong claims).

```r
# expr_mat here must be variance-stabilized (vst / log-CPM); for raw counts use the family=nb() + offset
# form above instead of this default-Gaussian fit, or the per-gene SEs and p-values are invalid.
results <- data.frame()
for (gene in rownames(expr_mat)) {
    df <- data.frame(expression = as.numeric(expr_mat[gene, ]), time = timepoints)
    fit <- gam(expression ~ s(time, k = 6), data = df, method = 'REML')
    s_tab <- summary(fit)$s.table
    results <- rbind(results, data.frame(gene = gene, edf = s_tab[, 'edf'],
                                         p_value = s_tab[, 'p-value']))
}
results$q_value <- p.adjust(results$p_value, method = 'BH')  # q<0.05: standard FDR floor
temporal_genes <- results[results$q_value < 0.05, ]
```

## tradeSeq (R/Bioconductor) -- off-label for bulk

tradeSeq is BUILT for single-cell pseudotime lineages, not bulk real-time. `fitGAM(counts, pseudotime, cellWeights, nknots)` expects a gene x cell count matrix, a cell x lineage pseudotime matrix, and cell x lineage soft-assignment weights, and fits an NB GAM per gene per lineage. Its tests (`associationTest`, `startVsEndTest`, `conditionTest`, `patternTest`) are keyed to pseudotime lineages.

```r
# Only reasonable when the design is genuinely a pseudobulk-over-a-lineage.
# For a standard bulk time-course with replicates at fixed timepoints, use mgcv directly:
# the same NB GAM, with full control over by= contrasts, offsets, and AR structure,
# and without single-cell scaffolding. nknots plays the same ceiling role as k (choose it with evaluateK()).
library(tradeSeq)
sce <- fitGAM(counts = count_mat, pseudotime = pt_mat, cellWeights = cw_mat, nknots = 6)
assoc_res <- associationTest(sce)   # matrix of Wald stat + df + p-value, NOT an SCE
```

## segmented (R) -- broken-line slope break

**Goal:** Locate a continuous change in SLOPE and test whether a break exists at all.

**Approach:** Pre-test with `davies.test` before fitting a break; estimate the breakpoint with `segmented()` from a starting value; limit to one break unless the data are dense.

```r
library(segmented)
lm_fit <- lm(expression ~ time, data = gene_df)

# davies.test H0: difference-in-slopes = 0 (no breakpoint). It searches candidate
# locations and corrects the minimum p for the search (slightly conservative).
# Decide WHETHER a break exists before interpreting one. pscore.test() is more powerful
# for a single break.
davies.test(lm_fit, seg.Z = ~time)

# segmented models a CONTINUOUS slope change (not a level jump); psi = starting value(s),
# NA auto-initializes. The estimator is a local search: sensitive to psi, fragile with
# multiple breaks (closely spaced breaks are non-identifiable).
seg_fit <- segmented(lm_fit, seg.Z = ~time, psi = NA)
summary(seg_fit)$psi   # breakpoint estimate + SE (break CI is approximate/often too narrow)
```

## ruptures (Python) -- level/distribution changepoints

**Goal:** Detect discrete times where the mean (or whole distribution) shifts.

**Approach:** Factorize as (search method) x (cost model) x (penalty); the penalty choice IS the number-of-changepoints choice; match the cost model to the shift type and estimate the noise variance, not the total variance.

```python
import numpy as np
import ruptures as rpt

signal = np.asarray(expression_values)

# model='l2' detects changes in MEAN (piecewise-constant level -- the natural choice for
#   step-like expression regimes). 'rbf' detects changes in the whole distribution
#   (mean AND variance) -- more general, hungrier for data, less interpretable.
# min_size=2: minimum segment length. Pelt is exact (O(n) via pruning); it takes a penalty
#   and RETURNS the number+location of breaks -- there is no separate 'how many' knob.
n = len(signal)
# Noise variance, NOT total variance: np.var(signal) includes between-regime variation, so
# a BIC-style penalty log(n)*np.var(signal) is too large and UNDER-detects. Estimate noise
# from lag-1 differences instead. BIC is derived for the l2/Gaussian-mean cost -- pairing it
# with model='rbf' is theoretically mismatched; use l2 with BIC, or calibrate rbf empirically.
sigma2 = np.var(np.diff(signal)) / 2.0
penalty = np.log(n) * sigma2
bkps = rpt.Pelt(model='l2', min_size=2).fit(signal).predict(pen=penalty)
# predict returns break indices INCLUDING the terminal index n:
n_changepoints = len(bkps) - 1   # bkps[:-1] are the actual break locations

# Binseg: greedy, approximate, fast; takes a KNOWN number of breaks (sanity check vs Pelt).
bkps_binseg = rpt.Binseg(model='l2', min_size=2).fit(signal).predict(n_bkps=2)
```

Guard against fabricated breaks: a liberal penalty always "finds" changepoints in a smooth ramp. Require a pre-test (a break exists) or compare a piecewise fit against a smooth-GAM fit by AIC -- if the smooth wins, the "changepoint" is a sampling-noise artifact. With few timepoints, be extremely skeptical of more than one break.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| p-values wrong / fitted curve predicts negative expression | Gaussian GAM on raw overdispersed counts | `family=nb()` + `offset(log(library_size))`, or fit Gaussian on vst/log-CPM |
| Many genes "significantly change over time" implausibly | residual autocorrelation inflates smooth-term significance | `gamm(..., correlation=corAR1(form=~time|subject))` or `bam(..., rho=, AR.start=)`; check residual ACF |
| Treating k as "the number of bends I want" | k is the flexibility CEILING, not realized complexity | set k generously (< #timepoints), let REML pick lambda; read `edf`, not k |
| Cranking k whenever k-index < 1 | low k-index can mean autocorrelation/heteroscedasticity, not low basis | double k and refit -- if edf jumps, raise k; if not, look at correlation/distribution |
| Reading p=1e-30 as thirty orders of certainty | smooth p-values are approximate (ignore full lambda uncertainty) | treat as categorical significant/not; apply BH FDR across genes |
| Unordered `by=` used "to test if curves differ" | it gives each group vs zero, not a divergence test | `as.ordered(condition)` -> the difference smooth's p-value IS the divergence test |
| `by=` smooth without the parametric main effect | centered smooths cannot carry the group level | include `condition +` alongside `s(time, by=condition)` |
| `tradeSeq::fitGAM` on a plain bulk time-course | tradeSeq is single-cell pseudotime machinery (needs cellWeights) | use `mgcv` directly for bulk real-time; reserve tradeSeq for pseudobulk lineages |
| ruptures under-detects real changepoints | `pen=log(n)*np.var(signal)` uses TOTAL variance -> penalty too large | estimate noise from `np.var(np.diff(signal))/2`; sweep the penalty |
| `model='rbf'` with a BIC (`log n * var`) penalty | BIC penalty is derived for the l2/Gaussian-mean cost | use `model='l2'` with BIC, or calibrate the rbf penalty empirically |
| Changepoints "found" in a clearly smooth ramp | a liberal penalty fabricates breaks in gradual data | pre-test with `davies.test`; compare piecewise vs smooth-GAM AIC; the smooth often wins |
| Fitted curve behaves wildly past the last timepoint | extrapolating a penalized spline beyond the data | predict only within `[min(time), max(time)]` |
| "The condition bands overlap, so no difference" | 1.96*SE bands are pointwise, not simultaneous | use the difference-smooth p-value or simultaneous (posterior-simulation) intervals |
| `segmented` with several `psi` gives unstable breaks | multiple breakpoints are weakly identifiable with few noisy points | limit to 1 break unless data are dense; supply good `psi` starts; check convergence |

## Related Skills

- temporal-clustering - group genes by trajectory shape after fitting
- circadian-rhythms - periodic (known-period) trajectory models rather than smooth trends
- periodicity-detection - discover unknown-period oscillation instead of a smooth trend
- differential-expression/timeseries-de - linear/spline model alternatives for temporal DE
- single-cell/trajectory-inference - single-cell pseudotime (latent inferred ordering), the case tradeSeq is built for

## References

- Wood SN. 2011. Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *J R Stat Soc B* 73(1):3-36. doi:10.1111/j.1467-9868.2010.00749.x. (REML smoothing-parameter selection, better-behaved than GCV.)
- Wood SN. 2013. On p-values for smooth components of an extended generalized additive model. *Biometrika* 100(1):221-228. doi:10.1093/biomet/ass048. (Smooth-term p-values are approximate; pointwise interval coverage.)
- Wood SN. 2017. *Generalized Additive Models: An Introduction with R*, 2nd ed. Chapman & Hall/CRC. ISBN 9781498728331. (Basis-penalty framework, gam.check, concurvity.)
- Pedersen EJ, Miller DL, Simpson GL, Ross N. 2019. Hierarchical generalized additive models in ecology: an introduction with mgcv. *PeerJ* 7:e6876. doi:10.7717/peerj.6876. (Global-plus-difference-smooth and factor-smooth condition-comparison structure.)
- Van den Berge K, Roux de Bezieux H, Street K, Saelens W, Cannoodt R, Saeys Y, Dudoit S, Clement L. 2020. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat Commun* 11(1):1201. doi:10.1038/s41467-020-14766-3. (tradeSeq: NB-GAM DE along pseudotime lineages, hence off-label for bulk.)
- Muggeo VMR. 2003. Estimating regression models with unknown break-points. *Stat Med* 22(19):3055-3071. doi:10.1002/sim.1545. (Broken-line estimator behind segmented and davies.test.)
- Killick R, Fearnhead P, Eckley IA. 2012. Optimal detection of changepoints with a linear computational cost. *J Am Stat Assoc* 107(500):1590-1598. doi:10.1080/01621459.2012.737745. (PELT exact O(n) penalized algorithm behind rpt.Pelt.)
- Truong C, Oudre L, Vayatis N. 2020. Selective review of offline change point detection methods. *Signal Processing* 167:107299. doi:10.1016/j.sigpro.2019.107299. (The cost x search x constraint taxonomy; the ruptures reference paper.)
