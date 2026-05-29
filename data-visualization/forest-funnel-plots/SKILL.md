---
name: bio-data-visualization-forest-funnel-plots
description: Build forest plots (HR, OR, RR, beta-coefficient summaries with CIs) and funnel plots (meta-analysis publication-bias diagnostics) using forestplot, metafor, ggforest, and MendelianRandomization with proper axis-scaling, summary-diamond placement, subgroup nesting, and Egger / trim-and-fill asymmetry tests. Use when summarizing effects across subgroups, trials, or instruments — meta-analysis, Mendelian randomization, subgroup HRs.
tool_type: r
primary_tool: metafor
---

## Version Compatibility

Reference examples tested with: metafor 4.4+, forestplot 3.1+, ggforestplot 0.1+ (subgroup forests), ggforest from survminer 0.4.9+, MendelianRandomization 0.10+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Forest and Funnel Plots

**"Summarize effects across studies / subgroups"** -> Render each effect estimate (HR, OR, RR, β) as a square (size = inverse variance / weight), horizontal bar (95% CI), and label, with an optional summary diamond at the bottom from a meta-analysis pool (fixed-effect or random-effects). The funnel plot diagnoses publication bias by plotting effect size vs precision; asymmetry indicates missing small-study-with-null-result publications (Egger 1997).

- R: `metafor::forest`, `metafor::funnel`, `forestplot::forestplot`, `survminer::ggforest` (Cox HR forests), `MendelianRandomization::mr_forest`

## The Single Most Important Modern Insight -- Heterogeneity Is the First Question

A pooled effect estimate is meaningless if the underlying studies are heterogeneous. The Higgins I² statistic (Higgins-Thompson 2002 *Stat Med* 21:1539) quantifies between-study variance; the conventional 25/50/75% interpretation tiers come from the Cochrane Handbook §10.10.2 (Higgins et al editors), NOT the original Higgins-Thompson paper which cautioned against rigid cutoffs. A meta-analysis with I² > 75% and a pooled effect must explain the heterogeneity (subgroup analysis, meta-regression) — pooling without explanation is statistically defensible but biologically unhelpful.

A forest plot's bottom must report: pooled estimate + 95% CI + I² + τ² (between-study variance) + Q-test p-value. Without these, the plot is a list of effects, not a meta-analysis.

## Decision Tree by Analysis Type

| Analysis | Tool | Pooling model | Forest method |
|----------|------|---------------|---------------|
| Single-trial subgroup HRs | survminer::ggforest | None (subgroup display) | Coxph object |
| Meta-analysis of binary outcomes | metafor::rma -> forest() | DerSimonian-Laird or REML random-effects | Standard forest |
| Meta-analysis of continuous outcomes | metafor::rma(yi, vi) | REML random-effects | Standard forest |
| Mendelian randomization | MendelianRandomization::mr_forest | Multiple MR methods | MR-specific forest |
| Subgroup forest with interaction p | metafor::rma + addpoly + interaction model | Subgroup REML | Nested forest |
| Network meta-analysis | netmeta::forest.netmeta | Bayesian or frequentist | Network forest |
| Cumulative meta-analysis (over time) | metafor::cumul + forest | – | Cumulative forest |

## Fixed-Effect vs Random-Effects Meta-Analysis

| Model | Assumption | When appropriate | Pooled estimate weight |
|-------|------------|-------------------|------------------------|
| Fixed-effect (Mantel-Haenszel, IVS) | All studies estimate the same true effect | Single mechanism, homogeneous design | 1 / within-study variance |
| Random-effects (DerSimonian-Laird, REML) | Studies estimate distinct true effects from a common distribution | Heterogeneous designs / populations | 1 / (within + between variance) |

**Use random-effects by default.** Fixed-effect assumes all studies estimate the *same* parameter, which is almost never true across multi-center trials with different populations. REML is the modern default (Viechtbauer 2005); DerSimonian-Laird is the older default still commonly seen.

## metafor::rma + forest -- The Reference Implementation

**Goal:** Pool study-level effect estimates with random-effects meta-analysis; render a forest plot with study weights, individual effect+CI, and pooled summary diamond.

**Approach:** Compute per-study yi (effect) and vi (sampling variance); fit REML random-effects model; pass to forest() with prediction interval if heterogeneity is non-trivial.

```r
library(metafor)

# Input: per-study effect (yi) and variance (vi)
# For OR: yi = log(OR), vi = SE(log(OR))^2
# For HR: yi = log(HR), vi = SE(log(HR))^2
res <- rma(yi = log_or, vi = log_or_se^2,
           data = studies, slab = paste(author, year),
           method = 'REML')

# I^2 and tau^2 in the summary
summary(res)
# I^2 (residual heterogeneity)
# tau^2 (estimated amount of (residual) heterogeneity)
# Q-test for heterogeneity

forest(res,
       atransf = exp,                           # display OR on natural scale
       at = log(c(0.25, 0.5, 1, 2, 4)),         # ticks at meaningful OR values
       refline = 0,                              # log(1) for OR/HR/RR
       xlab = 'Odds Ratio (95% CI)',
       header = c('Study', 'OR [95% CI]'),
       mlab = bquote(paste('RE Model (Q = ', .(round(res$QE, 2)),
                            ', df = ', .(res$k - 1),
                            ', p = ', .(format.pval(res$QEp, digits = 2)),
                            '; ', I^2, ' = ', .(round(res$I2, 1)), '%)')),
       addpred = TRUE)                          # prediction interval per Higgins 2009
```

`addpred = TRUE` adds a 95% prediction interval — where a new study's effect is expected to fall (Higgins-Thompson-Spiegelhalter 2009 *JRSS-A*). This is the most honest summary when I² > 30%.

## ggforest for Cox Subgroup Forests

```r
library(survminer)
fit <- coxph(Surv(time, status) ~ treatment + age + sex + stage, data = df)
ggforest(fit,
         data = df,
         main = 'Subgroup HRs',
         cpositions = c(0.02, 0.22, 0.4),
         fontsize = 0.7,
         refLabel = 'Reference',
         noDigits = 2)
```

ggforest produces a publication-ready subgroup forest from a coxph object. For pre-specified subgroup analyses (treatment × subgroup interaction), test interaction explicitly and annotate the p-value.

## Small-k Regime -- When Meta-Analysis Asymptotics Break

For k < 5 studies, the REML-based 95% CI from `metafor::rma()` is severely anti-conservative — it relies on chi-square asymptotics that fail with few studies. Use **Hartung-Knapp-Sidik-Jonkman (HKSJ)** adjustment (`test = 'knha'`):

```r
res_hksj <- rma(yi = log_or, vi = log_or_se^2, data = studies,
                method = 'REML', test = 'knha')           # HKSJ for k<5
```

HKSJ uses a t-distribution with k-1 degrees of freedom and adjusts SE via the Q-statistic — well-calibrated even at k=3. For k < 3 a meta-analysis is not advisable; report individual study effects in a forest plot without a pooled summary.

Also: I² is uninterpretable below k = 5 (Borenstein 2017 *Res Synth Methods* 8:5); the point estimate has wide CI dominated by k itself, not heterogeneity. **Do not report I² for k < 5.**

## Funnel Plot and Egger Test

**Goal:** Diagnose publication bias by visual asymmetry of effect size vs precision.

**Approach:** Plot effect (x) vs SE (inverted y); under no bias, points form a symmetric inverted funnel with the pooled estimate at the apex. Asymmetry suggests missing small-N null-result studies. Egger's regression test (Egger 1997 *BMJ* 315:629) formalizes the asymmetry.

```r
funnel(res,
       xlab = 'log(OR)',
       refline = res$b)

# Egger's test
regtest(res, model = 'lm', predictor = 'sei')
# significant p indicates asymmetry; suggests publication bias

# Trim-and-fill (Duval-Tweedie 2000) -- adjusts for asymmetry
res_tf <- trimfill(res)
forest(res_tf)
funnel(res_tf)
```

**Contour-enhanced funnel plot** (Peters 2008 *J Clin Epidemiol* 61:991) overlays significance contours (p < 0.10, < 0.05, < 0.01); asymmetry concentrated in "non-significant" regions indicates publication bias more specifically than generic asymmetry.

```r
funnel(res, level = c(90, 95, 99), shade = c('white', 'gray55', 'gray75'),
       refline = 0, legend = TRUE)
```

## Per-Method Failure Modes

### Pooling under high heterogeneity without explanation

**Trigger:** Random-effects meta-analysis pooled with I² > 75%; no subgroup or meta-regression.

**Mechanism:** Pooled estimate is a weighted average across substantively different effects; biologically meaningless.

**Symptom:** Pooled OR = 1.5 with 95% CI (1.2-1.8) but per-study effects range 0.3-5.0.

**Fix:** Run subgroup analysis or meta-regression to explain heterogeneity; report I², τ², and prediction interval; do NOT report a single pooled effect as the answer.

### Fixed-effect when studies are heterogeneous

**Trigger:** Default fixed-effect model on multi-population data.

**Mechanism:** Fixed-effect weights = 1/within-study variance, ignoring between-study variance.

**Symptom:** CI is misleadingly narrow; reviewer asks "why fixed effect with high I²?"

**Fix:** Switch to REML random-effects (`method = 'REML'`); document the choice.

### Egger test p-value over-interpreted

**Trigger:** k < 10 studies; significant Egger p taken as definitive publication bias.

**Mechanism:** Egger's test is underpowered with few studies; sensitive to single outliers.

**Symptom:** Conclusion "publication bias" from k=6 trials.

**Fix:** Egger requires k ≥ 10 (Sterne 2011 *BMJ* 343:d4002); for fewer studies, visual funnel + contour-enhanced funnel is more reliable.

### Trim-and-fill imputed studies presented as data

**Trigger:** Reporting trim-and-fill adjusted estimate as "the answer."

**Mechanism:** Trim-and-fill is a sensitivity analysis; imputed studies are hypothetical.

**Symptom:** Original pooled OR = 2.0; trim-and-fill adjusted to 1.5; report says "adjusted estimate is 1.5."

**Fix:** Present original AND trim-and-fill side-by-side; trim-and-fill is sensitivity, not primary.

### Subgroup forest without interaction test

**Trigger:** Subgroup HRs plotted; conclusion "treatment works in subgroup X."

**Mechanism:** Visual differences across subgroups don't establish significant interaction.

**Symptom:** Subgroup forest shows HR=0.5 in subgroup A, HR=1.0 in subgroup B; no formal test.

**Fix:** Add treatment × subgroup interaction term to the model; report interaction p; cite Brookes 2001 / Wang 2007 for subgroup analysis caveats.

### Forest plot axis on linear scale for ratios

**Trigger:** OR/HR/RR plotted with linear x-axis.

**Mechanism:** Ratios are multiplicatively symmetric; linear axis compresses < 1 effects.

**Symptom:** OR = 0.5 (halving) appears smaller than OR = 2 (doubling) on a linear scale, even though they are biologically equivalent.

**Fix:** Always log-scale the x-axis for ratios. metafor's `atransf = exp` + `at = log(c(0.25, 0.5, 1, 2, 4))` is the canonical pattern.

### Weights not visible (point sizes uniform)

**Trigger:** Default forestplot package without weight encoding.

**Mechanism:** Reader cannot tell study influence on pool.

**Symptom:** A 5-patient pilot looks visually equivalent to a 5000-patient trial.

**Fix:** Use metafor `forest()` which auto-encodes weight via box size. forestplot package needs `boxsize =` argument.

## Reconciliation: When Methods Disagree

| Pattern | Cause | Action |
|---------|-------|--------|
| Fixed-effect significant; random-effects n.s. | High heterogeneity inflates RE variance | Trust random-effects when I² > 30% |
| Egger n.s. but funnel looks asymmetric | k < 10 -> Egger underpowered | Trust visual; report contour-enhanced funnel |
| Trim-and-fill imputes many studies | Severe asymmetry | Caution; sensitivity, not primary |
| Subgroup forest suggests effect modification; interaction test n.s. | Visual difference does not establish formal interaction | Trust interaction test |
| MR forest shows divergent estimates across methods | Pleiotropy or weak instruments | Run MR-Egger, weighted median, mode-based (sensitivity); cite Bowden 2015 |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| I² substantial heterogeneity | > 50% | Cochrane Handbook §10.10.2 (Higgins et al editors) |
| I² considerable heterogeneity | > 75% | Cochrane Handbook §10.10.2 (Higgins et al editors) |
| Egger test min k | ≥ 10 | Sterne 2011 *BMJ* |
| Prediction interval (where new study lands) | report when I² > 30% | Higgins-Thompson-Spiegelhalter 2009 |
| Forest x-axis | log scale for ratios | Convention |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Pooled estimate with I² = 90% | No heterogeneity exploration | Subgroup / meta-regression |
| Linear x-axis for OR forest | Ratios should be log-symmetric | `atransf = exp, at = log(...)` |
| Uniform point sizes | Weights not encoded | metafor::forest auto-encodes; forestplot needs boxsize |
| Egger from k=5 | Underpowered | k ≥ 10 for Egger |
| Trim-and-fill as primary | Sensitivity, not primary | Present both; document |
| Subgroup effect without interaction test | Visual ≠ test | Add interaction term |
| MR forest with single method | Pleiotropy risk | Triangulate methods |

## References

- Bowden J, Davey Smith G, Burgess S. 2015. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *Int J Epidemiol* 44(2):512-525.
- Brookes ST, Whitley E, Peters TJ, et al. 2001. Subgroup analyses in randomised controlled trials: quantifying the risks of false-positives and false-negatives. *Health Technol Assess* 5(33):1-56.
- DerSimonian R, Laird N. 1986. Meta-analysis in clinical trials. *Control Clin Trials* 7(3):177-188.
- Duval S, Tweedie R. 2000. Trim and fill: a simple funnel-plot–based method of testing and adjusting for publication bias in meta-analysis. *Biometrics* 56:455-463.
- Egger M, Davey Smith G, Schneider M, Minder C. 1997. Bias in meta-analysis detected by a simple, graphical test. *BMJ* 315(7109):629-634.
- Higgins JPT, Thompson SG. 2002. Quantifying heterogeneity in a meta-analysis. *Stat Med* 21(11):1539-1558.
- Higgins JPT, Thompson SG, Spiegelhalter DJ. 2009. A re-evaluation of random-effects meta-analysis. *JRSS-A* 172(1):137-159.
- Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L. 2008. Contour-enhanced meta-analysis funnel plots help distinguish publication bias from other causes of asymmetry. *J Clin Epidemiol* 61(10):991-996.
- Sterne JAC, Sutton AJ, Ioannidis JPA, et al. 2011. Recommendations for examining and interpreting funnel plot asymmetry in meta-analyses of randomised controlled trials. *BMJ* 343:d4002.
- Viechtbauer W. 2010. Conducting meta-analyses in R with the metafor package. *J Stat Softw* 36(3):1-48.
- Hartung J, Knapp G. 2001. A refined method for the meta-analysis of controlled clinical trials with binary outcome. *Stat Med* 20(24):3875-3889.
- IntHout J, Ioannidis JPA, Borm GF. 2014. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. *BMC Med Res Methodol* 14:25.
- Borenstein M, Higgins JPT, Hedges LV, Rothstein HR. 2017. Basics of meta-analysis: I² is not an absolute measure of heterogeneity. *Res Synth Methods* 8(1):5-18.
- Higgins JPT, Thomas J, Chandler J, et al (editors). *Cochrane Handbook for Systematic Reviews of Interventions* (current version). Section 10.10.2 — interpretation tiers for I².

## Related Skills

- clinical-biostatistics/effect-measures - HR / OR / RR / NNT definitions
- clinical-biostatistics/subgroup-analysis - Interaction tests for subgroup HRs
- causal-genomics/mendelian-randomization - MR-specific forest + sensitivity
- clinical-biostatistics/trial-reporting - CONSORT and meta-analysis reporting
- data-visualization/color-palettes - Palette for multi-study or subgroup forests
