# Categorical Association Tests - Usage Guide

## Overview

Tests associations between categorical variables in clinical trial data with calibration to design (independent vs paired), sample size (asymptotic vs exact), and regulatory expectation. Covers Pearson chi-square (no continuity correction), Boschloo's exact (uniformly more powerful than Fisher per Mehta-Senchaudhuri 2003), CMH/Breslow-Day stratified analysis with Simpson's-paradox detection, modern McNemar variants (mid-p per Fagerland 2013), and calibrated CIs (Wilson, Newcombe, Miettinen-Nurminen).

## Prerequisites

```bash
pip install scipy statsmodels pingouin pandas numpy
```

R for production NI/regulatory work with Miettinen-Nurminen stratified CIs:

```r
install.packages(c('ratesci', 'exact2x2', 'Exact'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Test treatment-outcome association using Pearson chi-square (no Yates) with expected count check"
- "Boschloo's exact test for small 2x2 -- more powerful than Fisher's"
- "CMH stratified analysis with Breslow-Day homogeneity + stratum-specific forest plot to detect Simpson's paradox"
- "Mid-p McNemar's for paired binary endpoint (Fagerland 2013 recommended over exact conditional)"
- "Miettinen-Nurminen score CI for RD as the regulatory-standard alternative to Wald"

## Example Prompts

### Independent 2x2 testing

> "I have a 2x2 table with 100 patients per arm. Expected counts all > 5. Run Pearson chi-square with correction=False and report the effect size (phi/Cramer's V)."

> "My 2x2 has 12 patients per arm. Run Boschloo's exact test in scipy (more powerful than Fisher's; cite Mehta-Senchaudhuri 2003). Compare to Fisher's exact for transparency."

### Stratified analysis

> "I have 8 stratification sites in my RCT. Run CMH for the pooled OR with Breslow-Day homogeneity. Plot stratum-specific ORs in a forest plot to detect Simpson's paradox."

> "If Breslow-Day is non-significant, that doesn't mean homogeneity (low power). Add an LR interaction test from logistic regression as a second check."

### McNemar's paired binary

> "Pre/post treatment binary outcome on n=80 subjects. Discordant pairs = 18. Use asymptotic McNemar without continuity correction per Fagerland-Lydersen-Laake 2013."

> "Matched-pair NI trial for diagnostic device. Use Suissa-Shuster exact unconditional (R exact2x2) for 20-40% smaller sample size than exact conditional."

### Calibrated CIs

> "Compute Wilson score CI for a single proportion (60/100 responders) -- better coverage than Wald, especially for extreme p."

> "Newcombe-Wilson hybrid CI for the difference of proportions (45/60 vs 30/60) -- preferred over Wald."

> "Miettinen-Nurminen score CI for RD in my stratified NI margin assessment -- the regulatory standard."

### Effect sizes and multiplicity

> "Compute Cramer's V with bias correction via pingouin's chi2_independence."

> "Post-hoc pairwise comparisons across 4 outcome categories with Holm correction; explicit method='holm' (NOT default 'hs' = Holm-Sidak)."

## What the Agent Will Do

1. Inspect the contingency table; verify expected counts; decide chi-square vs Boschloo
2. For stratified designs, include stratification factors via CMH or logistic regression
3. For paired designs, choose asymptotic vs mid-p McNemar based on discordant pair count
4. Compute appropriate CI (Wilson for single proportion; MN/Newcombe for differences/ratios)
5. Compute effect sizes (phi, Cramer's V with bias correction)
6. Apply multiplicity correction (Holm/Bonferroni for confirmatory; BH-FDR for exploratory) with explicit method specification

## Tips

- **Fisher's exact is conservative because it double-conditions** on the margins (Mehta-Senchaudhuri 2003). Boschloo's exact uses Fisher's p-value as test statistic but conditions on only one margin -- uniformly more powerful at the same Type-I.
- **Yates' continuity correction is now discouraged** (D'Agostino 1988): inflates Type-II by ~10%. Correct fix for sparse 2x2 is Boschloo, not Yates'.
- **Cochran 1954 strict rule:** no expected cell < 1 AND no more than 20% of cells with expected < 5. Textbook "all >=5" is the conservative simplification.
- **Stratified randomisation factors MUST appear in analysis** (Kahan-Morris 2012 *Stat Med* 31:328). Ignoring inflates Type-I error up to 30%.
- **CMH pooled OR can mask Simpson's paradox** when stratum-specific ORs reverse direction. Always supplement with a forest plot of stratum ORs.
- **Breslow-Day with k=3 strata has ~40% power** for moderate heterogeneity. Non-significance does NOT prove homogeneity.
- **Fagerland-Lydersen-Laake 2013** showed exact conditional McNemar is over-conservative. Asymptotic without continuity correction (b+c >= 25) or mid-p are preferred.
- **Wald CIs have well-documented coverage failures** (Brown-Cai-DasGupta 2001). Wilson for single proportions; Newcombe-Wilson hybrid or Miettinen-Nurminen for differences/ratios.
- **`scipy.stats.fisher_exact` since 1.10 returns SAMPLE (unconditional) OR**, not conditional MLE. For R-matching conditional MLE, use `scipy.stats.contingency.odds_ratio(table, kind='conditional')`.
- **`statsmodels.stats.multitest.multipletests` default `method='hs'` is Holm-Sidak**, NOT Holm or Bonferroni. Always specify explicitly.
- **Table2x2 column ordering trap:** `pd.crosstab` orders columns alphabetically. If outcome is 0/1, the OR will be reciprocal. Reorder: `cross[[1, 0]]`.

## Related Skills

- clinical-biostatistics/effect-measures - Detailed OR/RR/RD with calibrated CIs
- clinical-biostatistics/logistic-regression - Covariate-adjusted alternative with marginal vs conditional
- clinical-biostatistics/subgroup-analysis - Stratified analysis with interaction terms
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs for confirmatory multiplicity
- clinical-biostatistics/trial-reporting - CONSORT 2025 and ICH E9(R1) reporting
- experimental-design/multiple-testing - General multiplicity correction
