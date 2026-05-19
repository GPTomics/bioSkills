# Subgroup Analysis - Usage Guide

## Overview

Performs subgroup and heterogeneous treatment effect (HTE) analyses for clinical trials with explicit distinction between confirmatory ("assessment") subgroups eligible for label claims (EMA 2019) and discovery/exploratory subgroups. Covers Mantel-Haenszel pooling with Simpson's-paradox detection, interaction tests, RERI for additive interaction, Gail-Simon qualitative interaction test, modern data-adaptive HTE methods (STEPP, SIDES, causal forests, X/R-learners), Bayesian shrinkage (Dixon-Simon, EXNEX), and graphical multiplicity (Bretz-Maurer via gMCP).

## Prerequisites

```bash
pip install statsmodels scipy numpy pandas matplotlib scikit-learn
```

R is recommended for modern HTE methods (most are R-only):

```r
install.packages(c('grf', 'policytree', 'causalToolbox', 'personalized',
                   'SIDES', 'stepp', 'quint', 'partykit', 'gMCP', 'RBesT', 'brms'))
```

Python alternatives for HTE:

```bash
pip install econml causalml dowhy
```

## Quick Start

Tell your AI agent what you want to do:
- "Single model with treatment-by-subgroup interaction (NOT comparing per-subgroup p-values)"
- "Mantel-Haenszel pooled OR + Breslow-Day + forest plot to detect Simpson's paradox"
- "Causal forest in grf with honest splitting + RATE/AUTOC test (Yadlowsky 2025)"
- "STEPP sliding-window analysis for continuous biomarker subgroups"
- "Bayesian shrinkage via Dixon-Simon for adjusted subgroup estimates"
- "Graphical multiplicity in gMCP allocating alpha across primary + key secondary + pre-specified subgroup"

## Example Prompts

### Interaction tests (the correct way)

> "Fit a single logistic regression with treatment * age_group interaction. Test interaction p-value as the formal evidence for HTE. Do NOT compare per-subgroup p-values (statistically invalid)."

> "Compute RERI for additive interaction on the treatment-by-sex effect. Delta-method or bootstrap CI for RERI."

### Stratified analysis

> "Run MH pooled OR + Breslow-Day across study sites. Forest plot stratum-specific ORs to detect Simpson's paradox."

### Modern HTE

> "Causal forest in grf for HTE on response with patient features X. Use honesty=TRUE; report calibration test result and RATE/AUTOC omnibus p (Yadlowsky 2025). Bias-correct estimates via grf's variable_importance + cross-fitting."

> "STEPP sliding-window analysis for continuous biomarker (HbA1c). Use permutation supremum test for pattern flatness -- not naive pointwise CIs."

> "SIDES recursive partitioning with permutation-adjusted base-vs-complement p (Lipkovich 2011) for subgroup discovery."

### Bayesian shrinkage

> "Apply Dixon-Simon Bayesian shrinkage to my forest of subgroup effects. Centre tau prior at HalfNormal(0, 0.25). Report posterior shrunken estimates."

> "EXNEX basket trial with 5 strata; 0.5 EX / 0.5 NEX mixture per Neuenschwander 2016. Sensitivity over weights 0.1-0.9."

### Multiplicity

> "Pre-specified 6 subgroups in confirmatory trial. Allocate 20% of primary alpha to subgroup family (Dane et al 2019 EFSPI recommendation). Implement as gMCP graph with weighted alpha propagation."

### Forest plot

> "Generate forest plot of subgroup ORs (sex, age, race, severity, region) with overall pooled estimate as reference. Log scale; CI shown."

## What the Agent Will Do

1. Distinguish pre-specified ("assessment") vs discovery ("exploratory") subgroups per EMA 2019
2. Fit single model with interaction term (NOT separate per-subgroup models with p-value comparison)
3. For stratified designs, include strata in CMH or logistic regression
4. For continuous biomarkers, use STEPP (sliding-window) with permutation supremum test
5. For data-adaptive HTE, use causal forests with honest splitting + RATE/AUTOC + bias-correction
6. For basket trials or borrowing, EXNEX or robust MAP with sensitivity over mixture weights
7. Apply graphical multiplicity per Bretz-Maurer 2009 with pre-specified alpha allocation
8. Forest plot stratum/subgroup-specific ORs as the primary visual

## Tips

- **Senn's foundational point** (Senn 2018 *Nature* 563:619): observed between-patient response variation is NOT evidence of patient-level HTE. Don't conflate noise with heterogeneity.
- **Brookes 2004 4x penalty:** detecting a treatment-by-subgroup interaction requires ~4x the n for the main effect of similar magnitude. Most non-significant interaction tests are underpowered, not null.
- **Comparing per-subgroup p-values is statistically invalid.** Separate models have different power; p differences confound effect size with sample size. ALWAYS use a single model with interaction term.
- **Breslow-Day with k=3 strata has ~40% power.** Non-significance does NOT prove homogeneity. Supplement with forest plot AND LR interaction test from logistic regression.
- **CMH pooled OR can mask Simpson's paradox.** Always plot stratum-specific ORs.
- **EMA 2019 guideline** distinguishes "assessment subgroups" (pre-specified, regulatory weight) from "discovery subgroups" (hypothesis-generating only). Interaction tests are "neither necessary nor sufficient" for credibility.
- **Causal forests:** honest splitting (`honesty=TRUE`) is mandatory for valid CIs. Run `test_calibration()` AND `rate_omnibus()` for diagnostics. Rehill 2025 audit found ~70% of applied papers skip these.
- **STEPP with naive pointwise CIs is wrong** -- overlapping windows give correlated estimates. Use permutation supremum test.
- **Bayesian shrinkage (Hemmings-Koch 2019):** appropriate for replication PLANNING, not signal GENERATION. Shrinkage pre-emptively damps the heterogeneity being searched for.
- **EXNEX default 0.5/0.5 weights** allows substantial borrowing; sensitivity over weights essential.
- **Yadlowsky RATE/AUTOC 2025** single-p test: whether CATE ranking has predictive (not just prognostic) value. Should replace `test_calibration` as primary HTE omnibus check.
- **Sun et al 2012 BMJ 11 credibility criteria** is the canonical academic framework: pre-specified + significant interaction + biological plausibility + consistency across endpoints + replication.
- **Winner's curse (Sun 2010):** median observed effect in significant subgroups is ~2.4x trial-overall effect. Bayesian shrinkage or honest cross-validation corrects.
- **Stratified randomisation factors MUST appear in subgroup analyses** (Kahan-Morris 2012). Ignoring inflates Type-I.

## Related Skills

- clinical-biostatistics/categorical-tests - CMH, Breslow-Day
- clinical-biostatistics/effect-measures - Forest plots and effect estimation
- clinical-biostatistics/logistic-regression - Interaction terms in regression
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs in depth
- clinical-biostatistics/bayesian-trials - MAP/EXNEX/Berry hierarchical in depth
- clinical-biostatistics/trial-reporting - CONSORT 2025 + EMA 2019 reporting
- experimental-design/multiple-testing - General methods
- machine-learning/biomarker-discovery - HTE for biomarker subgroups
