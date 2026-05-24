# Forest and Funnel Plots - Usage Guide

## Overview

Forest plots summarize per-study effect estimates with CIs, ordered vertically with a pooled summary diamond. Funnel plots diagnose publication bias by plotting effect vs precision; asymmetry suggests missing small-N null-result studies (Egger 1997). The first question for any meta-analysis is heterogeneity (I², τ²) - pooling without addressing high I² is biologically uninterpretable. metafor is the R reference; ggforest (survminer) covers Cox subgroup forests; MendelianRandomization for MR-specific forests.

## Prerequisites

```r
install.packages(c('metafor', 'forestplot', 'survminer', 'MendelianRandomization'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Meta-analysis of 12 trials with random-effects REML; forest plot with log-scale OR axis"
- "Cox subgroup forest from this model with interaction p-values"
- "Funnel plot + Egger test for publication-bias diagnosis"
- "Contour-enhanced funnel showing significance regions"
- "MR forest with IVW, weighted-median, and MR-Egger methods"

## Example Prompts

### Standard meta-analysis

> "Meta-analyze 15 RCTs of intervention X on binary outcome Y. Use REML random-effects. Forest plot with study weights, log-OR axis, summary diamond, 95% prediction interval, I² and Q-test annotated."

### Cox subgroup forest

> "Subgroup forest of treatment HR across age, sex, stage, and biomarker subgroups from coxph fit. Test interaction explicitly; annotate p-values."

### Mendelian randomization forest

> "MR forest with IVW (primary), weighted-median, weighted-mode, and MR-Egger sensitivity. Annotate per-method p-values."

### Funnel and bias diagnostics

> "Contour-enhanced funnel with 90/95/99 significance regions. Egger test if k ≥ 10. Trim-and-fill as sensitivity, NOT primary."

## What the Agent Will Do

1. Compute per-study effect (log-OR / log-HR / β) and variance.
2. Fit random-effects REML model via metafor::rma (default unless homogeneous design justifies fixed-effect).
3. Compute I², τ², Q-test, and 95% prediction interval.
4. Render forest with log-scale x-axis for ratios; auto-sized boxes encoding weights; summary diamond at bottom.
5. Annotate caption with I², τ², Q-p, prediction interval.
6. Funnel plot + Egger test if k ≥ 10; trim-and-fill as sensitivity.
7. For Cox subgroup forests, add treatment × subgroup interaction term to model and annotate interaction p.
8. For MR, run multiple methods (IVW, weighted median, MR-Egger) and overlay in one forest.

## Tips

- **I² is the first question.** Pooling under I² > 75% without subgroup or meta-regression is biologically uninterpretable.

- **Random-effects REML is the modern default.** Fixed-effect assumes all studies estimate the same parameter - rarely true.

- **Log-scale x-axis for ratios (OR/HR/RR).** Linear scale visually misrepresents fold changes < 1 vs > 1.

- **Add 95% prediction interval** (`addpred = TRUE`) when I² > 30% - shows where a new study's effect is expected to fall.

- **Egger test requires k ≥ 10** (Sterne 2011). With fewer studies, visual funnel + contour-enhanced funnel is more reliable.

- **Trim-and-fill is sensitivity, not primary.** Present original alongside adjusted; the imputed studies are hypothetical.

- **Subgroup forest needs interaction test.** Visual differences across subgroups don't establish significant interaction; add the interaction term to the model.

- **Contour-enhanced funnel** (Peters 2008) distinguishes publication bias from other asymmetry - significance contours show whether missing studies are in p < 0.05 or p ≥ 0.05 regions.

- **MR forest requires triangulation.** Single-method MR is vulnerable to pleiotropy; report IVW + weighted-median + MR-Egger.

- **For >20 studies the forest becomes unreadable.** Group by subtype, or use a funnel only, or a caterpillar plot (sorted by effect).

## Related Skills

- clinical-biostatistics/effect-measures - HR / OR / RR definitions
- clinical-biostatistics/subgroup-analysis - Interaction tests
- causal-genomics/mendelian-randomization - MR-specific methods
- clinical-biostatistics/trial-reporting - Meta-analysis reporting standards
- data-visualization/color-palettes - Subgroup palette
