# Statistical Annotation - Usage Guide

## Overview

Statistical annotation overlays pairwise p-value brackets and significance asterisks onto distribution plots. The visual is the same regardless of test; the underlying statistics are not. Choosing the appropriate test (parametric vs non-parametric, paired vs unpaired, nested vs independent) and adjusting for multiple comparisons are mandatory. ggpubr (R) and statannotations (Python) are the modern interfaces.

## Prerequisites

```r
install.packages(c('ggpubr', 'ggsignif', 'rstatix'))
```

```bash
pip install statannotations seaborn scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Add Wilcoxon p-values between control and treatment to my boxplot"
- "Pairwise tests for 4 groups with Holm adjustment"
- "Switch from t-test to Mann-Whitney because the data are non-normal"
- "Show asterisks instead of numeric p"
- "Use paired test because before-after measurements in same subjects"

## Example Prompts

### Two-group comparison

> "Wilcoxon p-value bracket between Control and Treatment on the boxplot. Show numeric p, not asterisks."

### Multi-group with adjustment

> "Pairwise Wilcoxon tests across 4 conditions with Holm adjustment. Display as a bracket pyramid above the boxes with asterisks for significance."

### Paired test

> "Pre/post measurements paired by subject_id. Use paired Wilcoxon. Indicate paired-test in the figure legend."

### Nested data

> "I have ~200 cells per patient × 8 patients. Don't pseudo-replicate; aggregate to per-patient mean or use a linear mixed model."

### Overall + post-hoc

> "Kruskal-Wallis overall test in the top corner; per-pair Dunn post-hoc as brackets between groups."

## What the Agent Will Do

1. Identify N per group and overall data structure (paired? nested? repeated?).
2. Test normality (Shapiro-Wilk) if borderline N for normal-based choice.
3. Pick the appropriate test: t for normal+independent; Wilcoxon for non-normal or small N; paired versions when matched; ANOVA / Kruskal for 3+ groups.
4. Compute pairwise comparisons; adjust with Holm (default) or BH.
5. For nested data, switch to LMM or aggregate before pairwise testing.
6. Render p-values as brackets with consistent height; show numeric p OR asterisks per request.
7. Report effect size in caption/supplementary (Cohen's d, Cliff's delta, or median difference + CI).

## Tips

- **Default ggpubr test is t-test - not always right.** Specify `method = 'wilcox.test'` for non-normal or small N.

- **Use Holm (`p.adjust.method = 'holm'`) as the default adjustment.** More powerful than Bonferroni, controls FWER.

- **Paired tests require matched subjects.** `paired = TRUE` (R) or `t-test_paired` (statannotations). Verify pairing before invoking.

- **Nested data must NOT be pairwise-tested.** Pseudoreplication inflates p by orders of magnitude. Use LMM (`lme4::lmer`) or aggregate to per-cluster median.

- **Always report effect size alongside p.** Cohen's d for normal; Cliff's delta for non-parametric; median difference with bootstrap CI.

- **Don't hide n.s. comparisons.** Pre-specify the comparisons tested in the caption; show all results including n.s.

- **For >5 groups, omnibus first.** Kruskal-Wallis or ANOVA; pairwise post-hoc only if overall significant.

- **Asterisks vs numeric:** `label = 'p.signif'` for asterisks (terse); `label = 'p.format'` for numeric (precise). Provide both via figure + supplementary table.

- **`map_signif_level = TRUE` in ggsignif** converts p to asterisks following the * < 0.05, ** < 0.01, *** < 0.001 convention.

- **p < 1e-50 from large N often means small effect with massive N.** Effect size is the biological story; p is the statistical question.

## Related Skills

- data-visualization/distribution-plots - Underlying boxplot / violin / raincloud
- clinical-biostatistics/categorical-tests - Tests for categorical outcomes
- clinical-biostatistics/effect-measures - Effect size companion to p
- experimental-design/multiple-testing - FWER and FDR methods
