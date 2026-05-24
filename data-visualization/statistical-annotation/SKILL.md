---
name: bio-data-visualization-statistical-annotation
description: Add p-value brackets, significance asterisks, and effect-size annotations to distribution plots using ggpubr, ggsignif, and statannotations with correct test selection (parametric vs non-parametric vs paired), multiple-testing adjustment, and rendering of negative results. Use when a boxplot/violin/raincloud needs in-figure statistical comparisons between groups.
tool_type: mixed
primary_tool: ggpubr
---

## Version Compatibility

Reference examples tested with: ggpubr 0.6+, ggsignif 0.6+, rstatix 0.7+, statannotations 0.6+ (Python), seaborn 0.13+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Statistical Annotation

**"Add p-values to my plot"** -> Render pairwise group comparisons as brackets with the correct statistical test (parametric vs non-parametric, paired vs unpaired, independent vs nested), adjusted for multiple testing, with rendering of significance as either numerical p OR asterisks. The choices that matter: which test is appropriate for the data, what multiple-testing adjustment applies, and whether to show n.s. (non-significant) results.

- R: `ggpubr::stat_compare_means`, `ggsignif::geom_signif`, `rstatix::t_test`/`wilcox_test`
- Python: `statannotations.Annotator`, `scipy.stats` directly

## The Single Most Important Modern Insight -- The Test Must Match the Data

Tool defaults are NOT data-appropriate. `ggpubr::stat_compare_means(method='t.test')` uses Welch's two-sample t-test assuming approximate normality and unequal variances. This is wrong when:

1. **Data are non-normal and N is small (<30):** use Mann-Whitney U (`method='wilcox.test'`).
2. **Data are paired:** use paired t-test or paired Wilcoxon (`paired = TRUE`).
3. **Comparing >2 groups:** ANOVA / Kruskal-Wallis with post-hoc, not all-pairs t-test (multiple-testing penalty).
4. **Data are nested (cells within patients, replicates within samples):** linear mixed model, NOT pairwise test.

The bracket-and-asterisk visual is the same; the underlying statistics are not. Choose the test deliberately.

## Decision Tree for Test Selection

| Question | Recommended test | Function |
|----------|------------------|----------|
| 2 unpaired groups, normal, N≥30 | Welch t-test | `t.test()`, `stat_compare_means(method='t.test')` |
| 2 unpaired groups, non-normal or small N | Mann-Whitney U (Wilcoxon rank-sum) | `wilcox.test()`, `stat_compare_means(method='wilcox.test')` |
| 2 paired groups | Paired t-test OR Wilcoxon signed-rank | `paired = TRUE` |
| 3+ groups, normal | One-way ANOVA + Tukey HSD post-hoc | `aov()`, `TukeyHSD()` |
| 3+ groups, non-normal | Kruskal-Wallis + Dunn post-hoc | `kruskal.test()`, `dunn.test()` |
| Nested data (cells in patients) | Linear mixed model | `lme4::lmer()` |
| Time-course / repeated measures | Repeated-measures ANOVA OR LMM | `nlme::lme()` |
| Two-way factorial | Two-way ANOVA + interaction term | `aov(y ~ a*b)` |
| Survival / time-to-event | Log-rank, NOT t-test | `survdiff()` |
| Categorical outcome | Chi-square OR Fisher exact | `chisq.test()`, `fisher.test()` |

## Multiple Testing Adjustment

For pairwise comparisons among K groups, there are K(K-1)/2 unadjusted p-values. Without adjustment, family-wise error rate inflates rapidly:
- 3 groups: 3 comparisons; α_FW = 14% at nominal 5%
- 4 groups: 6 comparisons; α_FW = 26%
- 6 groups: 15 comparisons; α_FW = 54%

```r
# rstatix supports per-comparison adjustment
library(rstatix)
df %>%
    pairwise_wilcox_test(value ~ group, p.adjust.method = 'bonferroni') %>%
    add_xy_position(x = 'group')
```

`p.adjust.method` options:
- `'bonferroni'` — strictest; controls FWER
- `'holm'` — stepwise Bonferroni; uniformly more powerful than Bonferroni
- `'BH'` (Benjamini-Hochberg) — FDR; less strict than FWER; standard for genomics
- `'fdr'` — alias for BH

For figure annotations, **holm** is the modern default — controls FWER and is more powerful than bonferroni. For a small number of pre-planned comparisons (≤3), Bonferroni is fine.

## ggpubr -- Standard ggplot2 Workflow

**Goal:** Add per-comparison p-value brackets between groups on a distribution plot, using a test appropriate to data shape and adjusting for multiple comparisons.

**Approach:** Build the base plot with `ggboxplot()`; add `stat_compare_means()` with explicit `method`, `comparisons`, `p.adjust.method`, and `label` arguments; render as asterisks (`p.signif`) for terse display or numeric (`p.format`) for precise display.

```r
library(ggpubr)

# Default boxplot + p-value bracket(s)
ggboxplot(df, x = 'group', y = 'value', color = 'group',
          add = 'jitter', palette = 'npg') +
    stat_compare_means(method = 'wilcox.test',           # explicit; default is t-test
                       comparisons = list(c('Control', 'Treatment'),
                                          c('Control', 'Vehicle'),
                                          c('Treatment', 'Vehicle')),
                       label = 'p.signif',               # 'p.signif' for asterisks; 'p.format' for numeric
                       p.adjust.method = 'holm',
                       method.args = list(alternative = 'two.sided'))
```

For an overall test plus pairwise:

```r
# Overall + per-comparison
ggboxplot(df, x = 'group', y = 'value', color = 'group') +
    stat_compare_means(method = 'kruskal.test',          # overall test
                       label.y = 1.05 * max(df$value)) +
    stat_compare_means(comparisons = pairs,
                       method = 'wilcox.test',
                       p.adjust.method = 'holm',
                       label = 'p.signif')
```

## ggsignif -- Lighter Alternative

```r
library(ggsignif)
ggplot(df, aes(group, value, fill = group)) +
    geom_boxplot() +
    geom_signif(comparisons = list(c('Control', 'Treatment')),
                test = 'wilcox.test',
                map_signif_level = TRUE,                  # asterisks vs numeric p
                step_increase = 0.1) +
    scale_fill_manual(values = c('#0072B2', '#D55E00'))
```

`map_signif_level = TRUE` converts p-values to asterisks per Wasserstein-Lazar 2016 convention:
- `***` p < 0.001
- `**` p < 0.01
- `*` p < 0.05
- `ns` p ≥ 0.05

For literal p-values, set `FALSE`.

## statannotations (Python)

```python
import seaborn as sns
from statannotations.Annotator import Annotator

ax = sns.boxplot(x='group', y='value', data=df, palette=['#0072B2', '#D55E00', '#009E73'])

pairs = [('Control', 'Treatment'),
         ('Control', 'Vehicle'),
         ('Treatment', 'Vehicle')]

annotator = Annotator(ax, pairs, data=df, x='group', y='value')
annotator.configure(test='Mann-Whitney',                    # 't-test_ind', 't-test_paired', 'Wilcoxon', etc
                    comparisons_correction='holm',
                    text_format='star',                     # 'star', 'simple', 'full'
                    line_height=0.02,
                    text_offset=0.5)
annotator.apply_and_annotate()
```

## Per-Method Failure Modes

### Default t-test on non-normal data

**Trigger:** `stat_compare_means(method='t.test')` (default) on log-distributed expression.

**Mechanism:** t-test assumes approximate normality; non-normal data with small N inflates Type-I.

**Symptom:** Significant p where rank test gives p > 0.05.

**Fix:** Switch to `method='wilcox.test'` for non-normal or small-N data. Verify normality with Shapiro-Wilk if borderline.

### Pairwise tests without adjustment

**Trigger:** Multiple bracket annotations with raw p-values.

**Mechanism:** K(K-1)/2 comparisons inflate FWER without adjustment.

**Symptom:** All-pairwise significant at nominal 0.05; doesn't replicate.

**Fix:** `p.adjust.method = 'holm'` (or 'bonferroni' or 'BH'). Document choice.

### Paired data tested as independent

**Trigger:** Before/after measurements in same subjects, tested with unpaired t-test.

**Mechanism:** Ignores within-subject correlation; loses power.

**Symptom:** Non-significant p where paired test gives significant.

**Fix:** `paired = TRUE` (R) or `t-test_paired` (statannotations). Verify subjects are correctly matched.

### Nested data tested with pairwise t

**Trigger:** Hundreds of cells per patient, tested as if each cell is independent.

**Mechanism:** Pseudoreplication — within-patient correlation ignored; p-values dramatically over-significant.

**Symptom:** p < 1e-50 from a dataset where the actual N is ~10 patients.

**Fix:** Linear mixed model (`lme4::lmer(value ~ group + (1|patient_id))`); aggregate to per-patient median first; or pseudobulk.

### Asterisks shown but p-values not reported anywhere

**Trigger:** `label = 'p.signif'` exclusively.

**Mechanism:** Reader cannot recover the actual p-value.

**Symptom:** Reviewer asks for exact p; not in figure or supplementary.

**Fix:** Either show numeric p (`label = 'p.format'`) or include test results table in supplementary.

### n.s. annotation hidden

**Trigger:** Showing only significant brackets, omitting non-significant.

**Mechanism:** Selective reporting biases interpretation.

**Symptom:** Reader assumes untested pairs were significant.

**Fix:** Either annotate all comparisons (with n.s. for non-significant) OR pre-specify which pairs are tested in the legend/caption.

### Reading effect size from p-value

**Trigger:** Conclusion "highly significant difference" from p = 1e-10 on a tiny effect.

**Mechanism:** Large N inflates significance for trivial differences.

**Symptom:** Effect size negligible despite extreme p.

**Fix:** Always report effect size (Cohen's d, Cliff's delta, median difference with CI) alongside p. The bracket should convey direction AND magnitude, not just significance.

## Reconciliation: When Tests Disagree

| Pattern | Cause | Action |
|---------|-------|--------|
| t-test significant; Wilcoxon n.s. | Outliers driving t-test; rank test robust | Trust Wilcoxon for non-normal data |
| Unpaired n.s.; paired significant | Within-subject correlation matters | Use paired if subjects are matched |
| Pairwise all-significant; ANOVA n.s. | Multiple-testing inflation in pairwise | ANOVA / Kruskal-Wallis is the overall test; pairwise post-hoc only after omnibus significant |
| Pseudo-replicated p < 1e-50; LMM p = 0.1 | Pseudoreplication | LMM is correct; pseudo-replicated p is meaningless |
| Bonferroni-adjusted n.s.; raw p < 0.05 | Adjustment correctly identified borderline | Trust adjusted; document the test family |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| α for FWER control | 0.05 family-wise | Standard |
| α for FDR control | 0.05 expected FDR (BH) | Benjamini-Hochberg 1995 |
| Asterisk convention | * <0.05, ** <0.01, *** <0.001 | Common practice |
| Bonferroni cutoff | 0.05 / K(K-1)/2 | Standard |
| Holm step-down | better than Bonferroni for all K | Holm 1979 |
| FDR (BH) | less strict than FWER | Genomics standard |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Reviewer asks "why t-test?" | Default not justified | Pre-justify test choice |
| Many pairwise-significant; doesn't replicate | No multiple-testing adjustment | Holm or BH |
| Effect "highly significant" but tiny | Large N inflates p | Report effect size |
| Asterisks only; no p-values | `label = 'p.signif'` exclusively | Show p.format OR provide table |
| Pseudoreplication inflated p | Cells treated as independent | LMM or pseudobulk |
| n.s. comparisons hidden | Selective reporting | Annotate all pre-specified pairs |
| Numeric p truncated to '<2.22e-16' | R default precision | Manual formatting or report as `< 2e-16` |

## References

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289-300.
- Dunn OJ. 1964. Multiple comparisons using rank sums. *Technometrics* 6(3):241-252.
- Holm S. 1979. A simple sequentially rejective multiple test procedure. *Scand J Stat* 6(2):65-70.
- Kassambara A. 2020. *Practical Statistics in R for Comparing Groups: Numerical Variables.* (ggpubr / rstatix tutorial).
- Wasserstein RL, Lazar NA. 2016. The ASA's statement on p-values: context, process, and purpose. *Am Stat* 70(2):129-133.

## Related Skills

- data-visualization/distribution-plots - Underlying box/violin/raincloud
- clinical-biostatistics/categorical-tests - Chi-square / Fisher tests for categorical outcomes
- clinical-biostatistics/effect-measures - Effect size to report alongside p
- experimental-design/multiple-testing - Methods for controlling FWER and FDR
