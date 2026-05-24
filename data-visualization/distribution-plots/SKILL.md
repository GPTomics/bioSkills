---
name: bio-data-visualization-distribution-plots
description: Plot per-group distributions of continuous data using boxplots, violins, beeswarms, quasirandom jitter, and raincloud plots with sample-size honesty (Weissgerber 2015), KDE-bandwidth awareness, and N-aware encoding choices. Use when comparing distributions across a small number of groups — expression per cluster, biomarker per arm, scores per condition — and the bar-of-mean default is misleading.
tool_type: mixed
primary_tool: ggplot2
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+, ggbeeswarm 0.7+, ggdist 3.3+, gghalves 0.1.4+, seaborn 0.13+, matplotlib 3.8+, ptitprince 0.3+ (Python raincloud).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Distribution Plots

**"Plot the distribution per group"** -> Render boxplot, violin, beeswarm, or raincloud calibrated to N per group, the underlying distribution shape, and the audience's ability to read each encoding. The default `geom_bar(stat='summary')` is the canonical misleading choice — Weissgerber 2015 *PLOS Biol* documented that 703 top physiology papers use bar-of-mean despite multiple distinct distributions producing identical bars.

- R: `ggplot2::geom_boxplot`, `ggplot2::geom_violin`, `ggbeeswarm::geom_quasirandom`, `ggdist::stat_halfeye`, `gghalves::geom_half_violin`
- Python: `seaborn.boxplot/violinplot/swarmplot/stripplot`, `ptitprince.RainCloud`

## The Single Most Important Modern Insight -- Bars of Means Lie

Weissgerber, Milic, Winham & Garovic 2015 *PLOS Biol* 13:e1002128 surveyed 703 papers in top physiology journals and found that bar-and-line graphs of means dominate, despite **many distinct distributions producing identical bar plots**. Bimodal data, skewed data, and data with outliers all collapse to the same bar height and error bar. The bar plot is a hypothesis test result rendered as visualization; the visualization should show the data.

The modern alternative is to **show every point** for n < 30, layer summary on top, and reserve summary-only plots for large N where points would overplot.

## Decision Tree by N per Group

| N per group | Recommended | Avoid |
|-------------|-------------|-------|
| 3-10 | Dot plot or jittered raw points + median bar | Bar of mean |
| 10-30 | Beeswarm OR quasirandom + box overlay | Bare boxplot (hides bimodality) |
| 30-200 | Raincloud (Allen 2019) OR box + jitter | Bare violin (default KDE bandwidth oversmooths) |
| 200-1000 | Letter-value plot (Hofmann 2017) OR violin with explicit bandwidth | Box alone (collapses tails) |
| >1000 | Density (KDE) or histogram + summary stats | Individual points (overplot) |

**Always annotate N** somewhere on the plot (caption, x-axis tick label, or stratum count).

## Box, Violin, Beeswarm, Raincloud -- The Four Standard Encodings

### Boxplot (Tukey 1977) -- summary only

```r
ggplot(df, aes(group, value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_manual(values = c('#0072B2', '#D55E00')) +
    labs(x = NULL, y = 'Expression') +
    theme_classic()
```

Box shows: median, IQR, 1.5×IQR whiskers, outliers. Hides: bimodality, sample size, density.

**Notched boxplot** (`notch = TRUE`): notches show 95% CI for median (±1.58·IQR/√n); non-overlapping notches roughly indicate distinct medians. Use with N ≥ 15.

### Violin -- density + summary

```r
ggplot(df, aes(group, value, fill = group)) +
    geom_violin(alpha = 0.7, trim = FALSE,
                bw = 'SJ') +                            # Sheather-Jones bandwidth
    geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
    scale_fill_manual(values = c('#0072B2', '#D55E00'))
```

**KDE bandwidth pitfall:** ggplot's default is Silverman's rule of thumb, which oversmooths bimodal data into a single mode. Use `bw = 'SJ'` (Sheather-Jones plug-in) for honest representation of multimodality.

**`trim = TRUE`** (default) cuts the violin at the data range — visually misleading because the violin's tails imply density extending beyond the data. `trim = FALSE` lets the KDE extend.

### Beeswarm / quasirandom -- every point shown deterministically

```r
library(ggbeeswarm)
ggplot(df, aes(group, value, color = group)) +
    geom_quasirandom(method = 'quasirandom', width = 0.3, alpha = 0.7) +
    scale_color_manual(values = c('#0072B2', '#D55E00')) +
    stat_summary(fun = median, geom = 'crossbar', width = 0.5, color = 'black')
```

Quasirandom (van der Corput sequence; Bostock implementation) gives reproducible jitter that fills space without random scatter. Beeswarm is similar but with collision avoidance. Both are deterministic — reruns produce identical layouts.

### Raincloud (Allen 2019) -- distribution + summary + raw

**Goal:** Show distribution (half-violin), summary (boxplot), and raw observations (jittered points) in a single per-group panel without occlusion.

**Approach:** Place a half-violin on one side, a thin boxplot in the middle, and jittered points on the other side via `gghalves::geom_half_violin` + `geom_boxplot` + `geom_half_point` with `position_nudge` offsets; flip to horizontal so the visual reads as a literal raincloud.

```r
library(gghalves)
ggplot(df, aes(group, value, fill = group, color = group)) +
    geom_half_violin(side = 'r', alpha = 0.7, position = position_nudge(x = 0.15)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7,
                 position = position_nudge(x = -0.05)) +
    geom_half_point(side = 'l', alpha = 0.5, size = 1.5, range_scale = 0.4,
                    position = position_nudge(x = -0.2)) +
    scale_fill_manual(values = c('#0072B2', '#D55E00')) +
    scale_color_manual(values = c('#0072B2', '#D55E00')) +
    coord_flip()                                          # horizontal "raincloud"
```

```python
import ptitprince as pt
import seaborn as sns
pt.RainCloud(x='group', y='value', data=df,
             palette=['#0072B2', '#D55E00'],
             bw='scott', cut=0,                          # bandwidth + trim
             width_viol=0.6, orient='h')
```

Raincloud = half-violin (distribution) + boxplot (summary) + jittered raw points. Allen 2019 *Wellcome Open Res* 4:63 — modern publication default for N 30-200.

### Letter-value plot (Hofmann-Wickham 2017)

```r
library(lvplot)
ggplot(df, aes(group, value, fill = group)) +
    geom_lv(k = 5, alpha = 0.7) +
    scale_fill_manual(values = c('#0072B2', '#D55E00'))
```

Extends Tukey's boxplot via additional letter-value quantiles (Hofmann, Wickham, Kafadar 2017 *J Comput Graph Stat* 26:469). For large N, the standard boxplot collapses tail structure; letter-value preserves it.

```python
sns.boxenplot(x='group', y='value', data=df,
              palette=['#0072B2', '#D55E00'])             # seaborn calls it boxenplot
```

### Stacked / split violin (paired comparisons)

```r
library(introdataviz)               # split-violin geom
ggplot(df, aes(group, value, fill = condition)) +
    geom_split_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, position = position_dodge(0.5), outlier.shape = NA)
```

For 2-condition comparison within each group, split-violin shows both densities back-to-back. More compact than dodged violins.

## Per-Method Failure Modes

### Bar of mean with SEM

**Trigger:** `geom_bar(stat = 'summary')` + `geom_errorbar(stat = 'summary', fun.data = mean_se)`.

**Mechanism:** Mean ± SEM collapses all distributional information; reader cannot assess bimodality, skew, or N.

**Symptom:** Reviewer asks to "show the data"; the figure must be redone.

**Fix:** Replace with raincloud, beeswarm, or boxplot+jitter. Show points for N < 30.

### Violin with default Silverman bandwidth oversmooths bimodality

**Trigger:** `geom_violin()` without specifying `bw`.

**Mechanism:** Silverman's rule of thumb assumes unimodal Gaussian; oversmooths bimodal data into a single peak.

**Symptom:** Single-cell expression bimodality (off / on) renders as a unimodal violin; biologically false.

**Fix:** `bw = 'SJ'` (Sheather-Jones plug-in) for honest bimodality. Note: `nrd0` IS Silverman; `nrd` (Scott) oversmooths less than Silverman but Sheather-Jones is preferred.

### Notched boxplot with too-small N

**Trigger:** `notch = TRUE` with N < 15 per group.

**Mechanism:** Notch can extend beyond Q1/Q3, producing visually-misleading "inside-out" notches.

**Symptom:** ggplot warning ("notch went outside hinges"); notches look weird.

**Fix:** Use notches only with N ≥ 15. For smaller N, show raw points instead.

### Boxplot hides outliers when jittered points are overlaid

**Trigger:** `geom_boxplot() + geom_jitter()` with default `outlier.shape = 19`.

**Mechanism:** Outliers render twice — once from boxplot (large dots), once from jitter (smaller dots) — visually duplicated.

**Symptom:** Some points appear bigger than others without reason.

**Fix:** `geom_boxplot(outlier.shape = NA)` when overlaying raw points.

### Trim = TRUE on violin misleads about tails

**Trigger:** `geom_violin()` default `trim = TRUE`.

**Mechanism:** Default trims violin at the data range; the visual still shows narrowing "tails" implying density extends slightly beyond the data.

**Symptom:** Reader infers density beyond observed range.

**Fix:** `trim = FALSE` to let KDE extend, OR explicitly cap with `coord_cartesian`. Document the choice.

### No N annotation

**Trigger:** Boxplot with no N per group reported.

**Mechanism:** Reader cannot assess statistical power; tiny N looks identical to large N at this encoding.

**Symptom:** Reviewer requests "show N per group."

**Fix:** Add N to x-axis tick label (`Control (n=12)`) or use `stat_summary(geom='text', fun.data = function(x) data.frame(label = paste('n=', length(x))))`.

### Wide raincloud at small N

**Trigger:** Raincloud applied with N = 5 per group.

**Mechanism:** KDE with N=5 is meaningless; violin shape is artifact of bandwidth.

**Symptom:** Smooth violin from 5 points; misleads about underlying distribution.

**Fix:** For N < 30, drop the violin half; use box + raw points only.

## Reconciliation: When Encodings Disagree

| Pattern | Cause | Action |
|---------|-------|--------|
| Bar of mean shows clear separation; raincloud shows overlapping distributions | Bars hide overlap | Use raincloud; bars exaggerate effect |
| Violin shows unimodal; histogram shows bimodal | Default Silverman bandwidth oversmooths | Re-render with `bw = 'SJ'` |
| Boxplot medians look distinct; t-test n.s. | Boxplot of small N is unreliable | Show raw points; rerun with appropriate non-parametric test |
| Notched boxplot notches non-overlap; rank test n.s. | Notch is an approximation, not a hypothesis test | Notches are heuristic only; use formal test |

**Operational rule:** for N < 30, show every point. For N 30-200, raincloud. For N > 200, letter-value or violin with explicit bandwidth. Always annotate N.

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| N for valid notched boxplot | ≥15 | Common practice |
| N to show raw points | <30 | Weissgerber 2015 |
| N where violin > box | >30 (with explicit bandwidth) | Visualization guidance |
| Whisker length (Tukey) | 1.5 × IQR | Tukey 1977 |
| Notch length (McGill 1978) | ±1.58 × IQR / sqrt(N) | McGill 1978 |
| KDE bandwidth (Sheather-Jones) | plug-in selector | Sheather-Jones 1991 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Bimodal data shown as unimodal violin | Default Silverman bandwidth | `bw = 'SJ'` |
| Duplicate large points on box + jitter | `outlier.shape` not suppressed | `geom_boxplot(outlier.shape = NA)` |
| Notches "inside-out" | N too small | Show raw points; remove notch |
| Raincloud looks smooth at N=5 | KDE meaningless at small N | Drop violin half; box + points only |
| No N visible | Default boxplot | Add `n=...` to x label or stat_summary text |
| Violin tails extend beyond data | `trim = TRUE` default + KDE bandwidth | `trim = FALSE` and cap with `coord_cartesian` |
| Bar of mean criticized in review | Weissgerber 2015 default failure | Replace with raincloud or box+jitter |

## References

- Allen M, Poggiali D, Whitaker K, Marshall TR, van Langen J, Kievit RA. 2019. Raincloud plots: a multi-platform tool for robust data visualization. *Wellcome Open Res* 4:63. doi:10.12688/wellcomeopenres.15191.1
- Hofmann H, Wickham H, Kafadar K. 2017. Letter-value plots: boxplots for large data. *J Comput Graph Stat* 26(3):469-477. doi:10.1080/10618600.2017.1305277
- McGill R, Tukey JW, Larsen WA. 1978. Variations of box plots. *Am Stat* 32(1):12-16.
- Sheather SJ, Jones MC. 1991. A reliable data-based bandwidth selection method for kernel density estimation. *J R Stat Soc B* 53(3):683-690.
- Streit M, Gehlenborg N. 2014. Points of view: Bar charts and box plots. *Nat Methods* 11(2):117.
- Tukey JW. 1977. *Exploratory Data Analysis.* Addison-Wesley.
- Weissgerber TL, Milic NM, Winham SJ, Garovic VD. 2015. Beyond bar and line graphs: time for a new data presentation paradigm. *PLOS Biol* 13(4):e1002128. doi:10.1371/journal.pbio.1002128

## Related Skills

- data-visualization/statistical-annotation - Add p-value brackets to distribution plots
- data-visualization/color-palettes - CVD-safe categorical palettes
- data-visualization/ggplot2-fundamentals - Grammar of graphics base
- single-cell/markers-annotation - Stacked / split violin for scRNA gene-by-cluster
- clinical-biostatistics/effect-measures - Effect size to accompany distribution
