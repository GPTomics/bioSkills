---
name: bio-data-visualization-upset-plots
description: Build UpSet plots to visualize set intersections beyond 4 sets (where Venn fails) using ComplexUpset (modern, ggplot2-grammar) or the unmaintained UpSetR, with explicit cardinality vs degree sorting, attribute panels, and query highlighting. Use when comparing overlap across many gene sets, peak sets, variant lists, or any set membership matrix where Venn diagrams become illegible.
tool_type: mixed
primary_tool: ComplexUpset
---

## Version Compatibility

Reference examples tested with: ComplexUpset 1.3+ (R, Krassowski), UpSetR 1.4.0 (last 2019 release; effectively unmaintained), upsetplot 0.9+ (Python).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# UpSet Plots

**"Show set intersections for 4+ sets"** -> Replace Venn diagrams (which become illegible past 4 sets) with UpSet (Lex 2014 *IEEE TVCG* 20:1983). The display: a matrix of dots indicating which sets participate in each intersection, with a vertical bar above each column showing intersection size and horizontal bars on the left showing per-set total size. Sort by intersection size (cardinality) for "biggest overlaps first" or by degree (number of sets) for grouped layout.

- R: `ComplexUpset::upset` (Krassowski; ggplot2-native, **recommended**), `UpSetR::upset` (Conway 2017; legacy, unmaintained)
- Python: `upsetplot.UpSet`

## The Single Most Important Modern Insight -- UpSetR Is Effectively Unmaintained

UpSetR (Conway, Lex, Gehlenborg 2017 *Bioinformatics* 33:2938) is the original R implementation but has had no CRAN release since v1.4.0 (2019). ComplexUpset (Krassowski; CRAN active through 2025-07) is the actively maintained ggplot2-grammar replacement. For new work in 2026, **prefer ComplexUpset**. Caveat: ggplot2 4.0 (mid-2025) broke ComplexUpset's `upset()` function (issue #213); pin to compatible versions until patched.

The Lex 2014 paper and underlying UpSet visualization concept are not affected — the visualization is the same; the difference is which R package implements it best in the current ecosystem.

## ComplexUpset (Modern Default)

**Goal:** Render a set-intersection plot with cardinality-sorted bars, optional metadata stacks (e.g., percent of intersection significant), and pre-specified queries highlighting biologically relevant intersections.

**Approach:** Convert set memberships to a long-format data frame with one row per element and binary columns per set; pass to `upset()` with `intersections='all'` or pre-specified subset; use `ComplexUpset::upset_query` to highlight intersections.

```r
library(ComplexUpset)
library(ggplot2)

# Convert from list of sets to long format
sets <- list(SetA = c('Gene1','Gene2','Gene3','Gene4'),
             SetB = c('Gene2','Gene3','Gene5','Gene6'),
             SetC = c('Gene1','Gene3','Gene6','Gene7'),
             SetD = c('Gene3','Gene4','Gene7','Gene8'))

# Long-format binary membership matrix
all_elements <- unique(unlist(sets))
df <- data.frame(element = all_elements)
for (s in names(sets)) df[[s]] <- df$element %in% sets[[s]]

# UpSet
upset(df,
      intersect = names(sets),                          # which columns are sets
      n_intersections = 20,                             # show top 20 intersections
      sort_intersections = 'descending',                # by cardinality
      sort_intersections_by = 'cardinality',            # 'cardinality' OR 'degree'
      base_annotations = list(
          'Intersection size' = intersection_size(
              counts = TRUE,
              text = list(size = 3))),
      themes = upset_modify_themes(
          list('Intersection size' = theme(panel.grid = element_blank()))))
```

## Sorting -- Cardinality vs Degree

**Cardinality sort** (default): intersections ordered by size (largest first). Reveals "the biggest overlap is A∩B."

**Degree sort**: intersections grouped by *number of sets they include* (1-set intersections, then 2-set, then 3-set, etc.). Reveals "how distributed are the overlaps across set counts?"

Choose based on the scientific question. Cardinality is the default for "find the biggest overlap"; degree is appropriate when comparing across "exclusive to 1 set" vs "shared by all."

## Pre-Specified Queries / Highlighting

```r
upset(df,
      intersect = names(sets),
      queries = list(
          upset_query(intersect = c('SetA', 'SetB'),
                       color = '#D55E00', fill = '#D55E00',
                       only_components = c('intersections_matrix', 'Intersection size')),
          upset_query(intersect = c('SetA', 'SetC', 'SetD'),
                       color = '#0072B2', fill = '#0072B2',
                       only_components = c('intersections_matrix', 'Intersection size'))))
```

## Attribute Panels (ComplexUpset Strength)

Unlike UpSetR's "boxplot.summary," ComplexUpset supports arbitrary ggplot annotations stacked above the intersection bars:

```r
upset(df,
      intersect = names(sets),
      annotations = list(
          'log2 FC' = ggplot(mapping = aes(x = intersection, y = log2FC)) +
                       geom_boxplot() + theme_classic(),
          'Significant fraction' = ggplot(mapping = aes(x = intersection, fill = significant)) +
                                    geom_bar(position = 'fill') +
                                    scale_fill_manual(values = c('TRUE' = '#D55E00', 'FALSE' = 'grey80')) +
                                    theme_classic()))
```

## upsetplot (Python)

```python
from upsetplot import from_contents, UpSet
import matplotlib.pyplot as plt

sets = {'SetA': ['Gene1','Gene2','Gene3','Gene4'],
        'SetB': ['Gene2','Gene3','Gene5','Gene6'],
        'SetC': ['Gene1','Gene3','Gene6','Gene7']}
data = from_contents(sets)

upset = UpSet(data,
              subset_size='count',
              show_counts=True,
              sort_by='cardinality',                    # 'cardinality' OR 'degree'
              sort_categories_by='cardinality',
              facecolor='#0072B2',
              element_size=40)
upset.style_subsets(present=['SetA', 'SetB'], facecolor='#D55E00')   # highlight specific intersection
fig = plt.figure(figsize=(8, 5))
upset.plot(fig=fig)
plt.savefig('upset.pdf', bbox_inches='tight')
```

## UpSetR (Legacy — Use Only for Reproducibility)

```r
library(UpSetR)
upset(fromList(sets),
      nsets = 4, nintersects = 20,
      order.by = 'freq',
      decreasing = TRUE,
      mb.ratio = c(0.6, 0.4),
      point.size = 3,
      line.size = 1,
      text.scale = c(1.5, 1.3, 1.3, 1, 1.5, 1.3))
```

UpSetR works but lacks ggplot2 grammar and active maintenance. Reproducing a paper's UpSetR figure is the main reason to use it in 2026.

## Per-Method Failure Modes

### Using UpSetR for new work in 2026

**Trigger:** Following older tutorials that default to UpSetR.

**Mechanism:** UpSetR has not had a CRAN release since 2019; integration with current ggplot2 / R ecosystem stale.

**Symptom:** Limited customization; ggplot2 layer not available; eventual breakage.

**Fix:** Switch to ComplexUpset for new figures. UpSetR is fine for reproducing old figures.

### ggplot2 4.0 broke ComplexUpset

**Trigger:** ggplot2 4.0 (mid-2025) introduced API changes; ComplexUpset's `upset()` errored.

**Mechanism:** Upstream ggplot2 changes affected ComplexUpset internals (issue #213).

**Symptom:** "Error in `upset()`: ..." after ggplot2 upgrade.

**Fix:** Pin compatible versions (`renv::install('ggplot2@3.5.2')`) until ComplexUpset patches. Check GitHub issues for fix status.

### Too many sets makes UpSet unreadable

**Trigger:** UpSet with 10+ sets and `n_intersections = Inf`.

**Mechanism:** Number of possible intersections is 2^N − 1; with 10 sets that's 1023 columns.

**Symptom:** Vertical bars too thin to read; matrix dots unrecognizable.

**Fix:** Set `n_intersections = 20` (or whatever fits); pre-filter to relevant intersections via `intersections = list(c('SetA','SetB'), c('SetA','SetC','SetD'))`.

### Single-set "intersections" obscure cross-set overlap story

**Trigger:** Default sort by cardinality puts "set exclusives" first (often largest).

**Mechanism:** "SetA only" is technically a 1-set intersection; usually larger than any 2+set overlap.

**Symptom:** First 4-5 bars are "exclusive to X," obscuring the cross-set story.

**Fix:** Filter via `intersections` argument to exclude 1-set; OR sort by degree to group; OR use `mode='intersect'` (vs `'distinct'`) for different counting.

### Element duplicate across sets in `fromList`

**Trigger:** Same element appears in multiple sets but stored as duplicate rows.

**Mechanism:** `fromList` expects each element appears once per set; duplicates inflate counts.

**Symptom:** Intersection counts don't sum to known totals.

**Fix:** `lapply(sets, unique)` before `fromList`.

### upsetplot from_contents vs from_indicators

**Trigger:** Wrong input format function used.

**Mechanism:** `from_contents` for dict of element lists; `from_indicators` for already-pivoted binary frame.

**Symptom:** TypeError or wrong intersections.

**Fix:** Check input shape; use the appropriate constructor.

## Reconciliation: When Implementations Differ

| Pattern | Cause | Action |
|---------|-------|--------|
| ComplexUpset and UpSetR show different intersection counts | Different element duplication handling | `lapply(sets, unique)`; verify both agree |
| ComplexUpset slow on >10 sets | 2^N intersections enumerated | Pre-specify relevant intersections; use `n_intersections` |
| upsetplot Python output differs from R | sort_by default differs | Set sort_by explicitly in both |
| Excluding 1-set intersections | mode='distinct' vs 'intersect' | `intersections` parameter; document |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Max sets for legible UpSet | 8-10 | Visualization practical |
| Show top intersections | 15-25 | Above this matrix too thin |
| When to use UpSet vs Venn | >3 sets | Lex 2014 |
| 2^N intersections | grows exponentially | Set n_intersections limit |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Intersection columns too thin | Too many intersections shown | `n_intersections = 20`; pre-filter |
| 1-set bars dominate | Default cardinality sort | Exclude 1-set OR sort by degree |
| Intersection counts wrong | Duplicate elements in fromList input | `lapply(sets, unique)` |
| UpSetR error after R upgrade | Unmaintained package | Switch to ComplexUpset |
| ComplexUpset breaks after ggplot2 update | ggplot2 4.0 issue #213 | Pin ggplot2 ≤ 3.5.2 |
| Python upsetplot mismatch with R | Different default sort | Standardize sort_by |

## References

- Conway JR, Lex A, Gehlenborg N. 2017. UpSetR: an R package for the visualization of intersecting sets and their properties. *Bioinformatics* 33(18):2938-2940.
- Krassowski M. 2020. ComplexUpset (R package). https://github.com/krassowski/complex-upset
- Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H. 2014. UpSet: visualization of intersecting sets. *IEEE Trans Vis Comput Graph* 20(12):1983-1992.

## Related Skills

- data-visualization/heatmaps-clustering - Alternative for smaller set membership (Venn alternative is OncoPrint-style)
- pathway-analysis/go-enrichment - Gene-set overlaps to visualize
- differential-expression/de-results - DE gene-list comparisons
- data-visualization/flow-and-transition-plots - Alluvial as alternative for membership flow
