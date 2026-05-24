# UpSet Plots - Usage Guide

## Overview

UpSet plots replace Venn diagrams when comparing more than 4 sets. The visualization (Lex 2014) shows a matrix of dots indicating which sets participate in each intersection, with bars above showing intersection size and bars on the left showing per-set total size. ComplexUpset (Krassowski) is the actively-maintained ggplot2-grammar implementation and the modern default; UpSetR (Conway 2017) is the original but has had no CRAN release since 2019.

## Prerequisites

```r
install.packages('ComplexUpset')
# Legacy: install.packages('UpSetR')
```

```bash
pip install upsetplot
```

## Quick Start

Tell your AI agent what you want to do:
- "Make an UpSet plot of overlap across 6 gene sets, sorted by intersection size"
- "Add a boxplot of log2FC stacked above the intersection bars"
- "Highlight the SetA∩SetB intersection in red"
- "Sort intersections by degree (number of sets) instead of cardinality"
- "Migrate this UpSetR figure to ComplexUpset"

## Example Prompts

### Standard ComplexUpset

> "Build an UpSet plot from 5 gene sets using ComplexUpset. Top 20 intersections by cardinality. Add intersection counts above each bar."

### With metadata stack

> "UpSet plus a stacked-bar annotation showing percent significant per intersection (significant column in input)."

### Highlight queries

> "Highlight the SetA-and-SetB intersection with #D55E00 and the SetA-SetC-SetD intersection with #0072B2 using upset_query."

### Migration from UpSetR

> "Convert this UpSetR figure to ComplexUpset, preserving sort order and color scheme."

### Python upsetplot

> "upsetplot.UpSet from a dict of contents, sort by cardinality, show counts; export PDF."

## What the Agent Will Do

1. Convert input from list-of-sets to long-format binary membership frame (ComplexUpset) or counts series (upsetplot).
2. Deduplicate elements per set; verify totals match expectations.
3. Choose sort criterion: cardinality (largest overlap first) or degree (grouped by set count).
4. Filter to top N intersections (default 20-25) to keep matrix readable.
5. Configure intersection-size bar styling and color.
6. Add metadata annotations via `annotations = list(...)` (ComplexUpset) or `add_stacked_bars` (upsetplot).
7. Apply pre-specified queries to highlight intersections of interest.
8. Export at 300 DPI with cairo_pdf (R) or Type-42 fonts (Python).

## Tips

- **ComplexUpset is the modern default.** UpSetR is effectively unmaintained (no CRAN release since 2019).

- **ggplot2 4.0 (mid-2025) broke ComplexUpset's `upset()`** (issue #213). Pin ggplot2 ≤ 3.5.2 until ComplexUpset patches.

- **Cap intersections at 20-25.** Above this, vertical bars too thin to read. Use `n_intersections = 20`.

- **Pre-filter via `intersections = list(...)`** for specific multi-set overlaps of interest; avoids showing every 2^N − 1.

- **Cardinality sort vs degree sort** answer different questions. Cardinality = "biggest overlap"; degree = "grouped by set count."

- **`lapply(sets, unique)`** before fromList - duplicates inflate counts and silently mislead.

- **Single-set intersections often dominate cardinality sort.** Either exclude (`intersections` parameter) or use degree sort to relegate them.

- **ComplexUpset annotations are ggplot layers** - full power of ggplot grammar above the intersection bars. boxplot, stacked bar, jitter, scatter all work.

- **upset_query() for highlighting** specific intersections; `only_components` limits which panels show the highlight.

- **For Python**, `upsetplot.from_contents(dict)` for dict-of-lists input; `from_memberships(list_of_lists)` for membership lists.

## Related Skills

- data-visualization/heatmaps-clustering - OncoPrint as alternative
- pathway-analysis/go-enrichment - Gene-set overlaps
- differential-expression/de-results - DE list comparisons
- data-visualization/flow-and-transition-plots - Alluvial alternative
