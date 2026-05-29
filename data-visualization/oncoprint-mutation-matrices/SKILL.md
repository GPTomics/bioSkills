---
name: bio-data-visualization-oncoprint-mutation-matrices
description: Build OncoPrint and co-mutation matrix plots from somatic-variant cohorts using ComplexHeatmap, maftools, and comut.py with alteration-type stacking, sample ordering by mutational burden, mutual-exclusivity overlays, and clinical annotation tracks. Use when visualizing per-sample mutation patterns across recurrent driver genes, comparing alteration classes, or identifying mutually-exclusive / co-occurring driver pairs.
tool_type: mixed
primary_tool: ComplexHeatmap
---

## Version Compatibility

Reference examples tested with: ComplexHeatmap 2.18+, maftools 2.18+, comut 0.0.3+, MAFtools requires R 4.0+; comut.py requires pandas 2.0+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# OncoPrint and Mutation Matrix Plots

**"Plot mutations across a cohort"** -> Render a gene-by-sample matrix where each cell stacks colored rectangles encoding alteration class (missense, truncating, splice, copy-gain, copy-loss, fusion). Sort samples by burden, optionally split by clinical group, and overlay co-mutation / mutual-exclusivity annotations. OncoPrint (Cerami 2012 *Cancer Discov* 2:401; canonical at cBioPortal) is the genre-defining visualization.

- R: `ComplexHeatmap::oncoPrint`, `maftools::oncoplot`
- Python: `comut.CoMut`, `cbioportal`-style implementations

## The Single Most Important Modern Insight -- Cell Stacking Encodes Multiple Alterations Per Cell

OncoPrint differs from a generic heatmap because each cell can encode multiple alterations simultaneously through *stacked* rectangles. A patient with both a missense and a copy-gain in TP53 shows one cell with two overlapping colored rectangles (e.g., green diamond inside red square). This stacking is the whole point — it preserves the multi-modal alteration landscape that flattening to a single category destroys.

In ComplexHeatmap's `oncoPrint`, the `alter_fun` argument is the rendering specification: a named list of functions, one per alteration class, each drawing its rectangle inside the cell. Get this right and the figure works; get it wrong and overlapping alterations are invisible.

## Decision Tree by Cohort and Question

| Question | Sort by | Display |
|----------|---------|---------|
| Which genes are most altered? | Gene frequency (default) | Bar above samples (sample TMB); bar right of genes (gene frequency) |
| Per-patient burden patterns | Sample burden | TMB bar on top; sample-name labels |
| Subtype-driver enrichment | Clinical group then burden | `column_split` by group; per-group frequency right bar |
| Mutual exclusivity (BRAF vs NRAS) | Custom (alphabetic-by-mutation pattern) | Memo sort; overlay log10(OR) heatmap |
| Co-occurrence (TP53 + MYC) | Custom | Same pattern; positive OR coloring |
| Driver vs passenger comparison | Two panels | Concatenate two oncoPrints horizontally |

## ComplexHeatmap::oncoPrint -- Canonical Implementation

**Goal:** Render a cohort mutation matrix with stacked alteration-class encoding, sample annotations, and a sample-sorted, gene-frequency-ranked layout.

**Approach:** Convert the MAF/variant table to a gene-by-sample matrix of `;`-delimited alteration strings; define `alter_fun` rendering one rectangle per class; pass to `oncoPrint()` with column annotations.

```r
library(ComplexHeatmap)
library(circlize)

# Input: matrix where each cell is a string like 'Missense;Amp' or '' for no alteration
# Rows = genes; columns = samples

# Color per alteration class
col <- c('Missense'   = '#56B4E9',
         'Truncating' = '#000000',
         'Splice'     = '#CC79A7',
         'Amp'        = '#D55E00',
         'HomDel'     = '#0072B2',
         'Fusion'     = '#009E73')

# alter_fun -- one function per class, each drawing inside the cell
alter_fun <- list(
    background = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h - unit(0.5, 'mm'),
                  gp = gpar(fill = '#EEEEEE', col = NA)),
    Amp = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h - unit(0.5, 'mm'),
                  gp = gpar(fill = col['Amp'], col = NA)),
    HomDel = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h - unit(0.5, 'mm'),
                  gp = gpar(fill = col['HomDel'], col = NA)),
    Missense = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h * 0.5,
                  gp = gpar(fill = col['Missense'], col = NA)),
    Truncating = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h * 0.33,
                  gp = gpar(fill = col['Truncating'], col = NA)),
    Splice = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h * 0.25,
                  gp = gpar(fill = col['Splice'], col = NA)),
    Fusion = function(x, y, w, h)
        grid.points(x, y, pch = 17, size = unit(2, 'mm'),
                    gp = gpar(col = col['Fusion'])))

# Clinical column annotation
ha_clin <- HeatmapAnnotation(
    Subtype = clinical$subtype,
    Stage   = clinical$stage,
    col = list(Subtype = c(Luminal='#0072B2', Basal='#D55E00', HER2='#009E73'),
               Stage   = c(I='#FFFFCC', II='#FED976', III='#FD8D3C', IV='#BD0026')))

oncoPrint(mat,
          alter_fun = alter_fun,
          col = col,
          top_annotation = ha_clin,
          column_title = 'TCGA-BRCA mutation landscape',
          row_names_gp = gpar(fontsize = 8),
          pct_gp = gpar(fontsize = 7),
          show_pct = TRUE,
          remove_empty_columns = FALSE,
          remove_empty_rows = FALSE)
```

## maftools::oncoplot -- Faster Onboarding

For TCGA-style MAF files, `maftools::oncoplot` is the lower-friction option:

```r
library(maftools)
maf <- read.maf(maf = 'tcga.maf', clinicalData = clinical)
oncoplot(maf = maf,
         top = 20,                            # top 20 mutated genes
         clinicalFeatures = c('Subtype', 'Stage'),
         annotationColor = list(Subtype = c(Luminal='#0072B2', Basal='#D55E00'),
                                 Stage = c(I='#FFFFCC', IV='#BD0026')),
         sortByAnnotation = TRUE,
         removeNonMutated = FALSE)
```

maftools defaults handle alteration-class colors, sample sorting, and percentage bars automatically. Customization is more limited than ComplexHeatmap.

## Mutual Exclusivity and Co-Occurrence

```r
# maftools provides somaticInteractions
si <- somaticInteractions(maf = maf, top = 20,
                          pvalue = c(0.05, 0.01),
                          fontSize = 0.7)
# Plot returns a matrix of -log10(p) with sign by direction (+ co-occur, - mutex)
```

Mutual-exclusivity testing on small cohorts (N < 50) is underpowered; reported "significant" mutex on n=20 with 2 mutations each is uninterpretable. Aggregate to larger cohorts (TCGA + ICGC pan-cancer) or report effect size with CI rather than p-value.

**Fisher exact vs DISCOVER:** standard 2x2 Fisher tests sample-mutation pairs, ignoring per-gene mutation rate background. DISCOVER (Canisius 2016 *Genome Biol* 17:261) models per-tumor mutation probability and is preferred for pan-cancer analyses where mutation rate varies 100× across samples.

## comut.py -- Python Equivalent

```python
import comut
import pandas as pd

# Long-format: columns = sample, category (gene), value (alteration class)
toy_comut = comut.CoMut()
toy_comut.add_categorical_data(
    data=mutation_long_df,
    name='Mutations',
    category_order=top_genes,
    value_order=['Truncating', 'Missense', 'Splice', 'Amp', 'HomDel'],
    mapping={'Truncating': '#000000', 'Missense': '#56B4E9',
             'Splice': '#CC79A7', 'Amp': '#D55E00', 'HomDel': '#0072B2'})

toy_comut.add_categorical_data(
    data=clinical_long_df,
    name='Subtype',
    mapping={'Luminal': '#0072B2', 'Basal': '#D55E00'})

toy_comut.add_continuous_data(
    data=tmb_long_df,
    name='TMB',
    mapping='viridis',
    value_range=(0, 30))

toy_comut.plot_comut(figsize=(12, 8))
toy_comut.figure.savefig('comut.pdf', dpi=300, bbox_inches='tight')
```

## Per-Method Failure Modes

### Alterations flattened to a single class

**Trigger:** Reducing each cell to a single most-severe alteration, losing the stack.

**Mechanism:** Loses the multi-alteration biology (e.g., MYC amp + missense in TP53).

**Symptom:** OncoPrint looks like a simple heatmap; co-occurring multi-class events invisible.

**Fix:** Build the cell as `;`-separated alteration string; define `alter_fun` for each class.

### Sample sort by gene 1 frequency only

**Trigger:** Default `oncoPrint` sorts samples by altered-gene-1 status; weakens the "memo sort" pattern.

**Mechanism:** True OncoPrint uses memoSort (Cerami 2012) which sorts by the binary altered-or-not pattern across the top genes.

**Symptom:** Samples with the same alteration profile are not adjacent; "staircase" pattern lost.

**Fix:** ComplexHeatmap `oncoPrint` uses memoSort by default; do NOT override `column_order` unless intentional.

### Showing only mutated samples (`remove_empty_columns = TRUE`)

**Trigger:** Default in some implementations.

**Mechanism:** Drops samples with no mutations in the displayed genes — but those samples ARE part of the cohort.

**Symptom:** Sample count differs from cohort N; denominator-based percentages wrong.

**Fix:** `remove_empty_columns = FALSE` to preserve all samples; percentages now reflect true cohort fraction.

### Hypermutators dominate visual

**Trigger:** Cohort with 1-2 POLE-mutant or MSI-H samples; TMB bar saturates.

**Mechanism:** Hypermutator TMB is 10-100× the typical sample.

**Symptom:** All other samples' TMB bars are invisible; one column dominates.

**Fix:** Log-transform the TMB annotation: `anno_barplot(log10(tmb + 1))`; OR cap with `ylim`.

### Mutex/co-occurrence p-values overinterpreted on small cohorts

**Trigger:** Fisher exact test on N < 50 with low mutation counts.

**Mechanism:** With 2 mutations vs 3 mutations in 20 samples, all p-values are dominated by noise.

**Symptom:** "Significant mutex" claim from a tiny pilot.

**Fix:** Aggregate to ≥100 samples for credible mutex; use DISCOVER (Canisius 2016) instead of Fisher when mutation rate varies 100× across samples.

## Small-Cohort Regime (N = 20-50)

For rare-cancer cohorts where N < 50, the standard OncoPrint + Fisher mutex pipeline is statistically uninterpretable:

| Action | What to do |
|--------|------------|
| Report per-gene frequencies | Use exact-binomial CI (Clopper-Pearson via `binom.test`) — Wald CI is invalid at low frequency |
| Do NOT report mutex p-values | Fisher exact on 2x2 with cell counts ≤ 5 has no power; the "significant" mutex finding is noise |
| Hypothesis generation only | Pool with TCGA Pan-Cancer + ICGC for credible mutex; treat the cohort as the *replication* not the discovery |
| Co-occurrence reporting | OR with Haldane-Anscombe 0.5 correction for zero cells; report alongside cohort N |

Show the OncoPrint for visual transparency, but the per-gene-frequency *table* (with exact-binomial CIs) is the load-bearing scientific output, not the mutex test.

## Reconciliation: When Implementations Differ

| Pattern | Cause | Action |
|---------|-------|--------|
| ComplexHeatmap and maftools show different sample orders | Different memoSort defaults | Specify `sortByAnnotation` explicitly; report sort criterion in caption |
| Percentage labels differ | `remove_empty_columns = TRUE` vs FALSE | Document denominator (cohort-N vs altered-N) |
| Some alterations missing from a sample | Filtering: silent SNVs, low VAF | Document filtering criteria upstream |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Cohort N for valid mutex | ≥100 (pan-cancer); ≥50 (single-cohort with effect-size focus) | Common practice |
| Display top genes | 10-25 in single panel | More creates visual clutter |
| Sample N for OncoPrint | 50-1000 (above: switch to summary panel) | Visualization practical |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Co-occurring multi-class events invisible | Single-class flattening | Use `;`-separated cells + alter_fun list |
| Sample order doesn't show staircase | column_order override | Trust default memoSort |
| Sample count differs from cohort | `remove_empty_columns = TRUE` | Set to FALSE |
| TMB bar dominated by 1-2 samples | Hypermutators on linear scale | log10 + 1 transform |
| Mutex p-values on N=20 | Underpowered | Aggregate cohorts; use DISCOVER |
| Gene frequency right-bar mismatches percentages | Denominator definition | Document cohort-N vs altered-N |

## References

- Canisius S, Martens JWM, Wessels LFA. 2016. A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence. *Genome Biol* 17:261.
- Cerami E, Gao J, Dogrusoz U, et al. 2012. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov* 2(5):401-404.
- Gao J, Aksoy BA, Dogrusoz U, et al. 2013. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. *Sci Signal* 6(269):pl1.
- Gu Z, Eils R, Schlesner M. 2016. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics* 32(18):2847-2849.
- Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Res* 28(11):1747-1756.

## Related Skills

- data-visualization/heatmaps-clustering - Generic heatmap underlying oncoPrint
- data-visualization/lollipop-protein-maps - Per-gene mutation maps on protein domains
- data-visualization/color-palettes - Alteration-class palette selection
- clinical-databases/variant-prioritization - Filter variants before OncoPrint
- variant-calling/variant-annotation - Annotate consequences upstream
- copy-number/cnv-annotation - Integrate CNV calls into the oncoprint
