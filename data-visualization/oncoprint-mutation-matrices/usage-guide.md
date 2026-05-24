# OncoPrint and Mutation Matrix Plots - Usage Guide

## Overview

OncoPrint is the canonical cancer-cohort mutation visualization (Cerami 2012; cBioPortal). Each cell in a gene-by-sample grid stacks colored rectangles per alteration class (missense, truncating, splice, copy-gain, copy-loss, fusion), preserving multi-modal alteration patterns that a flat heatmap would lose. ComplexHeatmap::oncoPrint is the most customizable R implementation; maftools::oncoplot is faster to onboard; comut.py covers Python.

## Prerequisites

```r
install.packages(c('circlize'))
BiocManager::install(c('ComplexHeatmap', 'maftools'))
```

```bash
pip install comut pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Build an OncoPrint of top 20 recurrently mutated genes across the cohort"
- "Add clinical annotation tracks for subtype and stage"
- "Show co-occurrence and mutual exclusivity between top driver genes"
- "Split samples by molecular subtype"
- "Log-transform the TMB bar so hypermutators don't dominate"

## Example Prompts

### TCGA-style cohort

> "OncoPrint from a TCGA MAF: top 20 recurrently mutated genes, stacked by alteration class (missense, truncating, splice, amp, homdel, fusion). Add subtype and stage clinical annotations. Default memoSort for samples."

### Subtype split

> "OncoPrint split into 3 column panels by molecular subtype. Per-panel gene frequency right bar."

### Mutex / co-occurrence

> "Run somaticInteractions (maftools) to identify mutually exclusive and co-occurring driver pairs among the top 20 genes."

### Custom alter_fun

> "Build the cells as ;-separated alteration strings; render Amp as full red rectangle, HomDel as full blue, Missense as green half-height, Truncating as black quarter-height, Fusion as triangle marker."

### Hypermutator-aware TMB

> "Add a TMB bar on top of the OncoPrint; log10-transform so 1-2 POLE-mutant samples don't saturate."

## What the Agent Will Do

1. Load MAF / variant table; map per-row Variant_Classification to alteration class.
2. Optionally integrate CNV calls (GISTIC2 or per-gene Amp/HomDel calls).
3. Pivot to gene-by-sample matrix with `;`-separated alteration strings.
4. Compute top N most-frequently-altered genes.
5. Define alter_fun list (one render function per alteration class).
6. Build HeatmapAnnotation for clinical metadata with CVD-safe palette.
7. Render with `oncoPrint(mat, alter_fun, col, top_annotation = ha)`; use default memoSort.
8. Optionally compute somaticInteractions for mutex/co-occurrence and add as side annotation.
9. Export PDF at 300 DPI; rasterize the cell layer for cohorts >500 samples.

## Tips

- **OncoPrint cells stack multiple alterations** - do NOT flatten to single class. The whole point is showing co-occurring multi-class events per gene per sample.

- **Default memoSort** (Cerami 2012) sorts samples by binary altered/not patterns across top genes, producing the canonical staircase. Override only if a clinical-group sort is needed.

- **`remove_empty_columns = FALSE`** preserves all cohort samples; otherwise sample count and percentage denominators don't reflect the cohort.

- **Mutex/co-occurrence on small cohorts is underpowered.** Need ≥100 samples for credible Fisher; DISCOVER (Canisius 2016) when mutation rate varies 100× across samples.

- **Hypermutators saturate linear TMB bars.** Log10 + 1 transform; or cap the y-axis.

- **maftools::oncoplot is faster to onboard** for TCGA MAFs; ComplexHeatmap::oncoPrint is more customizable.

- **`column_split` by clinical subtype** is the standard layout for subtype-driver-enrichment displays.

- **`alter_fun` rendering order** = list order. Foreground alterations should appear last in the list.

- **Percentage bars on the right** show per-gene altered fraction in the cohort. Document the denominator (cohort N vs altered-only N).

- **For cohorts >500 samples**, sample-name labels become unreadable; show only TMB and clinical tracks.

## Related Skills

- data-visualization/heatmaps-clustering - Generic heatmap underlying oncoPrint
- data-visualization/lollipop-protein-maps - Per-gene mutation maps on protein
- data-visualization/color-palettes - Alteration-class palette
- clinical-databases/variant-prioritization - Pre-OncoPrint variant filtering
- variant-calling/variant-annotation - Variant consequence annotation
- copy-number/cnv-annotation - Integrate CNVs into OncoPrint cells
