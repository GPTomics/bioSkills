# Volcano and MA Plots - Usage Guide

## Overview

Volcano and MA plots are the standard visualizations of per-feature differential analysis (RNA-seq, ChIP-seq, ATAC-seq, proteomics, methylation). The volcano places shrunken log2 fold change on x and significance on y; the MA plot replaces the significance axis with mean abundance to expose count-dependent variance. The decision space spans which effect estimate to use (raw vs shrunken), which significance measure (raw p vs adjusted vs s-value), how to color categories, and how to handle the extreme tail of p-values.

## Prerequisites

```r
install.packages(c('ggplot2', 'ggrepel', 'dplyr'))
BiocManager::install(c('DESeq2', 'EnhancedVolcano', 'apeglm', 'ashr'))
```

```bash
pip install matplotlib numpy pandas adjustText
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot a volcano for my DESeq2 results with shrunken LFC and FDR < 0.05"
- "Make an MA plot diagnosing why my low-count genes have extreme fold changes"
- "Label TP53, MYC, BRCA1 on the volcano and the top 10 hits by combined rank"
- "Cap the y-axis at -log10(p) = 50 because one gene has p ~ 1e-200"
- "Replace the raw-p threshold line with a padj = 0.05 line"

## Example Prompts

### Volcano with shrunken LFC

> "Generate a publication-quality volcano plot from a DESeq2 results table. Apply apeglm shrinkage. Color significant up/down with Okabe-Ito palette. Label the top 10 hits by combined rank (-log10 p × |LFC|) plus a list of genes of interest if provided."

### MA plot as diagnostic

> "Plot an MA diagnostic from my DESeq2 results. Show whether extreme LFCs are concentrated at low baseMean (i.e. inadequate shrinkage)."

### Handling extreme p-values

> "Volcano with y-axis capped at -log10(p) = 50; points beyond the cap should remain visible at the edge."

### Single-cell pseudobulk DE

> "Volcano plot for pseudobulk DESeq2 of cluster vs rest with >20000 features; rasterize the scatter for PDF size."

### Cross-tool comparison

> "Run lfcShrink with apeglm and ashr and overlay the two volcanos. Identify genes where shrinkage materially changes the rank."

## What the Agent Will Do

1. Load the DESeq2 / edgeR / limma / MSstats results object, or a generic feature × {pvalue, padj, log2FoldChange} table.
2. Apply LFC shrinkage if the input is DESeq2 results and shrinkage has not already been applied - apeglm by default, ashr if a `contrast=` is needed.
3. Compute the categorical significance variable from padj AND |LFC| thresholds.
4. Select labels: explicit gene list if provided, otherwise top-N by combined rank.
5. Render with ggplot2 + ggrepel (R) or matplotlib + adjustText (Python). Use Okabe-Ito for Up/Down/NS categories.
6. Draw threshold lines at the actual padj-equivalent y-value and `±log2(FC threshold)` x-values.
7. If extreme p-values dominate the y-axis, cap with `coord_cartesian` or `sqrt` transform.
8. Rasterize the point layer for >5000 features so PDF stays under ~2 MB.
9. Export at 300 DPI with `pdf.fonttype=42` (matplotlib) or `cairo_pdf` (ggsave) for journal compliance.

## Tips

- **Plot shrunken LFC, not raw LFC.** Raw LFC at low counts is meaningless. apeglm (Zhu 2019) is the DESeq2 default since v1.28 and preserves large real effects while deflating noise. Run `lfcShrink(dds, coef=..., type='apeglm')` before plotting.

- **Use padj as the y-axis, not raw p.** A threshold line at `-log10(0.05)` on a raw-p axis crosses many non-FDR-significant points and looks identical to a line on padj - so the plot is interpreted as FDR-controlled when it isn't.

- **Use s-values for sign-aware ranking.** Stephens 2017 introduced s-value (local false sign rate) as a stronger ranking criterion than padj for the question a volcano implicitly asks: do we know the sign? Use `ashr` to compute s-values.

- **Color categorically, not by gradient.** Up/Down/NS is a 3-category decision. A continuous color gradient implies that p-value is interpretable on a continuous scale - but the threshold decision is binary by design.

- **Label by combined rank, not by smallest p.** With large N, the smallest p-values are housekeeping genes with tiny effects. Use `rank_score = -log10(p) * abs(log2FoldChange)` or pre-specify labels.

- **Cap extreme y-values.** Cancer datasets routinely produce p < 1e-200. Use `coord_cartesian(ylim = c(0, 50))` (keeps points, just clips the display) or `sqrt(-log10(p))` transform.

- **Set `ggrepel::max.overlaps = Inf`** or `options(ggrepel.max.overlaps = Inf)`. The default of 10 silently drops labels and only emits a warning - invisible in non-interactive renders.

- **`EnhancedVolcano::selectLab` filters by thresholds.** Listed genes that fail pCutoff or FCcutoff are silently unlabeled. Use a manual ggrepel layer when this matters.

- **Rasterize the point layer above ~5000 features.** Vector scatter at this scale crashes Illustrator and bloats PDFs. Axis labels and threshold lines should remain vector.

- **MA plots diagnose shrinkage failure.** A fan-shaped MA with extreme |LFC| concentrated at low baseMean means LFC is unshrunken. A horizontal stripe of significant genes means batch is confounded with treatment.

- **Independent filtering produces `padj = NA`.** DESeq2 sets padj to NA for low-mean genes (Bourgon-Gentleman-Huber 2010). EnhancedVolcano drops them silently. Confirm the count via `sum(is.na(res$padj))` if the displayed point count differs from the DESeq2 summary.

## Related Skills

- differential-expression/de-results - Filter and rank DE results before plotting
- differential-expression/deseq2-basics - Run DESeq2 to produce the input
- data-visualization/distribution-plots - Boxplot / raincloud follow-up for individual genes
- data-visualization/color-palettes - Okabe-Ito and CVD-safe palettes
- pathway-analysis/go-enrichment - Functional enrichment from the gene lists
