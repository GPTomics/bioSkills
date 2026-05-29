---
name: bio-data-visualization-volcano-and-ma-plots
description: Build volcano and MA plots from differential-expression / association results with LFC shrinkage, FDR-adjusted thresholds, sensible label placement, and axis-truncation conventions. Covers EnhancedVolcano, ggplot2, matplotlib, and the apeglm/ashr/normal shrinkage decision. Use when visualizing differential-expression results (RNA-seq, ChIP-seq, ATAC-seq, proteomics) or any per-feature effect-size + p-value table.
tool_type: mixed
primary_tool: ggplot2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, EnhancedVolcano 1.20+, ggplot2 3.5+, ggrepel 0.9.5+, matplotlib 3.8+, numpy 1.26+, adjustText 1.1+, apeglm 1.28+, ashr 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Volcano and MA Plots

**"Plot differential-expression results"** -> Place per-feature shrunken effect estimate on the x-axis and a significance measure (-log10 padj or -log10 p) on the y-axis. The decision space spans which effect estimate (raw vs shrunken), which significance measure (raw p vs adjusted vs s-value), how to encode categories (color by direction, not by gradient), how to label (top-N is rarely informative), and how to handle the tail (extreme p compresses the plot).

- R: `EnhancedVolcano::EnhancedVolcano()`, `ggplot2 + ggrepel`, `DESeq2::plotMA()`
- Python: `matplotlib.scatter` with `adjustText`, `sanbomics.tools.volcano`, custom seaborn

## The Single Most Important Modern Insight -- Plot Shrunken LFC

Raw log2 fold change from DESeq2 / edgeR is the maximum-likelihood estimate and **inflates wildly at low counts**. A gene with 2 vs 0 reads gets log2FC = Inf; one with 4 vs 1 gets log2FC = 2 with a huge standard error. A naive volcano labels these as "top hits" purely because they have extreme estimates, not because they have a real signal.

`DESeq2::lfcShrink()` applies an empirical-Bayes prior to pull noisy low-count LFCs toward zero while leaving well-estimated genes essentially untouched. Since DESeq2 v1.28, the default shrinkage is `type='apeglm'` (Zhu, Ibrahim, Love 2019 *Bioinformatics* 35:2084) which uses a Cauchy prior — heavy enough to preserve large real effects, sharp enough at zero to deflate noise. Plot the shrunken LFC. The unshrunken LFC is a misleading effect estimate for ranking, labeling, or thresholding.

A complementary modern alternative is the **s-value** (Stephens 2017 *Biostatistics* 18:275): the local false sign rate — probability that the sign of the effect is wrong. s-values rank genes by *how confident the sign is*, which is what a volcano plot is implicitly trying to communicate. Where padj answers "is the effect non-zero," s answers "do we know which direction."

## Shrinkage Method Selection

| Method | Prior | Best for | Fails when |
|--------|-------|----------|------------|
| `apeglm` (Zhu 2019) | Cauchy | Default; preserves large LFCs, deflates noise | Requires `coef=`; no support for `contrast=` |
| `ashr` (Stephens 2017) | Mixture of normals | Comparisons requiring `contrast=`; supports s-values | Slightly more aggressive shrinkage of medium effects |
| `normal` (DESeq2 original) | Zero-centered normal | Legacy reproducibility only | Over-shrinks large effects; **deprecated** in current DESeq2 vignette |
| Unshrunken MLE | None | NEVER for volcano/MA plots | Low-count genes dominate the tails with no real signal |

```r
library(DESeq2)
dds <- DESeq(dds)
res_apeglm <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
res_ashr <- lfcShrink(dds, contrast = c('condition', 'treated', 'control'), type = 'ashr')
# ashr also returns svalue column (Stephens 2017 local false sign rate)
```

For edgeR users: `topTags()` already provides moderated p-values but does not shrink LFC. Use `glmTreat()` for a moderated test against a non-zero LFC threshold; this is the edgeR equivalent of the shrunken-LFC philosophy.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bulk RNA-seq with DESeq2 | `lfcShrink(type='apeglm')` + plot on padj threshold | apeglm is the default since DESeq2 v1.28 |
| Non-default contrast required | `lfcShrink(type='ashr')` | apeglm requires `coef=`, ashr accepts `contrast=` |
| Single-cell pseudobulk DE | `lfcShrink(type='apeglm')` + raster the plot | sc datasets often >20k genes; raster prevents PDF lag |
| Proteomics (limma/MSstats) | Already moderated; plot logFC vs adj.P.Val directly | limma's empirical-Bayes already shrinks |
| Microarray (limma) | `topTable()` adj.P.Val + logFC | Same — limma is the original shrunken-LFC method |
| ATAC/ChIP differential peaks | `lfcShrink(type='apeglm')` on DESeq2/DiffBind | Treat peaks as features identically to genes |
| Want to rank by sign-confidence | ashr's svalue, NOT padj | Stephens 2017 — sign-aware ranking |
| Many comparisons in one figure | Faceted MA plot with shared y-axis | MA scales better than volcano for >6 panels |

## Volcano with ggplot2 + ggrepel

**Goal:** Plot shrunken LFC vs -log10 p-value, color by significance class, label genes of interest with non-overlapping repulsion.

**Approach:** Compute a categorical significance variable from padj AND |LFC| thresholds; pre-select labels (genes of interest OR top-N by combined rank) before plotting; use `ggrepel::geom_text_repel` with `max.overlaps = Inf` to guarantee every selected label appears.

```r
library(ggplot2)
library(ggrepel)
library(dplyr)

volcano_plot <- function(res, fdr = 0.05, lfc_threshold = 1, label_genes = NULL, top_n = 10) {
    res <- as.data.frame(res) %>%
        tibble::rownames_to_column('gene') %>%
        mutate(
            significance = case_when(
                is.na(padj) ~ 'NS',
                padj < fdr & log2FoldChange > lfc_threshold ~ 'Up',
                padj < fdr & log2FoldChange < -lfc_threshold ~ 'Down',
                TRUE ~ 'NS'
            ),
            neg_log10_p = -log10(pvalue)
        )

    if (is.null(label_genes)) {
        label_genes <- res %>%
            filter(significance != 'NS') %>%
            mutate(rank_score = -log10(pvalue) * abs(log2FoldChange)) %>%
            arrange(desc(rank_score)) %>%
            head(top_n) %>%
            pull(gene)
    }
    res$label <- ifelse(res$gene %in% label_genes, res$gene, '')

    okabe_ito <- c(Up = '#D55E00', Down = '#0072B2', NS = '#999999')

    ggplot(res, aes(log2FoldChange, neg_log10_p, color = significance)) +
        geom_point(alpha = 0.6, size = 1.3) +
        scale_color_manual(values = okabe_ito, name = NULL) +
        geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                   linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
        geom_hline(yintercept = -log10(fdr), linetype = 'dashed',
                   color = 'grey40', linewidth = 0.3) +
        geom_text_repel(aes(label = label), color = 'black', size = 3,
                        max.overlaps = Inf, box.padding = 0.4, segment.size = 0.2,
                        min.segment.length = 0) +
        labs(x = expression(log[2]~'fold change (shrunken)'),
             y = expression(-log[10]~italic(p))) +
        theme_classic(base_size = 10) +
        theme(panel.grid = element_blank())
}
```

Key design choices encoded above:
- Colors from Okabe-Ito (Wong 2011 *Nat Methods* 8:441) — CVD-safe categorical palette
- `max.overlaps = Inf` because the ggrepel default of 10 silently drops labels with no error (see [[api_gotchas]])
- Threshold line on `-log10(fdr)` matches the *adjusted* p threshold; drawing it on raw p creates a meaningless line
- `rank_score = -log10(pvalue) * abs(log2FoldChange)` selects labels that are both significant AND have non-trivial effect; pure top-N-by-p selects high-count genes with tiny effects

## EnhancedVolcano -- Production Use and Its Gotchas

```r
library(EnhancedVolcano)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',                    # use padj NOT pvalue for the threshold line
    pCutoff = 0.05,
    FCcutoff = 1,
    selectLab = c('TP53', 'MYC', 'BRCA1'),
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    maxoverlapsConnectors = Inf,
    colAlpha = 0.6,
    pointSize = 1.5,
    labSize = 3,
    col = c('grey60', '#0072B2', '#56B4E9', '#D55E00'),
    legendPosition = 'right')
```

**Gotcha 1: `selectLab` filters through pCutoff AND FCcutoff.** Genes explicitly listed but failing thresholds appear unlabeled, with no warning. This is the most common "why is the gene missing" failure. To force-label regardless of thresholds, pre-shrink the input so the genes pass the cutoff, or switch to a manual ggplot layer.

**Gotcha 2: `y = 'pvalue'` vs `y = 'padj'`.** Many tutorials use raw `pvalue` for the y-axis, then draw the threshold line at `-log10(0.05)` — that line corresponds to a raw p < 0.05, not an FDR. Use `y = 'padj'` so the threshold is meaningful.

**Gotcha 3: Asymmetric x-limits hide the diverging null distribution.** Use `xlim = c(-max(abs(LFC)), max(abs(LFC)))` so the plot is symmetric around zero.

## MA Plot -- The Underused Diagnostic

The MA plot (log2-mean vs log2-fold-change) is the original RNA-seq diagnostic (Dudoit 2002 *JASA*). It exposes the abundance-dependent variance structure that the volcano hides:

- A fan-shaped MA plot with extreme LFCs concentrated at low baseMean indicates **inadequate shrinkage**
- A horizontal "stripe" of significant genes at one LFC value indicates **batch confound with treatment**
- An asymmetric distribution (more Up than Down) at the low-count end indicates **library-size normalization failure**

```r
library(DESeq2)
plotMA(res_apeglm, alpha = 0.05, ylim = c(-5, 5))
# alpha colors significant points; ylim clips for readability without losing the gene
```

```python
import matplotlib.pyplot as plt
import numpy as np

def ma_plot(res, fdr=0.05, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    sig = (res['padj'] < fdr) & res['padj'].notna()
    ax.scatter(np.log10(res.loc[~sig, 'baseMean']), res.loc[~sig, 'log2FoldChange'],
               c='#999999', s=4, alpha=0.4, rasterized=True)
    ax.scatter(np.log10(res.loc[sig, 'baseMean']), res.loc[sig, 'log2FoldChange'],
               c='#D55E00', s=6, alpha=0.7, rasterized=True)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xlabel(r'$\log_{10}$ mean normalized count')
    ax.set_ylabel(r'$\log_2$ fold change (shrunken)')
    return ax
```

`rasterized=True` is critical for >5000 points — vector scatter creates 5MB+ PDFs that crash Illustrator.

## Per-Method Failure Modes

### Unshrunken LFC plotted as volcano

**Trigger:** Calling `results(dds)` and plotting `log2FoldChange` directly, without `lfcShrink()`.

**Mechanism:** ML estimate has infinite-variance tails at low counts; one read difference produces log2FC = Inf.

**Symptom:** "Top hits" by |LFC| are all genes with baseMean < 5; biologically interesting genes with moderate LFC are hidden in the noise cloud.

**Fix:** `lfcShrink(dds, coef=..., type='apeglm')` for the default case; `type='ashr'` if `contrast=` is needed.

### Raw p threshold line drawn on adjusted axis

**Trigger:** Drawing `geom_hline(yintercept = -log10(0.05))` on a plot whose y-axis is `-log10(padj)`.

**Mechanism:** padj < 0.05 corresponds to FDR < 5% control, NOT to raw p < 0.05. The drawn line is at the wrong y-value relative to the data.

**Symptom:** Visible "significant" points sit below the FDR line; the legend says "FDR < 0.05" but the line doesn't separate them correctly.

**Fix:** Be explicit: if y-axis is padj, use `-log10(fdr)` as the threshold line value AND label it "FDR threshold." If y-axis is raw p, an FDR threshold cannot be drawn as a horizontal line — the FDR threshold moves per gene.

### Top-N-by-p selects low-effect-size hits

**Trigger:** `head(arrange(res, pvalue), 20)` to choose labels.

**Mechanism:** With large N, the smallest p-values belong to high-count, low-variance, biologically-boring genes (housekeeping). Effect size and statistical confidence are not the same thing.

**Symptom:** Labels are GAPDH, ACTB, B2M — never the gene that drives the biology.

**Fix:** Rank by `-log10(p) * abs(log2FoldChange)` (geometric average of the two axes), OR pre-specify labels of interest from prior knowledge.

### ggrepel `max.overlaps` silently drops labels

**Trigger:** Default `max.overlaps = 10`; 30 genes labeled in code; only 10 render.

**Mechanism:** ggrepel emits a warning ("18 unlabeled data points (too many overlaps)") but no error. In a Quarto/Rmd render the warning is buried in the log.

**Symptom:** Reviewer asks "where is gene X?"; the label was specified in code but did not render.

**Fix:** `geom_text_repel(..., max.overlaps = Inf)` or `options(ggrepel.max.overlaps = Inf)` at the top of the script.

### Extreme p-values compress the upper axis

**Trigger:** Genes with p = 1e-200 or smaller (common in cancer datasets) push the y-axis maximum to 200; all biologically meaningful genes pile up at the bottom.

**Mechanism:** -log10 expands the tail; one ultra-significant gene visually dominates.

**Symptom:** Volcano looks like an Eiffel Tower with most genes squished near y = 0-20.

**Fix:** Cap the y-axis (`ylim = c(0, 50)` and use `coord_cartesian` so capped points stay in the data but render at the edge) OR transform with `sqrt(-log10(p))` to compress the tail OR split the y-axis with `ggbreak`.

### `lfcShrink(type='normal')` on a modern DESeq2

**Trigger:** Following old tutorials that pre-date DESeq2 v1.28 when apeglm became the default.

**Mechanism:** The `'normal'` prior over-shrinks large real effects toward zero.

**Symptom:** Volcano looks "too clean" — genuine 8-fold changes appear as 2-3 fold.

**Fix:** Use `type='apeglm'` (default since v1.28) or `type='ashr'`. Vignette removed `'normal'` from recommendations.

### EnhancedVolcano's `selectLab` filters by thresholds

**Trigger:** Listing genes in `selectLab` that have `padj > pCutoff`.

**Mechanism:** Source code applies `pCutoff` AND `FCcutoff` filter to `selectLab` membership; silently drops any that don't pass.

**Symptom:** Specific genes requested in `selectLab` do not appear; no warning.

**Fix:** Pre-shrink (so genes pass) or build the labeled subset manually with `ggrepel` and add as an annotation layer.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| apeglm and ashr give different "top hits" | Different prior shapes; medium-effect, medium-count genes are most sensitive | Both are valid; pick one and document. apeglm is the DESeq2 default and the published recommendation |
| EnhancedVolcano shows fewer points than ggplot | EnhancedVolcano drops `padj = NA` (DESeq2 independent filtering) | Confirm by counting `is.na(res$padj)`; to include them, set NA padj to 1 before plotting |
| Volcano has many "significant" genes but MA plot shows them all at low baseMean | Unshrunken LFC; the volcano is showing fold-change noise | Re-plot with shrunken LFC; the MA-plot fan is the diagnostic |
| Forest of horizontal stripes in MA at integer LFC | Pseudocount-induced quantization in low-count genes | Increase normalization-method aggressiveness OR filter low-count genes upstream |
| Half the genes have padj = NA | DESeq2 independent filtering (Bourgon-Gentleman-Huber 2010 *PNAS*) excluded them as low-mean | This is correct behavior; do NOT set `independentFiltering = FALSE` to hide it. Report the NA count |

**Operational rule:** the volcano plot is read in this priority order — (1) is LFC shrunken? (2) is the y-axis padj or raw p? (3) does the threshold line match the axis? (4) are labels selected by combined rank? If any of these is wrong, the plot is misleading regardless of how attractive it looks.

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Default LFC cutoff for "biologically relevant" | \|log2FC\| > 1 (2-fold) | Convention; sensitive analyses use 0.58 (1.5-fold) for subtle effects |
| Default FDR cutoff | padj < 0.05 | Benjamini-Hochberg 1995 *JRSS-B* 57:289 |
| Stricter cutoff for unbiased screens | padj < 0.01 | Reduces false positives in unbiased genome-wide analyses |
| Relaxed cutoff for exploratory / hypothesis-generating | padj < 0.10 or 0.20 | Acceptable for follow-up enrichment, NOT for "hits" |
| s-value cutoff (Stephens 2017) | s < 0.005 corresponds approximately to padj < 0.05 | Stephens 2017 *Biostatistics* 18:275 |
| Raster threshold | >5000 points | Vector PDF crashes Illustrator; raster scatter, keep axes vector |
| ggrepel max.overlaps | Set to Inf for publication | Default 10 silently drops labels |
| Volcano y-axis cap | -log10(p) > 50 typically warrants capping | Visual compression of biologically meaningful genes |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| All "top hits" have baseMean < 5 | Unshrunken LFC | `lfcShrink(type='apeglm')` |
| Threshold line doesn't separate colored from grey points | y = pvalue but threshold drawn at FDR | Switch y to padj OR redraw line at FDR-equivalent p |
| Labeled gene does not appear in EnhancedVolcano | `selectLab` filtered by pCutoff/FCcutoff | Pre-shrink or build labels manually |
| Volcano renders as a flat horizontal cloud | Extreme p (e.g., 1e-200) dominates y-axis | Cap with `coord_cartesian(ylim = c(0, 50))` or use sqrt transform |
| PDF crashes Illustrator | Vector scatter of >10000 points | Set `rasterized = TRUE` in geom_point or matplotlib scatter |
| ggrepel labels 10 of 30 selected genes | Default `max.overlaps = 10` | `geom_text_repel(max.overlaps = Inf)` |
| Volcano "significant" gene count differs from DESeq2 summary | EnhancedVolcano drops `padj = NA` | Set NA padj to 1 or document the discrepancy |
| Up and Down counts asymmetric for a balanced experiment | Library-size normalization failure | Re-run DESeq2 with `estimateSizeFactors(type='poscounts')` for sparse data |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why is this LFC shrunken? Show me the unshrunken." | Shrunken LFC is the recommended estimate for ranking and visualization (Zhu 2019). Unshrunken LFC inflates at low counts and gives misleading rank. Unshrunken is available in the supplementary table |
| "Why padj < 0.05 not p < 0.05?" | padj controls FDR via Benjamini-Hochberg. Raw p < 0.05 across 20000 genes yields ~1000 false positives by chance; padj < 0.05 caps the expected false-positive rate at 5% of called hits |
| "Why are X gene and Y gene not labeled?" | Labels selected by combined rank (-log10(p) * \|LFC\|) or pre-specified gene list. List of all significant genes is in supplementary table T1 |
| "The volcano looks too clean / too sparse." | Color encodes 3 categories (Up/Down/NS), not a gradient. Gradient encoding implies a continuous interpretation of significance which is invalid — significance is a threshold decision |
| "Why is the x-axis asymmetric?" | Asymmetric x-axis reflects the asymmetry of the data. If symmetry is preferred for visual interpretation, use `xlim = c(-X, X)` with X = max(\|LFC\|) |

## References

- Anders S, Huber W. 2010. Differential expression analysis for sequence count data. *Genome Biol* 11:R106.
- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS-B* 57:289-300.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107:9546-9551.
- Dudoit S, Yang YH, Callow MJ, Speed TP. 2002. Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. *Stat Sin* 12:111-139.
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550.
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294. doi:10.1093/biostatistics/kxw041
- Wong B. 2011. Points of view: Color blindness. *Nat Methods* 8(6):441. doi:10.1038/nmeth.1618
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092. doi:10.1093/bioinformatics/bty895

## Related Skills

- differential-expression/de-results - Filter and rank DE result tables before plotting
- differential-expression/deseq2-basics - Run DESeq2 to produce the input results object
- differential-expression/de-visualization - DESeq2 / edgeR built-in plot helpers
- data-visualization/distribution-plots - Boxplot / raincloud follow-up for specific gene panels
- data-visualization/color-palettes - Okabe-Ito and CVD-safe palette selection
- data-visualization/ggplot2-fundamentals - Underlying grammar of graphics
- pathway-analysis/go-enrichment - Functional enrichment from the gene lists produced
