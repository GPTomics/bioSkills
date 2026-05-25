---
name: bio-differential-expression-de-visualization
description: Creates DE-specific diagnostic and result visualizations using DESeq2/edgeR built-in functions and lightweight ggplot2 wrappers. Covers MA plot (with the shrunken-LFC compression effect), volcano (with the apeglm caveat that p-values are unchanged), PCA on VST/rlog (never raw counts), sample distance heatmaps, top-DE-gene heatmaps with the row-scaling trap, dispersion / BCV plot interpretation, p-value histogram diagnostics, plotCounts for individual genes, blind=TRUE vs FALSE rationale, and the n=3 visualization stake. Use when generating DE diagnostic plots, choosing VST vs rlog for visualization, troubleshooting suspicious plot patterns (shifted MA cloud, batch-dominated PCA, anti-conservative p-value histogram), or building a standard QC figure panel.
tool_type: r
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, limma 3.58+, ggplot2 3.5+, pheatmap 1.0+, RColorBrewer 1.1+, ggrepel 0.9+, EnhancedVolcano 1.20+, matrixStats 1.2+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# DE Visualization

**"Make the standard DE figure panel"** -> Use built-in functions or thin wrappers to produce diagnostic plots (dispersion, p-value histogram, PCA, sample distance) and result plots (MA, volcano, heatmap of top DE genes, per-gene counts), interpreted as diagnostics of the underlying model.

## Scope

This skill covers DE-specific built-in plots and immediate wrappers. For richer customization:
- Custom volcano/MA with apeglm-shrunken LFC and ggrepel labelling -> `data-visualization/volcano-and-ma-plots`
- PCA / UMAP / t-SNE customization -> `data-visualization/dimensionality-reduction-plots`
- Heatmap customization and ComplexHeatmap recipes -> `data-visualization/heatmaps-clustering`

## The Single Most Important Modern Insight -- A volcano with shrunken LFC compresses the cloud, but the p-values are unchanged

`lfcShrink()` pulls noisy estimates toward zero. On the volcano, that pulls genes horizontally toward the center. But the y-axis (`-log10(pvalue)`) is the **unshrunken Wald p-value** -- shrinkage does NOT recompute p-values (Zhu, Ibrahim, Love 2019 *Bioinformatics* 35:2084). A naive reader sees fewer extreme dots and concludes "fewer genes are significant". Wrong: the same genes are significant; the effect sizes are smaller and more honest.

Always label the volcano x-axis "shrunken log2 fold change (apeglm)" and note the y-axis comes from the unshrunken Wald test. The whole point of the apeglm volcano is the honest effect-size axis; if a publication shows an unshrunken volcano, it is showing inflated effects from low-count noise.

The MA plot has its own version of this: shrinkage flattens the left side (low-mean, formerly extreme LFC) and barely touches the right (high-mean, well-estimated LFC). That asymmetry is the visual signature of working shrinkage.

## Plot Taxonomy

| Plot | Diagnostic OR result | Built-in function | What it tests |
|------|---------------------|-------------------|---------------|
| Dispersion plot | Diagnostic | `plotDispEsts(dds)` (DESeq2), `plotBCV(y)` (edgeR) | Mean-dispersion trend fit quality |
| p-value histogram | Diagnostic | None; use ggplot2 | Null calibration, hidden batch, over-correction |
| PCA on VST/rlog | Diagnostic + result | `plotPCA(vsd, intgroup=...)` (DESeq2), `plotMDS()` (edgeR via limma) | Sample clustering, batch effects, outliers |
| Sample distance heatmap | Diagnostic | `pheatmap` on `dist(t(assay(vsd)))` | Within-group consistency, sample swaps |
| MA plot | Diagnostic + result | `plotMA(res)` (DESeq2), `plotMD(qlf)` (edgeR) | Normalization sanity, LFC vs mean |
| Volcano | Result | ggplot2 wrapper; `EnhancedVolcano` | Top-effect, top-significance gene story |
| Top-DE heatmap | Result | `pheatmap` on `assay(vsd)[sig_genes,]` | Per-gene pattern across conditions |
| `plotCounts` per gene | Result | `plotCounts(dds, gene, intgroup)` | Per-gene biology |

## Decision Tree by Scenario

| Scenario | Recommended approach |
|----------|---------------------|
| PCA for unbiased QC | `vst(dds, blind = TRUE)`; ask "do samples group as expected without design influence?" |
| PCA for results figure | `vst(dds, blind = FALSE)`; design is settled, accept its influence on dispersion |
| n < 30, library sizes vary >4x | `rlog(dds, blind = FALSE)` instead of vst |
| n > 30 | `vst()`; rlog impractical |
| Volcano | Plot shrunken LFC on x, unshrunken p-value on y; label both axes |
| Sample distance heatmap | `vst(blind = TRUE)`; tells if a sample is the wrong group regardless of design |
| Top-DE heatmap, want to see PATTERN | `scale = 'row'` (z-score per gene) |
| Top-DE heatmap, want to see ABSOLUTE LEVEL | `scale = 'none'` on `assay(vsd)`; otherwise weak signal looks strong |
| Top-variable-gene selection | `matrixStats::rowMads(assay(vsd))` instead of `rowVars` -- MAD is outlier-robust |
| n = 3, top genes in volcano | Note Schurch 2016 finding: 20-40% of true positives missed; treat as exploratory |
| Many groups, comparing DE sets | UpSet plot (Lex 2014); Venn drowns above 3 sets |

## Dispersion Diagnostic (Run This First)

**Goal:** Verify the dispersion-mean trend was fit acceptably before trusting any results.

**Approach:** `plotDispEsts(dds)` (DESeq2) or `plotBCV(y)` (edgeR) shows gene-wise (black/blue), fitted trend (red), and final shrunken (blue) dispersions vs mean.

```r
plotDispEsts(dds)

plotBCV(y)
```

| Pattern | Meaning | Action |
|---------|---------|--------|
| Cloud follows trend; final shrunken estimates pulled toward red curve | Healthy fit | Proceed |
| Red trend nowhere near the gene-wise cloud | Parametric trend failed | `DESeq(dds, fitType = 'local')` or `fitType = 'mean'` |
| Many gene-wise dispersions FAR ABOVE the trend | Outlier or unmodeled batch genes | Inspect rather than trust QL F-test alone |
| Final estimates much lower than gene-wise everywhere | Excessive shrinkage; sample too small or trend too flat | Check `useEM`, robust hyperparameter setting |
| BCV decreases monotonically with mean | Correct in edgeR | Default trend |

A plot inspected before trusting results is worth a hundred lines of statistical safeguards.

## P-value Histogram (Run This Second)

**Goal:** Detect model misspecification or hidden batch before reporting any gene list.

**Approach:** Histogram of raw p-values; under a correctly specified null, uniform with a spike near zero.

```r
library(ggplot2)
ggplot(res_df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = 'steelblue', color = 'white') +
    labs(x = 'P-value', y = 'Frequency', title = 'P-value distribution') +
    theme_bw()
```

| Shape | Meaning | Action |
|-------|---------|--------|
| Uniform + spike at 0 | Correctly specified | Proceed |
| U-shape (spikes at 0 AND 1) | Anti-conservative; hidden batch or unmodeled covariate | Add the missing covariate; re-fit |
| Depleted near 0, spike near 1 | Conservative; over-modeled or wrong dispersion | Simplify model; check dispersion plot |
| Spike only at p = 1 | Discrete artifact from very-low-count genes | Pre-filter more aggressively |

## MA Plot (LFC vs Mean)

**Goal:** Inspect the relationship between LFC and mean expression for normalization correctness and shrinkage effect.

**Approach:** `plotMA` (DESeq2) or `plotMD` (edgeR). Always pick `ylim` deliberately; default can flatten the signal.

```r
plotMA(res, ylim = c(-5, 5), main = 'MA plot (unshrunken)')

res_apeglm <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
plotMA(res_apeglm, ylim = c(-5, 5), main = 'MA plot (apeglm-shrunken)')

plotMD(qlf, main = 'edgeR MD plot')
abline(h = c(-1, 1), col = 'blue', lty = 2)
```

| Pattern | Meaning |
|---------|---------|
| Symmetric cloud centered at LFC = 0 | Correct normalization |
| Cloud median clearly above or below 0 | Normalization failed (TMM/RLE assumption violated) -- see normalization skill |
| Funnel widening at low mean | Expected (low counts noisier) |
| Dramatic up/down asymmetry | Possibly real (large biological perturbation), possibly normalization failure -- cross-check |
| Discrete horizontal bands at low mean | Low-count artifacts; pre-filter more aggressively |

The apeglm-shrunken MA visually flattens the left side; the post-shrinkage cloud should be tighter at low means.

## Volcano with Shrunken LFC

**Goal:** Show effect size vs significance with honest fold changes.

**Approach:** Use a built-in renderer (EnhancedVolcano for quick publication-quality output) on shrunken LFCs. Always plot shrunken LFC; always set `max.overlaps = Inf` when labeling >10 genes -- the ggrepel default (10) silently drops labels. EnhancedVolcano accepts `max.overlaps` directly in 1.12+; version 1.10-1.11 has the older `maxoverlapsConnectors` argument (default 15); for either, falling back to `options(ggrepel.max.overlaps = Inf)` at the top of the script also works. For full ggplot2 customization (color schemes, faceting, label-set engineering), see `data-visualization/volcano-and-ma-plots`.

```r
library(EnhancedVolcano)

res_apeglm <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

EnhancedVolcano(res_apeglm,
    lab = rownames(res_apeglm),
    x = 'log2FoldChange', y = 'pvalue',
    pCutoff = 0.05, FCcutoff = 1,
    title = 'Treatment vs Control',
    subtitle = 'Shrunken LFC (apeglm); unshrunken Wald p',
    max.overlaps = Inf)
```

## PCA on VST/rlog (Never on Raw Counts)

**Goal:** Show sample clustering by condition; detect batch effects, swaps, outliers.

**Approach:** Variance-stabilize first (VST or rlog), THEN PCA. Raw counts make PC1 = library size; log(counts+1) makes PC1 = mean expression. Neither carries biological signal until variance is stabilized.

```r
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c('condition', 'batch'))

pca_df <- plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, 'percentVar'))

library(ggplot2)
ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch)) +
    geom_point(size = 4) +
    xlab(paste0('PC1: ', percentVar[1], '% variance')) +
    ylab(paste0('PC2: ', percentVar[2], '% variance')) +
    theme_bw()

library(limma)
plotMDS(cpm(y, log = TRUE), col = as.numeric(group), pch = 16)
```

`blind=TRUE` (default for `vst()`) re-estimates dispersions ignoring the design -- appropriate for unbiased QC ("are samples consistent independent of design?"). `blind=FALSE` uses the fitted dispersions -- appropriate for downstream visualization where the design is settled. Modern DESeq2 vignette recommends `blind=FALSE` for any plot after the model is fit.

| PCA pattern | Interpretation | Action |
|-------------|----------------|--------|
| Clear separation by condition on PC1 or PC2 | Strong biological signal | Proceed |
| Separation by batch, not condition | Batch effect dominates | Include batch in design; DO NOT subtract before DE (see batch-correction Nygaard 2016) |
| One sample far from its group | Outlier or swap | Check library QC; sex check; somalier |
| Condition signal on PC3+, not PC1-PC2 | Subtle effect | May still find DE; review dispersion plot |
| Two distinct sample clusters not explained by metadata | Hidden covariate | Investigate processing date, lane, machine |

## Sample Distance Heatmap (for QC)

```r
library(pheatmap)
vsd <- vst(dds, blind = TRUE)
sd <- dist(t(assay(vsd)))
mat <- as.matrix(sd)
ann <- data.frame(condition = colData(dds)$condition,
                  row.names = colnames(dds))
pheatmap(mat, annotation_col = ann, annotation_row = ann,
         clustering_distance_rows = sd, clustering_distance_cols = sd,
         color = colorRampPalette(c('white', 'steelblue'))(100),
         main = 'Sample distance (vst blind)')
```

The diagonal should be dark; within-group samples should cluster. A within-group sample distant from its peers is a candidate for sample swap.

## Top-DE Heatmap and the Row-Scaling Trap

**Goal:** Show expression patterns of significant genes across samples for results figure.

**Approach:** Use `vst(blind=FALSE)`, select top genes (by adjusted p-value or MAD-robust variance), choose scaling deliberately.

```r
library(pheatmap)

sig <- rownames(subset(res, padj < 0.01))[1:50]
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig, ]

mat_scaled <- t(scale(t(mat)))

ann_col <- data.frame(condition = colData(dds)$condition,
                      batch     = colData(dds)$batch,
                      row.names = colnames(mat))

pheatmap(mat_scaled, annotation_col = ann_col,
         show_rownames = FALSE,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         color = colorRampPalette(c('blue', 'white', 'red'))(100),
         main = 'Top 50 DE genes (z-scored per gene)')
```

`scale='row'` (z-score per gene) is the conventional choice for "show me patterns". It DESTROYS absolute expression level information -- a gene at 5-7 with mean 6 looks identical to a gene at 10-1000. For pattern detection: correct. For QC heatmaps showing batch shifts: WRONG -- use `scale='none'` on `assay(vsd)`.

Top-variable-gene selection robustness:

```r
library(matrixStats)
vars_mad <- rowMads(assay(vsd))
top500 <- order(vars_mad, decreasing = TRUE)[1:500]
```

`rowMads` (median absolute deviation) is outlier-robust; `rowVars` is dominated by single-outlier-sample genes. For exploratory PCA of "top variable genes", MAD selection avoids artifacts.

## Per-gene Plot

```r
plotCounts(dds, gene = 'GENE_NAME', intgroup = 'condition')

d <- plotCounts(dds, gene = 'GENE_NAME', intgroup = c('condition','batch'),
                returnData = TRUE)
library(ggplot2)
ggplot(d, aes(x = condition, y = count, color = batch)) +
    geom_jitter(width = 0.1, size = 3) +
    scale_y_log10() +
    ggtitle('GENE_NAME') +
    theme_bw()
```

With n=3, the boxplot is misleading (3 points per box). Prefer `geom_jitter` over `geom_boxplot` at small n.

## UpSet for Multi-set Comparisons

For >3 DE gene sets (e.g., contrasts treated_drugA, treated_drugB, treated_drugC each vs control), Venn diagrams become unreadable. UpSet (Lex et al. 2014 *IEEE Trans Vis Comput Graph* 20:1983) scales:

```r
library(UpSetR)
upset(fromList(list(drugA = sig_drugA, drugB = sig_drugB, drugC = sig_drugC)))
```

## Per-Method Failure Modes

### Volcano with unshrunken LFC -- inflated story

**Trigger:** `ggplot(res_df, aes(x=log2FoldChange, ...))` without `lfcShrink()`; extreme dots at the corners are low-count genes.

**Mechanism:** Unshrunken MLE LFCs are dominated by very-low-count genes whose log ratios are noisy. The visual top-left and top-right corners look impressive but are artifacts.

**Symptom:** Top genes by abs(LFC) are obscure low-count genes; reviewer asks "why are these the top hits?"

**Fix:** `res_apeglm <- lfcShrink(dds, coef=..., type='apeglm')`; plot from `res_apeglm`. Label axis "shrunken log2 fold change (apeglm)".

### ggrepel `max.overlaps` silently drops labels

**Trigger:** `geom_text_repel(data = top30, aes(label = gene))`; only 10 labels render.

**Mechanism:** Default `max.overlaps = 10`; warning printed but easily missed in a knitr/Quarto render.

**Symptom:** Reviewer asks "where is gene X?"; it was in `top30` but did not render.

**Fix:** `geom_text_repel(..., max.overlaps = Inf)` or `options(ggrepel.max.overlaps = Inf)` at top of script.

### PCA shows batch, not condition

**Trigger:** `plotPCA(vsd, intgroup='batch')` cleanly separates batches; `intgroup='condition'` does not separate.

**Mechanism:** Batch variance exceeds condition variance.

**Symptom:** Treatment effect looks weak; DE p-values inflated if batch not in design.

**Fix:** Include batch in design (`design = ~ batch + condition`). DO NOT use `removeBatchEffect` then re-do DE on corrected counts (Nygaard 2016 cardinal sin -- see `batch-correction`). For VISUALIZATION only, `removeBatchEffect` is OK.

### Heatmap row-scaling hid a sample-level shift

**Trigger:** QC heatmap with `scale='row'` looks consistent within group; downstream PCA shows clear sample outlier.

**Mechanism:** z-score per gene removes per-sample additive shifts. A sample that's globally inflated 1.5x looks identical to peers after row scaling.

**Symptom:** "The heatmap looked fine but PCA shows a problem."

**Fix:** For QC heatmaps, use `scale = 'none'` on `assay(vsd)` directly. For result heatmaps after QC is clean, `scale = 'row'` is the appropriate choice for pattern emphasis.

### Top-N-by-rowVars dominated by single-outlier-sample genes

**Trigger:** "Top 500 variable genes" PCA shows a striped pattern, one or two samples driving the spread.

**Mechanism:** `rowVars` is squared-deviation; one outlier sample of one gene inflates that gene's "variance" massively.

**Symptom:** Top variable gene list includes many genes where N-1 samples are flat and one sample is extreme.

**Fix:** `matrixStats::rowMads()` for MAD-based selection; or `genefilter::rowQ()`.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `plotPCA` reports only 2 PCs | DESeq2 `plotPCA` is hard-coded to PC1/PC2 | Use `prcomp(t(assay(vsd)))` and plot any pair |
| PCA cloud collapses to one point | Forgot to log-transform; raw counts plotted | `vst(dds)` first |
| All MA-plot points red | `alpha` set too high or sig-flag bug | Verify `alpha`; check `padj` vs `pvalue` in flag |
| `pheatmap` complains "infinite values" | NA / Inf in scaled matrix; gene with zero variance | Remove zero-variance rows before scaling |
| Volcano axis labels obscured | Default ggplot theme too compact | `theme_bw(base_size = 14)` |
| `plotCounts` says gene not found | Wrong ID type (symbol vs Ensembl) | Match `rownames(dds)` exactly |
| `vst()` errors with very low gene count post-filter | Default `nsub=1000` exceeds available genes | Lower `nsub` (e.g., `vst(dds, nsub=500)`) |

## References

- Anders S, Huber W. 2010. Differential expression analysis for sequence count data. *Genome Biol* 11(10):R106. doi:10.1186/gb-2010-11-10-r106
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092. doi:10.1093/bioinformatics/bty895
- Robinson MD, McCarthy DJ, Smyth GK. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26(1):139-140. doi:10.1093/bioinformatics/btp616
- Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H. 2014. UpSet: Visualization of Intersecting Sets. *IEEE Trans Vis Comput Graph* 20(12):1983-1992. doi:10.1109/TVCG.2014.2346248
- Schurch NJ et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22(6):839-851. doi:10.1261/rna.053959.115
- Nygaard V, Rødland EA, Hovig E. 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17(1):29-39. doi:10.1093/biostatistics/kxv027

## Related Skills

- deseq2-basics - Generates the `dds` / `res` objects plotted here; `vst`/`rlog` choice
- edger-basics - Generates `y` / `qlf` for plotMD, plotBCV, plotMDS
- de-results - p-value histogram, padj=NA diagnosis informs what to plot
- batch-correction - removeBatchEffect for visualization only (never as DE input)
- expression-matrix/normalization - VST vs rlog vs log-CPM mechanics
- data-visualization/volcano-and-ma-plots - Full custom volcano/MA with apeglm + ggrepel
- data-visualization/dimensionality-reduction-plots - PCA, UMAP, t-SNE customization
- data-visualization/heatmaps-clustering - pheatmap and ComplexHeatmap recipes
- data-visualization/upset-plots - UpSet plot customization
