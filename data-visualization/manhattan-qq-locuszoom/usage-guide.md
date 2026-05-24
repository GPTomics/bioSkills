# Manhattan, QQ, and Locuszoom Plots - Usage Guide

## Overview

GWAS / TWAS / PWAS / QTL summary statistics need four canonical visualizations: Manhattan (per-variant -log10 p across the genome), QQ (expected vs observed p quantiles + λGC), Miami (two traits mirrored), and locuszoom-style regional plot (LD-colored zoom into a locus with gene track). The decisions that matter: which significance threshold matches the analysis, how to handle extreme-tail peaks, when λGC indicates confounding vs polygenicity, and which LD reference population matches the GWAS.

## Prerequisites

```r
install.packages(c('qqman', 'CMplot', 'locuszoomr', 'LDlinkR'))
BiocManager::install(c('EnsDb.Hsapiens.v86'))
```

```bash
pip install matplotlib pandas numpy scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Make a Manhattan plot from my GWAS summary statistics with significance at 5e-8"
- "Compute λGC and plot QQ with the value in the title"
- "Cap the Manhattan y-axis at 25 with arrows for capped points"
- "Mirror two traits as a Miami plot for shared-locus discovery"
- "Build a locuszoom regional plot around TCF7L2 with EUR LD coloring"

## Example Prompts

### Standard GWAS

> "Manhattan + QQ from a common-variant GWAS summary table. Significance at 5e-8 (Pe'er 2008). Cap y-axis at 25; mark capped points with arrows; report λGC in QQ title."

### TWAS

> "Manhattan for TWAS results from ~20000 genes. Bonferroni-correct threshold (2.5e-6) instead of 5e-8. Per-chromosome two-color alternation."

### Multi-ancestry meta-analysis

> "Plot a trans-ancestry GWAS with 5e-9 threshold. Match LD reference to the largest contributing population for any regional follow-up."

### Locuszoom

> "Zoom into a 1 Mb window around TCF7L2 with locuszoomr. Color SNPs by LD r² to the lead SNP using 1000G EUR. Overlay recombination rate and gene track."

### Inflation diagnostic

> "Compute λGC for these summary stats. If > 1.1, run LDSC and report whether intercept (confounding) or slope (polygenicity) is elevated."

## What the Agent Will Do

1. Load summary statistics, verify columns (CHR, BP, P, SNP) and sort by chromosome + position.
2. Compute λGC = median(chi-square)/0.4549 (median of chi-square_1, Devlin-Roeder 1999); compute sample-size adjusted λ_1000 if requested.
3. Render QQ plot with diagonal reference; annotate λGC in title.
4. Render Manhattan with two-color chromosome alternation; significance line at the analysis-appropriate threshold.
5. Cap extreme-tail peaks with marker indicators OR use a split axis.
6. Label lead SNPs (smallest p per locus, with a configurable distance window for "lead").
7. For locuszoom: match LD reference to ancestry; overlay gene track and recombination rate.
8. Export at 300 DPI; vector PDF for axes + raster for the scatter layer to keep file size sane.

## Tips

- **5e-8 is for common-variant European GWAS only** (Pe'er 2008). Use ~5e-9 for WGS EUR (Pulit 2017; Xu 2014), Bonferroni / n_genes for TWAS/PWAS, 5e-9 for trans-ancestry meta.

- **λGC alone does not distinguish confounding from polygenicity.** Run LDSC: intercept > 1 = confounding; slope > 0 = polygenicity (Bulik-Sullivan 2015).

- **λ > 1.05 - investigate. λ > 1.10 - almost certainly confounded.** Add 10-20 PCs or switch to LMM (BOLT-LMM, GEMMA, SAIGE).

- **Sample-size adjusted λ_1000 = 1 + (λ-1) × 1000/n.** Compare λ across studies of different N.

- **Always sort by CHR, BP before Manhattan.** qqman plots in input order; unsorted produces vertical scatter at wrong x.

- **Two-color chromosome alternation** (`c('#0072B2', '#56B4E9')`). 22 distinct hues add no information.

- **Cap extreme y with indicator.** Mark capped points with `^` arrow; annotate axis "(capped at 25)". Numerical value in caption.

- **Match LD reference to GWAS population.** EUR 1000G LD applied to Japanese GWAS gives wrong r² coloring in locuszoom.

- **Suggestive threshold 1e-5** is for follow-up prioritization, NOT for claims of significance.

- **QQ plot deflation below diagonal** usually indicates conservative test statistic or computational bug - investigate.

- **Lead SNP labeling**: pick the smallest p per LD block / per ±500 kb window. Avoid labeling 5 SNPs at the same locus.

## Related Skills

- population-genetics/association-testing - Run the GWAS
- workflows/gwas-pipeline - End-to-end QC + GWAS + plots
- causal-genomics/fine-mapping - Within-locus fine-mapping after locuszoom
- causal-genomics/colocalization-analysis - Two-trait colocalization
- phasing-imputation/imputation-qc - Pre-GWAS QC affects QQ plot
