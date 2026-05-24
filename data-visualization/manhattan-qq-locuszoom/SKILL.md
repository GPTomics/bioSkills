---
name: bio-data-visualization-manhattan-qq-locuszoom
description: Build Manhattan, Miami, QQ, and locuszoom-style regional plots from GWAS, TWAS, PWAS, and QTL summary statistics with correct genomic-inflation diagnostics, multi-trait overlays, lead-SNP labeling, and LD-aware regional rendering. Use when visualizing association results across the genome, comparing two traits, computing genomic inflation lambda, or zooming into a locus with LD coloring.
tool_type: mixed
primary_tool: qqman
---

## Version Compatibility

Reference examples tested with: qqman 0.1.9 (R), CMplot 4.5+ (R), matplotlib 3.8+, pandas 2.2+, scipy 1.12+, plinkQC 0.3+. For locuszoom-style: locuszoomr 0.3+ (R) or pyranges + matplotlib.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Manhattan, QQ, and Locuszoom Plots

**"Plot my GWAS results"** -> Render per-variant -log10(p) across the genome (Manhattan), compare expected vs observed p quantiles (QQ + λGC), overlay two traits with mirrored axes (Miami), and zoom into a locus with LD-colored points + recombination rate + gene track (locuszoom). The choices that matter: significance thresholds, axis truncation for ultra-significant peaks, lead-SNP labeling, and LD reference selection for regional plots.

- R: `qqman::manhattan` / `qqman::qq` (Turner 2018), `CMplot::CMplot`, `locuszoomr::locus_plot`
- Python: `matplotlib` + `pandas` for custom; `assocplots` for ready-made

## The Single Most Important Modern Insight -- The Threshold Is Always Conditional

The "genome-wide significant" line at `p < 5e-8` (Pe'er 2008 *Genet Epidemiol* 32:381) is calibrated for **European-ancestry common-variant GWAS** assuming ~1M effectively independent tests. It is the wrong threshold for:

- **Whole-genome sequencing** including rare variants (~5e-9 EUR, ~1e-9 AFR; Pulit 2017 *Genet Epidemiol* 41:145; Xu 2014 *Genet Epidemiol* 38:281)
- **Non-European ancestry** with different LD structure (typically more stringent)
- **TWAS / PWAS** with ~20,000 tested genes (Bonferroni 2.5e-6)
- **Multi-ethnic meta-analysis** (5e-9 by convention for trans-ancestry)
- **Burden / SKAT rare-variant tests** (per-gene; ~2.5e-6)
- **Locus-wise fine-mapping** (within-locus testing, no genome-wide correction needed)

A Manhattan plot's significance line is a contract with the reader about which multiple-testing regime applies. Mismatched thresholds over- or under-report hits.

## Decision Tree by Analysis

| Analysis | Genome-wide threshold | Suggestive threshold | Reference |
|----------|----------------------|---------------------|-----------|
| Common-variant GWAS (Eur) | 5e-8 | 1e-5 | Pe'er 2008 *Genet Epidemiol* 32:381 |
| Whole-genome sequencing (all variants, EUR) | 5e-9 | 5e-8 | Pulit 2017 *Genet Epidemiol* 41:145; Xu 2014 *Genet Epidemiol* 38:281 |
| Non-European ancestry (empirical per pop) | ~3.24e-8 AFR; ~9.26e-8 EAS | – | Kanai 2016 *J Hum Genet* 61:861 |
| TWAS (~20k genes) | 2.5e-6 (Bonferroni) | 1e-4 | Standard practice |
| PWAS (~5k proteins) | 1e-5 | 1e-4 | Standard practice |
| eQTL trans (genome-wide per probe) | Bonferroni over genes × variants | Per-tissue | GTEx convention |
| eQTL cis (within 1Mb) | nominal p < 1e-5 with permutation | – | GTEx FastQTL |
| Rare-variant gene burden | 2.5e-6 | 1e-4 | Bonferroni 20k genes |
| Trans-ancestry meta-analysis | 5e-9 | – | Convention |

## Genomic Inflation (λGC) -- The Mandatory QC Step

**Goal:** Quantify whether observed p-values are inflated relative to the chi-square null, indicating cryptic population structure, relatedness, or technical artifacts.

**Approach:** Convert observed p to chi-square; compute median chi-square divided by 0.4549 (the median of chi-square_1; Devlin-Roeder 1999); plot expected vs observed quantiles (QQ plot).

```r
library(qqman)
chisq <- qchisq(1 - df$P, df = 1)
lambda <- median(chisq) / 0.4549
# lambda = 1.0 -> no inflation
# lambda > 1.1 -> investigate; could indicate confounding
# lambda > 1.2 -> almost certainly confounded; principal components or LMM needed

qq(df$P, main = paste('QQ plot (lambda =', round(lambda, 3), ')'))
```

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

chisq = stats.chi2.isf(df['P'], df=1)
lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)

# QQ plot
expected = -np.log10(np.arange(1, len(df)+1) / (len(df) + 1))
observed = -np.log10(np.sort(df['P']))
fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(expected, observed, s=2)
ax.plot([0, max(expected)], [0, max(expected)], 'r--')
ax.set_xlabel(r'Expected $-\log_{10}(p)$')
ax.set_ylabel(r'Observed $-\log_{10}(p)$')
ax.set_title(f'QQ ($\\lambda_{{GC}} = {lambda_gc:.3f}$)')
```

**Interpretation**:
- λ = 1.00 ± 0.02 — well-calibrated
- λ > 1.05 — possible inflation; consider sample-size adjustment (`λ_1000 = 1 + (λ - 1) * 1000/n`)
- λ > 1.10 — confounded; population structure not removed; rerun with PC adjustment or LMM (BOLT-LMM, GEMMA, SAIGE)
- λ < 1.00 — deflation; usually a bug (wrong test statistic, conservative p-values)

Inflation can also be **legitimate polygenic signal** (Yang 2011 *Eur J Hum Genet* 19:807). Distinguish via LD-score regression intercept: confounding inflates intercept; polygenic signal inflates slope.

## Small-N and Rare-Variant Regimes -- When Standard Asymptotics Break

Standard logistic regression / Wald test p-values are anti-conservative when (a) case count < 200, (b) per-variant minor-allele count < 20, (c) case-control ratio is severely unbalanced (typical in EHR-derived cohorts). λGC may look normal but per-variant p-values are inflated independently — a Manhattan plot of these p-values is misleading regardless of inflation diagnostics.

| Regime | Test choice | Tool |
|--------|-------------|------|
| Balanced case-control, N>5000, MAF>0.01 | Standard logistic / linear regression | PLINK, REGENIE |
| Unbalanced (case fraction <10%), large N | SPA-corrected logistic regression | SAIGE, REGENIE Firth/SPA |
| Small N (<5000) | Penalized regression with bias correction | SAIGE Firth, REGENIE |
| Rare variants (MAC <20) | Gene-burden or SKAT-O | STAAR, REGENIE burden |

Specifically: **SAIGE** (Zhou 2018 *Nat Genet* 50:1335) and **REGENIE** (Mbatchou 2021 *Nat Genet* 53:1097) implement saddlepoint-approximation (SPA) and Firth-bias correction. Use them for any cohort with severe case-control imbalance; the Manhattan / QQ output is then defensibly calibrated.

## Manhattan Plot -- Canonical Layout

```r
library(qqman)
manhattan(df,
          chr = 'CHR', bp = 'BP', p = 'P', snp = 'SNP',
          col = c('#0072B2', '#56B4E9'),
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          ylim = c(0, max(-log10(df$P)) * 1.1),
          highlight = lead_snps,
          annotatePval = 5e-8,
          annotateTop = TRUE)
```

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(df, p_threshold=5e-8, suggestive=1e-5, y_cap=None):
    df = df.sort_values(['CHR', 'BP']).copy()
    df['neg_log10_p'] = -np.log10(df['P'])
    if y_cap:
        df['neg_log10_p'] = df['neg_log10_p'].clip(upper=y_cap)

    df['x'] = np.arange(len(df))
    chr_ticks = df.groupby('CHR')['x'].median()
    chr_colors = ['#0072B2', '#56B4E9']

    fig, ax = plt.subplots(figsize=(10, 4))
    for i, (chrom, group) in enumerate(df.groupby('CHR')):
        ax.scatter(group['x'], group['neg_log10_p'],
                   c=chr_colors[i % 2], s=3, rasterized=True)
    ax.axhline(-np.log10(p_threshold), color='red', linestyle='--', lw=0.5)
    ax.axhline(-np.log10(suggestive), color='grey', linestyle='--', lw=0.5)
    ax.set_xticks(chr_ticks)
    ax.set_xticklabels(chr_ticks.index, rotation=0)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(r'$-\log_{10}(p)$')
    return fig
```

## Extreme Tail Handling -- The Y-Axis Cap

Genome-wide significant peaks routinely reach -log10(p) = 100+ (e.g., GWAS of BMI at FTO). The visual effect: one peak fills the y-axis, all other signal is invisible.

**Fixes (ordered by preference):**

1. **Cap and indicate:** `y_cap = 25`; clip points above the cap; mark capped points with `^` arrow at the top of the panel. Use Y-axis label "−log10(P), capped at 25"
2. **Split y-axis** via `ggbreak::scale_y_break()` (R) or `axes_grid1.divider` (matplotlib)
3. **Two-panel plot** with full y range in top, zoomed range in bottom
4. **Use sqrt or asinh transform:** less intuitive but preserves all data; rarely chosen for Manhattan

## Miami Plot -- Two-Trait Comparison

```r
# CMplot supports Miami natively
library(CMplot)
CMplot(list(trait1_df, trait2_df),
       plot.type = 'm',
       multraits = TRUE,
       threshold = 5e-8,
       threshold.col = 'red',
       col = list(c('#0072B2','#56B4E9'), c('#D55E00','#E69F00')),
       file = 'jpg', file.output = TRUE)
```

Miami plot mirrors trait 1 above the x-axis, trait 2 below. Useful for shared-locus discovery (mirrored peaks at the same locus = pleiotropy candidate).

## Locuszoom-Style Regional Plot

A locuszoom plot zooms into a ~1Mb window around a lead SNP, colors SNPs by LD r² to the lead, overlays recombination rate (cM/Mb), and shows the gene track. The canonical tool is locuszoom.org (Pruim 2010 *Bioinformatics* 26:2336); the R package is `locuszoomr` (Lai 2024).

```r
library(locuszoomr)
loc <- locus(gene = 'TCF7L2',
             flank = 5e5,
             ens_db = 'EnsDb.Hsapiens.v86',
             data = gwas_df,
             snp = 'SNP', chrom = 'CHR', pos = 'BP', p = 'P', labs = 'SNP')
# LD computed via LDlinkR / 1000G reference
loc <- link_LD(loc, pop = 'EUR', token = ldlink_token)
locus_plot(loc, labels = c('index', 'top'))
```

**LD reference choice:** ALWAYS match the GWAS population. Using a 1000G European LD reference for a Japanese GWAS produces wrong LD colorings and misleads fine-mapping.

## Per-Method Failure Modes

### Inflation diagnosed as polygenic signal

**Trigger:** λGC = 1.15; analyst concludes "polygenic," moves on.

**Mechanism:** Inflation can be confounding (population structure, relatedness, technical) OR polygenicity. LD-score regression separates them: intercept = confounding; slope = polygenicity.

**Symptom:** Top hits replicate poorly in independent cohorts.

**Fix:** Run LDSC (`ldsc.py --h2`); intercept significantly > 1 indicates confounding. Adjust with 10-20 PCs or switch to LMM.

### Wrong significance threshold for the analysis

**Trigger:** Plotting Bonferroni-naive 5e-8 line on a TWAS or rare-variant burden plot.

**Mechanism:** 5e-8 is calibrated for ~1M independent common-variant tests; other analyses have different effective test counts.

**Symptom:** Reviewer flags "this gene doesn't pass Bonferroni" but the plot's red line is at 5e-8.

**Fix:** Match the threshold to the analysis (table above).

### Cap with no indication

**Trigger:** Y-axis clipped at 25 without arrow markers for capped points.

**Mechanism:** Reader cannot tell how high the true peak is.

**Symptom:** Reviewer asks for actual p-value; the reported 1e-100 conflicts with the displayed cap at 25.

**Fix:** Mark capped points with `^` symbol; annotate axis "(capped at 25)"; include unclipped numerical value in caption.

### Manhattan with random chromosome colors

**Trigger:** `col = rainbow(22)` for chromosome alternating.

**Mechanism:** 22 random hues add no information; visual chaos.

**Symptom:** Reader cannot quickly identify which chromosome a peak is on.

**Fix:** Two-color alternation (`c('#0072B2', '#56B4E9')`). Chromosome boundaries are clear from spacing alone.

### LD reference mismatched to GWAS population

**Trigger:** 1000G EUR LD reference used to color a Japanese / African GWAS regional plot.

**Mechanism:** LD differs by population; r² is population-specific.

**Symptom:** Locuszoom shows uncorrelated SNPs in red (high r²) or vice versa.

**Fix:** Match LD reference to GWAS population; for trans-ancestry GWAS, show per-population panels or use largest-N population reference and annotate the discrepancy.

### qqman::manhattan ignores BP order within chromosome

**Trigger:** Unsorted input data frame.

**Mechanism:** qqman plots in input row order, not coordinate order.

**Symptom:** Peaks render as vertical scatter at wrong x-position.

**Fix:** `df <- df %>% arrange(CHR, BP)` before plotting.

## Reconciliation: When QC Metrics Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| λGC > 1.1 but LDSC intercept ~1 | Polygenic signal | Document; no action needed |
| λGC > 1.1 AND LDSC intercept > 1 | Confounding | Add PCs / use LMM |
| QQ plot "S-shaped" | Severe inflation or non-additive model misspecification | Inspect; possibly model misspecified |
| QQ plot deflated below diagonal | Conservative p (e.g., score test); or wrong test stat | Review test statistic computation |
| Top SNP genome-wide but small effect | Likely true; or relatedness | Verify in unrelated subset |
| Replication fails for top hits | Confounding (winner's curse); or true heterogeneity | Trans-ancestry meta-analysis or LMM rerun |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Common-variant GWAS sig | 5e-8 | Pe'er 2008 |
| WGS sig (all variants, EUR) | 5e-9 | Pulit 2017; Xu 2014 |
| Empirical pop-specific (e.g., EAS) | ~9.26e-8 EAS | Kanai 2016 |
| TWAS / PWAS Bonferroni | 0.05 / n_genes | Standard |
| λGC well-calibrated | 1.00 ± 0.02 | Standard |
| λGC investigate | >1.05 | Common practice |
| λGC confounded | >1.10 | Common practice |
| Sample-size adjusted λ | λ_1000 = 1 + (λ - 1) × 1000/n | Standard scaling |
| Suggestive threshold | 1e-5 | Pe'er 2008 |
| Y-cap typical | 25-50 -log10(p) | Visualization choice |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Peaks at wrong x position | Data not sorted by CHR, BP | `arrange(CHR, BP)` upstream |
| λGC reported as conclusion alone | Confounding vs polygenicity not separated | Run LDSC for intercept vs slope |
| Threshold line at 5e-8 on TWAS | Wrong multiple-testing regime | Use Bonferroni-correct threshold |
| Y-axis crushed by one peak | No cap, no split | Cap at 25-50 with arrow markers OR split axis |
| Regional plot LD colors look wrong | LD reference mismatched to GWAS pop | Match LD reference to ancestry |
| QQ plot deflated | Conservative test or wrong stat | Verify test statistic |
| Manhattan with 22 distinct hues | Cosmetic clutter | Two-color alternation |
| Lead SNPs unlabeled | Default labeling off | `annotatePval = 5e-8, annotateTop = TRUE` |

## References

- Kanai M, Tanaka T, Okada Y. 2016. Empirical estimation of genome-wide significance thresholds based on the 1000 Genomes Project data set. *J Hum Genet* 61:861-866.
- Lai R. 2024. locuszoomr: an R/Bioconductor package for locus visualization. (CRAN package documentation)
- Pulit SL, de With SAJ, de Bakker PIW. 2017. Resetting the bar: statistical significance in whole-genome sequencing-based association studies of global populations. *Genet Epidemiol* 41(2):145-151.
- Xu C, Tachmazidou I, Walter K, et al. 2014. Estimating genome-wide significance for whole-genome sequencing studies. *Genet Epidemiol* 38(4):281-290.
- Pe'er I, Yelensky R, Altshuler D, Daly MJ. 2008. Estimation of the multiple testing burden for genomewide association studies of nearly all common variants. *Genet Epidemiol* 32:381-385.
- Pruim RJ, Welch RP, Sanna S, et al. 2010. LocusZoom: regional visualization of genome-wide association scan results. *Bioinformatics* 26:2336-2337.
- Pearson TA, Manolio TA. 2008. How to interpret a genome-wide association study. *JAMA* 299(11):1335-1344.
- Turner SD. 2018. qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. *J Open Source Softw* 3(25):731.
- Yang J, Weedon MN, Purcell S, et al. 2011. Genomic inflation factors under polygenic inheritance. *Eur J Hum Genet* 19(7):807-812.
- Bulik-Sullivan BK, Loh PR, Finucane HK, et al. 2015. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nat Genet* 47:291-295.
- Devlin B, Roeder K. 1999. Genomic control for association studies. *Biometrics* 55(4):997-1004.
- Mbatchou J, Barnard L, Backman J, et al. 2021. Computationally efficient whole-genome regression for quantitative and binary traits. *Nat Genet* 53(7):1097-1103.
- Zhou W, Nielsen JB, Fritsche LG, et al. 2018. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. *Nat Genet* 50(9):1335-1341.

## Related Skills

- population-genetics/association-testing - Run the GWAS that produces the summary stats
- workflows/gwas-pipeline - End-to-end GWAS workflow including QC
- causal-genomics/fine-mapping - Within-locus fine-mapping post-locuszoom
- causal-genomics/colocalization-analysis - Two-trait shared-causal-variant analysis
- data-visualization/color-palettes - Two-color chromosome alternation
- phasing-imputation/imputation-qc - Pre-GWAS imputation QC affects QQ
