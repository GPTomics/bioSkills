---
name: bio-proteomics-differential-abundance
description: Tests for differentially abundant proteins between conditions with limma/DEqMS empirical-Bayes moderation, proDA/msqrob2/MSstats missingness modeling, and Python Welch+BH alternatives. Frames missing values as left-censored MNAR (model, do not impute), makes variance moderation the load-bearing step at n=3-5, and prefers feature/peptide-level testing. Use when identifying proteins with significant abundance changes between experimental groups. Summarization and normalization mechanics are proteomics/quantification; volcano and MA plots are data-visualization/volcano-and-ma-plots; pathway enrichment of the hit list is pathway-analysis/go-enrichment.
tool_type: mixed
primary_tool: limma
---

## Version Compatibility

Reference examples tested with: limma 3.58+, DEqMS 1.20+, proDA 1.20+, ashr 2.2+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Protein Abundance -- Moderated Testing on a Log-Intensity Matrix with Honest Missingness

**"Find differentially abundant proteins between my conditions"** -> Moderated statistical testing on a normalized log-intensity matrix, carrying missingness in the likelihood instead of filling it in -- because the missing values are low BECAUSE the protein is low, and imputing them manufactures false positives.
- R: `limma::eBayes(fit, trend=TRUE, robust=TRUE)` for empirical-Bayes moderated t-tests (the protein-level workhorse)
- R: `DEqMS::spectraCounteBayes()` when PSM/peptide counts are available (preferred over limma-trend when quant depth varies)
- R: `proDA::test_diff()` / `msqrob2` / `MSstats` when missing values are extensive (model the dropout, no imputation)
- Python: `scipy.stats.ttest_ind(equal_var=False)` + `statsmodels` BH (large n only; no moderation)

Scope: this skill owns the statistical TEST -- design/contrast construction, variance moderation, missingness handling, multiple-testing correction, minimum-fold-change testing, and fold-change shrinkage. Peptide-to-protein summarization and normalization mechanics -> proteomics/quantification. Volcano/MA plots -> data-visualization/volcano-and-ma-plots. Enrichment of the hit list -> pathway-analysis/go-enrichment. OUT OF SCOPE: how MaxLFQ/TMP/IRS produce the matrix (quantification); how to draw a volcano (data-visualization).

## The Single Most Important Modern Insight -- Model the Missingness, Moderate the Variance, Test at the Feature Level

1. **Missing values in label-free MS are left-censored MNAR -- missing BECAUSE the intensity is low -- and imputing them (especially Perseus/MaxQuant downshift) manufactures SYSTEMATIC false positives.** Downshift draws each missing value from a narrow Gaussian (mean = mu - 1.8*sigma, SD = 0.3*sigma). For an on/off protein (seen in all of group A, missing in all of B) the t-statistic numerator is inflated by construction (mean_B fixed ~1.8 sigma below observed, deterministic) and the denominator is artificially deflated (all imputed B values from one 0.3-sigma Gaussian -> collapsed within-group SD) -> enormous t -> tiny p. Because every on/off protein is treated identically, the false positives are systematic: the volcano-plot "anchor/wing" artifact (rigid near-vertical streaks of pinned points far out on both x-axis sides). The honest statement is "undetected in group B", not "20x lower, p=1e-6". The correct approach is to MODEL the dropout in the likelihood -- proDA (probabilistic dropout), msqrob2, MSstats-AFT -- NOT fill it (Lazar 2016; Ahlmann-Eltze & Anders 2019).
2. **At the n=3-5 replicates proteomics actually uses, per-protein variance has only 2-4 residual df and is unusable raw -- variance moderation is the load-bearing element, not optional.** limma borrows a prior d0 across all proteins so a 4-replicate design tests on ~10 df instead of 6; `trend=TRUE` makes the prior a function of mean intensity (effectively mandatory for label-free, where a single global prior mis-calibrates FDR across the abundance range); `robust=TRUE` Winsorizes outlier variances (Phipson 2016). DEqMS makes the prior a function of PSM/peptide count and generally outperforms limma-trend when quantification depth varies across proteins (Zhu 2020).
3. **Feature/peptide-level modeling beats summarize-then-test.** Summarizing first (one number per protein per run) discards the within-protein between-peptide variance and the correct degrees of freedom: 12 consistent peptides deserve a smaller SE than 12 disagreeing ones, but after summarization both look equally certain, and a protein with 30 observations looks as informative as one with 3. msqrob2/MSstats keep every peptide as a degree of freedom; this is why the same data gives different answers (Goeminne 2016; Sticker 2020; Choi 2014).

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---------------|----------|------------------|------|
| limma | Ritchie 2015; Phipson 2016 | EB moderated t; posterior variance blends a prior d0 with the per-protein estimate; `trend` ties the prior to mean intensity, `robust` Winsorizes outliers | protein-level summaries, small n, the default workhorse |
| DEqMS | Zhu 2020 | prior variance = loess of log-variance vs log2(count); precision follows quantification DEPTH not just intensity | TMT (count=PSM) and label-free DDA (count=peptide); quant depth varies; preferred over limma-trend |
| proDA | Ahlmann-Eltze & Anders 2019 (preprint) | probabilistic dropout: missing = left-censored, integrated under a per-sample sigmoid dropout curve; EB on location and variance; no imputation | label-free DDA with many MNAR missing values, small n, proteins absent in one group |
| msqrob2 | Sticker 2020; Goeminne 2016 | peptide-level robust ridge: Huber M-estimation downweights outlier peptides, ridge shrinks effects from few observations, EB variance moderation | label-free DDA, outlier-peptide / unbalanced-coverage risk; best FDR in hard spike-in regimes |
| MSstats | Choi 2014 | feature-level linear mixed model (group fixed + feature + run/subject random); AFT censored handling for missing | SRM/PRM/DIA, technical replicates, nested/repeated-measures, labeled designs |
| Welch t-test + BH | -- | per-protein two-sample t with `equal_var=False` + Benjamini-Hochberg | large n (>10/group), Python-only; no moderation, unusable at n=3-5 |
| ashr | Stephens 2017 | mixture prior with a point mass at zero; posterior means shrink uncertain effects toward zero | recovering "which proteins truly changed and by how much" (not for GSEA ranking) |
| volcano / MA plot | -- | (route OUT) | visualization -> data-visualization/volcano-and-ma-plots |
| enrichment of hits | -- | (route OUT) | functional interpretation -> pathway-analysis/go-enrichment |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Small n (3-5/group), protein-level summary matrix | limma `eBayes(trend=TRUE, robust=TRUE)` | EB borrows variance across proteins; the trend calibrates FDR across abundance |
| PSM/peptide counts available (TMT or label-free DDA) | DEqMS `spectraCounteBayes` | prior keyed on quant depth removes single-PSM false positives limma admits |
| Label-free with many MNAR missing values, on/off proteins | proDA `test_diff` | models the censored dropout; never imputes; correct verdict for "undetected in one group" |
| Outlier-peptide risk, unbalanced peptide coverage | msqrob2 (peptide-level robust ridge) | keeps feature df; Huber downweights bad peptides; best FDR in spike-in benchmarks |
| Technical replicates, nested/repeated-measures, labeled (SRM/PRM/DIA) | MSstats (feature-level mixed model) | random effects capture run/subject structure summarize-then-test discards |
| Batch present | batch as a covariate in the design (`~ batch + condition`) | `removeBatchEffect` is visualization-only; never feed its output to `lmFit` |
| Minimum biologically meaningful fold change | `treat()` + `topTreat()` (or SAM s0) | tests |log2FC|>c against the moderated null; a post-hoc FC+significance double filter inflates FDR |
| Large n (>10/group), Python-only | Welch t-test + BH | variance estimates reliable; no moderation needed at large n |

Default when uncertain: protein-level summary matrix at n=3-5 -> limma `eBayes(trend=TRUE, robust=TRUE)`; if PSM/peptide counts exist, escalate to DEqMS; if missingness is extensive and intensity-dependent, escalate to proDA.

## limma Workflow (R)

**Goal:** Identify differentially abundant proteins using moderated statistics that borrow information across all proteins.

**Approach:** Build the design (batch as a covariate when present), fit the linear model and contrast, apply EB moderation with the intensity trend and robust fitting, then extract BH-corrected results. Never feed `removeBatchEffect` output to `lmFit`.

```r
library(limma)

design <- model.matrix(~0 + condition + batch, data = sample_info)  # batch in the model, not removed first
colnames(design)[1:2] <- levels(factor(sample_info$condition))

fit <- lmFit(protein_matrix, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)  # trend mandatory for label-free; robust Winsorizes outliers

results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
# columns: logFC, AveExpr, t, P.Value, adj.P.Val, B  (adj.P.Val is the BH p; there is no $FDR)
```

### Minimum-fold-change testing

**Goal:** Call proteins whose effect exceeds a biologically meaningful threshold, not merely differ from zero.

**Approach:** Use `treat()` against the moderated null and read `topTreat()`. NEVER `topTable(lfc=...)` nor a post-hoc volcano double filter (`abs(logFC) > 1 & adj.P.Val < 0.05`); conditioning on both the FC and the p-value selects for high-variance nulls (a collider effect) and inflates realized FDR above 50% (Ebrahimpoor & Goeman 2021).

```r
LFC_THRESHOLD <- log2(1.2)  # 1.2-fold floor; treat tests against this null, no double-filter FDR inflation
fit2 <- treat(fit2, lfc = LFC_THRESHOLD)
results <- topTreat(fit2, coef = 1, number = Inf)  # topTreat omits the B column
```

## DEqMS Workflow (R)

**Goal:** Improve on limma by tying each protein's prior variance to its quantification depth -- proteins measured by more PSMs/peptides are more precise.

**Approach:** Run limma through `eBayes`, attach the count vector, then apply DEqMS's count-aware EB. Use PSM count for TMT (quant at MS2) and peptide count for label-free DDA; for multi-batch TMT use the MINIMUM count across batches (the bottleneck batch sets precision).

```r
library(DEqMS)

# fit2 is the limma fit through eBayes (above)
fit2$count <- psm_count_per_protein[rownames(fit2$coefficients)]  # PSM for TMT, peptide for LFQ; min across batches
fit3 <- spectraCounteBayes(fit2)

results <- outputResult(fit3, coef_col = 1)
# adds sca.t, sca.P.Value, sca.adj.pval (the count-adjusted statistics; use these, not the limma columns)
```

## proDA Workflow (R)

**Goal:** Test proteins with extensive MNAR missingness, including on/off proteins, without imputing a single value.

**Approach:** Fit the probabilistic-dropout model directly on the log-intensity matrix; missing values contribute as left-censored observations under a per-sample dropout curve. Then test the contrast against zero.

```r
library(proDA)

fit <- proDA(protein_matrix, design = ~condition, col_data = sample_info,
             reference_level = 'Control')
result_names(fit)  # list testable coefficients first
results <- test_diff(fit, conditionTreatment - conditionControl)
# columns: name, pval, adj_pval, diff (log2FC), t_statistic, se
```

## Python Workflow

**Goal:** Run the full pipeline in Python when no R is available and n is large enough that moderation is unnecessary.

**Approach:** Log2-transform, median-normalize, run per-protein Welch t-tests, apply Benjamini-Hochberg. This has NO variance moderation and should not be used at n=3-5 -- escalate to limma/DEqMS for small n.

```python
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

def preprocess(intensities):
    log2_data = np.log2(intensities.replace(0, np.nan))  # zeros -> NaN to avoid -inf
    sample_medians = log2_data.median(axis=0)
    return log2_data - sample_medians + sample_medians.median()

def differential_abundance(normalized, case_cols, ctrl_cols):
    rows = []
    for protein in normalized.index:
        case, ctrl = normalized.loc[protein, case_cols].dropna(), normalized.loc[protein, ctrl_cols].dropna()
        if len(case) >= 2 and len(ctrl) >= 2:
            _, pval = stats.ttest_ind(case, ctrl, equal_var=False)  # Welch; scipy defaults to Student's True
            rows.append({'protein': protein, 'log2fc': case.mean() - ctrl.mean(), 'pvalue': pval})
    df = pd.DataFrame(rows)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]  # default is Holm-Sidak; pass fdr_bh explicitly
    return df
```

## Fold-Change Reporting

**Goal:** Hand the right effect estimate to the right consumer.

**Approach:** Report the RAW fold change (the best unbiased point estimate) for GSEA/pathway ranking and meta-analysis -- those need the full continuous distribution or FC+SE pairs. Apply shrinkage (ashr) only when recovering "which proteins truly changed and by how much"; it fits a mixture prior with a point mass at zero and shrinks uncertain effects smoothly toward zero. This is preferred over hard-thresholding (zeroing FCs at padj 0.05), which creates an arbitrary step function. No mature Python ashr equivalent exists.

```r
library(ashr)

se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]
shrunk <- ash(fit2$coefficients[, 1], se, mixcompdist = 'normal')
shrunken_fc <- shrunk$result$PosteriorMean  # report alongside raw logFC, not as a replacement for GSEA
lfsr <- shrunk$result$lfsr
```

## Per-Method Failure Modes

### Downshift / any imputation feeding a variance-based test
**Trigger:** Perseus/MaxQuant downshift (or MinDet/MinProb/QRILC) fills NAs, then limma/t-test runs on the filled matrix.
**Mechanism:** Imputed values come from one narrow Gaussian -> fabricated low within-group variance + deterministic mean offset -> inflated t.
**Symptom:** Volcano "anchor/wing" -- rigid near-vertical streaks of pinned on/off proteins at high significance; realized FDR far above nominal.
**Fix:** Model the missingness instead (proDA / msqrob2 / MSstats-AFT); report on/off proteins as "undetected in group X".

### kNN imputation on left-censored data
**Trigger:** kNN/mean imputation applied to label-free data with MNAR dropout.
**Mechanism:** Mean-reverting -- pulls a truly-low (missing because low) value UP toward the mean.
**Symptom:** Real down-regulation is compressed; down hits weakened or lost.
**Fix:** Only valid under MCAR/MAR; for MNAR model the dropout. Under uncertainty Lazar 2016 shows the milder MCAR error beats MNAR-imputers slamming random highs to the floor.

### removeBatchEffect before testing
**Trigger:** `removeBatchEffect()` output fed to `lmFit`.
**Mechanism:** Subtracts the fitted batch component with no uncertainty propagation -> understated residual variance, inflated EB df; if batch is confounded with biology it deletes real signal.
**Symptom:** Anticonservative p-values; lost true effects when cases/controls split by batch.
**Fix:** Include batch as a covariate in the SAME model (`~ batch + condition`); use `removeBatchEffect` only for PCA/visualization.

### eBayes(trend=FALSE) on intensity data
**Trigger:** Plain `eBayes` (trend off) on a log-intensity matrix.
**Mechanism:** A single global prior over-shrinks high-abundance and under-shrinks low-abundance proteins.
**Symptom:** Mis-calibrated FDR across the abundance range.
**Fix:** `eBayes(trend = TRUE, robust = TRUE)`; escalate to DEqMS when quant depth varies.

### Wrong DEqMS count column
**Trigger:** Razor+unique counts vs MS2-level PSMs, or total-across-batches vs minimum-across-batches.
**Mechanism:** The variance-vs-count prior is fit on the wrong precision proxy.
**Symptom:** Mis-ranked proteins; the count moderation helps the wrong ones.
**Fix:** PSM count for TMT, peptide count for label-free; minimum count across batches for multi-batch TMT.

### proDA on MCAR missingness
**Trigger:** proDA applied where dropout is random (e.g. a TMT channel lost at random), not detection-limited.
**Mechanism:** The left-censored dropout model is mis-specified.
**Symptom:** Biased estimates; the model fits a dropout curve that does not exist.
**Fix:** proDA needs intensity-dependent missingness; for MCAR use limma/DEqMS on the observed values.

### FC + significance double filter
**Trigger:** `abs(logFC) > 1 & adj.P.Val < 0.05` applied after the test.
**Mechanism:** |logFC| is large for a true effect OR a large SE; filtering on both the FC and the p (both depend on SE) selects high-variance nulls (collider effect).
**Symptom:** Realized FDR above 50% at nominal 5% (Ebrahimpoor & Goeman 2021).
**Fix:** `treat()`+`topTreat()` or SAM s0, which sit inside the statistic before selection.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| n=3-5 replicates -> 2-4 residual df | -- | raw per-protein variance unusable; moderation is mandatory, not optional |
| limma adds prior d0 (~4) df | Ritchie 2015 | a 4-replicate design tests on ~10 df vs 6; the borrowed df is the benefit |
| downshift mean = mu - 1.8*sigma, SD = 0.3*sigma | Perseus default | 1.8 places imputed mass ~3.6th percentile; 0.3 gives only 30% of real spread -> manufactured false positives |
| `trend=TRUE` effectively mandatory for label-free | Ritchie 2015 | a single global prior mis-calibrates FDR across abundance |
| min-FC floor log2(1.2) (1.2-fold) via treat() | -- | example floor; common alternatives 1.5-fold (~0.58) or 2-fold (1.0); set by biology, tested against the moderated null |
| BH adjusted p < 0.05 | Benjamini-Hochberg | controls FDR over the WHOLE rejection set, not subsets carved out afterward |
| DEqMS multi-batch TMT: minimum count across batches | Zhu 2020 | the bottleneck batch sets the realized precision |
| realized FDR > 50% from FC+significance double filter | Ebrahimpoor & Goeman 2021 | top-100 at n=12 exceeded 50% FDR at nominal 5% |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `results$FDR` is NULL | limma `topTable`/`topTreat` have no `$FDR` column | use `adj.P.Val` (BH-adjusted p) |
| `topTreat` row has no `B` | `topTreat` omits `B` (a `topTable` column) | read `logFC, AveExpr, t, P.Value, adj.P.Val` |
| FDR mis-calibrated across abundance | `eBayes` with `trend=FALSE` on intensity data | `eBayes(fit, trend = TRUE, robust = TRUE)` |
| min-FC test inflates FDR | `topTable(lfc=...)` or post-hoc volcano double filter | `treat(fit, lfc=log2(1.2))` then `topTreat()` |
| anticonservative p after batch correction | `removeBatchEffect` output fed to `lmFit` | put batch in the design: `~ batch + condition` |
| DEqMS columns missing | forgot `fit$count` or read limma columns | set `fit$count`, run `spectraCounteBayes`, read `sca.adj.pval` from `outputResult` |
| Student's t instead of Welch | `scipy.stats.ttest_ind` defaults `equal_var=True` | pass `equal_var=False` |
| p-values look like Holm-Sidak | `statsmodels` `multipletests` defaults to `'hs'` | pass `method='fdr_bh'` |
| volcano "anchor/wing" streaks | downshift/imputation feeding the test | model dropout (proDA/msqrob2/MSstats-AFT); report on/off proteins as undetected |

## References

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43(7):e47.
- Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. 2016. Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Ann Appl Stat* 10(2):946-963.
- Zhu Y, Orre LM, Zhou Tran Y, et al. 2020. DEqMS: a method for accurate variance estimation in differential protein expression analysis. *Mol Cell Proteomics* 19(6):1047-1057.
- Ahlmann-Eltze C, Anders S. 2019. proDA: probabilistic dropout analysis for identifying differentially abundant proteins in label-free mass spectrometry. *bioRxiv* 661496 (preprint; cite `citation("proDA")`, never a journal).
- Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek O. 2014. MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments. *Bioinformatics* 30(17):2524-2526.
- Goeminne LJE, Gevaert K, Clement L. 2016. Peptide-level robust ridge regression improves estimation, sensitivity, and specificity in data-dependent quantitative label-free shotgun proteomics. *Mol Cell Proteomics* 15(2):657-668.
- Sticker A, Goeminne L, Martens L, Clement L. 2020. Robust summarization and inference in proteome-wide label-free quantification. *Mol Cell Proteomics* 19(7):1209-1219.
- Lazar C, Gatto L, Ferro M, Bruley C, Burger T. 2016. Accounting for the multiple natures of missing values in label-free quantitative proteomics data sets to compare imputation strategies. *J Proteome Res* 15(4):1116-1125.
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294.
- Ebrahimpoor M, Goeman JJ. 2021. Inflated false discovery rate due to volcano plots: problem and solutions. *Brief Bioinform* 22(5):bbab053.

## Related Skills

- quantification - peptide-to-protein summarization, normalization, and IRS that produce the matrix this skill tests
- proteomics-qc - quality control and batch-effect assessment before testing
- protein-inference - razor/shared-peptide ambiguity that drives which protein group gets the quantity
- ptm-analysis - site-level differential testing for modified peptides
- differential-expression/de-results - analogous empirical-Bayes interpretation for RNA-seq DE
- data-visualization/volcano-and-ma-plots - volcano and MA plots of the result table
- pathway-analysis/go-enrichment - functional enrichment of the significant protein hit list
- machine-learning/biomarker-discovery - building predictive panels from differential proteins
- workflows/proteomics-pipeline - end-to-end pipeline that calls this skill as the testing stage
