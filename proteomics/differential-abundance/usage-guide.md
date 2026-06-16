# Differential Abundance - Usage Guide

## Overview
Identify proteins with significantly different abundance between experimental conditions. The central decisions are how to handle missing values (model the left-censored MNAR dropout, do not impute it), how to moderate per-protein variance at the small sample sizes proteomics uses, and at what level to test (protein summary vs feature/peptide). The skill covers limma, DEqMS, proDA, msqrob2, MSstats, and a Python Welch+BH fallback, plus minimum-fold-change testing and fold-change shrinkage.

## Prerequisites
```bash
pip install numpy pandas scipy statsmodels
```
```r
BiocManager::install(c("limma", "DEqMS", "proDA", "msqrob2", "MSstats", "ashr"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differentially abundant proteins between treatment and control in my intensity matrix"
- "Run limma with empirical-Bayes moderation on my small-n protein data"
- "Use DEqMS because I have PSM counts per protein from a TMT experiment"
- "My label-free data has 30% missing values and some on/off proteins -- test it without imputing"
- "Test for at least a 1.5-fold change instead of just nonzero, without inflating FDR"

## Example Prompts

### Choosing a Method
> "I have 4 replicates per group of protein-level LFQ intensities. Pick and run the right differential test."

> "Analyze my TMT proteomics data for differential abundance. I have PSM counts per protein, so use DEqMS, and remember it is multi-batch."

> "My label-free data has many proteins detected in one group but missing in the other. Test these honestly instead of imputing."

### Missingness and Batch
> "Set up a limma model with condition and batch as covariates -- do not remove the batch effect before testing."

> "Some proteins are missing because they are below the detection limit. Use a method that models this dropout."

### Effect Size and Thresholds
> "Test whether fold changes exceed 1.5-fold using treat and topTreat, not a post-hoc filter."

> "Report raw fold changes for GSEA and shrunken estimates for the figure."

## What the Agent Will Do
1. Load the protein (or peptide) intensity matrix and sample metadata, and assess missingness structure (how much, and whether it is intensity-dependent).
2. Confirm the matrix is log2-transformed and normalized (summarization/normalization mechanics live in proteomics/quantification).
3. Select the test: limma for small-n protein summaries, DEqMS when PSM/peptide counts exist, proDA/msqrob2/MSstats when missingness is extensive or feature-level structure matters, Welch+BH only at large n.
4. Build the design matrix with batch as a covariate (never `removeBatchEffect` before testing) and define contrasts.
5. Fit with empirical-Bayes variance moderation (`trend=TRUE, robust=TRUE`; DEqMS count prior when available).
6. Apply Benjamini-Hochberg correction over the whole rejection set.
7. For a minimum effect size, use `treat()`+`topTreat()` rather than a FC+significance double filter.
8. Choose fold-change reporting: raw for GSEA/meta-analysis, ashr-shrunk for effect-size recovery.
9. Produce the results table and hand off to volcano/heatmap and enrichment skills.

## Statistical Method Selection

| Method | Best for | Key advantage |
|--------|----------|---------------|
| limma | Small n (3-5), protein-level summaries | Borrows variance across proteins via empirical Bayes |
| DEqMS | PSM/peptide counts available | Prior keyed on quantification depth; dominates limma-trend |
| proDA | Label-free with extensive MNAR missing values | Models dropout in the likelihood; no imputation |
| msqrob2 | Outlier-peptide / unbalanced coverage | Peptide-level robust ridge; keeps feature df |
| MSstats | Technical replicates, nested/labeled designs | Feature-level mixed models |
| Welch t-test + BH | Large n (>10), Python-only | Simple; no moderation, unsuitable at small n |

## Missing-Value Handling

The dominant statistical problem is missingness, and in label-free MS it is left-censored MNAR -- missing because the intensity is low. Imputing it injects directional bias:
- Perseus/MaxQuant downshift manufactures systematic false positives (the volcano "anchor/wing" artifact) by giving on/off proteins an inflated mean offset and a collapsed within-group variance.
- kNN/mean imputation is mean-reverting and compresses real down-regulation.
- The correct approach models the dropout: proDA (probabilistic dropout), msqrob2, or MSstats with AFT censoring. Report on/off proteins as "undetected in group X", not as a giant fold change.

## Fold-Change Reporting

- GSEA / pathway analysis: use raw fold changes for all proteins; these methods rank by the full continuous distribution. Do not threshold before GSEA.
- Effect-size recovery ("which proteins truly changed and by how much?"): apply ashr shrinkage in R for posterior-mean estimates.
- Reporting tables: report raw FC with adjusted p-value and confidence interval.
- Meta-analysis: use raw FCs with standard errors; shrinkage is study-specific and applied after pooling, not before.

## Tips
- Confirm intensities are log2-transformed and normalized before testing; summarization/normalization mechanics belong to proteomics/quantification.
- Use `eBayes(trend = TRUE, robust = TRUE)` -- the trend is effectively mandatory for label-free intensity data.
- Use `equal_var=False` in `scipy.stats.ttest_ind` (the default is Student's), and pass `method='fdr_bh'` to `multipletests` (the default is Holm-Sidak).
- limma `topTable`/`topTreat` return `adj.P.Val` (there is no `$FDR` column); `topTreat` omits `B`.
- For a minimum effect size use `treat()`+`topTreat()`; never `topTable(lfc=...)` nor a post-hoc `abs(logFC)>1 & adj.P.Val<0.05` double filter (it inflates realized FDR above 50%).
- Include batch as a covariate in the design; use `removeBatchEffect()` only for PCA/visualization, never as input to `lmFit`.
- For DEqMS use PSM count for TMT and peptide count for label-free, and the minimum count across batches for multi-batch TMT.
- Check volcano symmetry: rigid near-vertical streaks of pinned points signal imputation artifacts, not biology.

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
