# Differential CpG Testing - Usage Guide

## Overview

Per-CpG (per-site) differential methylation testing - finding individual CpGs (DMC/DMP) that differ between groups. The decisive first choice is the data object: bisulfite sequencing yields integer (methylated, total) counts whose coverage carries precision, so it routes to a beta-binomial / overdispersion-corrected count model (DSS, methylKit); arrays and any continuous beta matrix have no coverage, so they route to limma moderated-t on M-values. Throughout: test on the M-value scale, report effect on the beta scale (delta-beta), and gate hits on the intersection of FDR and effect size. Also covers scanning for differential VARIABILITY (variance, not mean) with DiffVar/iEVORA.

## Prerequisites

```bash
pip install numpy pandas scipy statsmodels
```

```r
BiocManager::install(c('limma', 'DSS', 'methylKit', 'missMethyl'))
```

Conceptual prerequisites and notes:
- Input from sequencing is per-CpG (methylated, total) counts (Bismark coverage / cytosine report -> methylation-calling). Input from arrays is a beta-value (or M-value) probe matrix.
- The minimum-coverage floor must hold in ALL compared samples; one low-coverage sample makes a CpG unusable for a group test.
- Cell-type composition is the dominant confounder in bulk-tissue EWAS - have cell-fraction estimates ready to add to the design matrix.
- Region-level questions (DMRs) are a separate, complementary step downstream of per-site testing.

## Quick Start

Tell your AI agent what you want to do:
- "Test each CpG for differential methylation between my case and control groups"
- "I have WGBS counts with replicates - use DSS beta-binomial per-site testing"
- "Run limma on M-values for my EPIC array, with cell-fraction covariates"
- "Recompute delta-beta from raw betas and gate hits on FDR < 0.05 and |delta-beta| > 0.2"
- "Scan for differentially variable CpGs, not just mean differences"

## Example Prompts

### Sequencing counts (count model)
> "I have per-CpG methylated and total read counts for 3 control and 3 treated WGBS samples. Run a beta-binomial per-site test with DSS, then keep CpGs with FDR < 0.05 and |delta-beta| >= 0.1."

> "Use methylKit to test my RRBS replicates per CpG with overdispersion correction and an F-test, then report hyper- and hypo-methylated sites at 25% difference and q < 0.01."

### Array / continuous matrix (limma)
> "I have an EPIC array beta matrix for 6 cases and 6 controls plus estimated cell fractions. Convert to M-values, fit limma with trend and robust empirical Bayes, include the cell fractions as covariates, and report adj.P.Val with delta-beta from the raw betas."

> "Run a quick continuous per-CpG Welch test on my high-coverage methylation matrix and apply BH FDR - I understand this loses coverage information."

### Effect size, multiple testing, variability
> "Filter my differential methylation results to the intersection of FDR < 0.05 and |delta-beta| >= 0.2, and flag the top discovery effect sizes as upward-biased."

> "Scan my methylation data for differentially variable CpGs with DiffVar on M-values, and co-report the mean delta-beta for each hit."

## What the Agent Will Do

1. Determine the data object (sequencing counts vs array/continuous beta) - this dictates the model.
2. For sequencing, assemble (chr, pos, M, Cov) per sample; for arrays, take the beta/M-value matrix.
3. Apply coverage filtering for sequencing (minimum 10x in every sample, 99.9th-percentile upper cap).
4. Route to the right test: DSS or methylKit (overdispersion-corrected) on counts; limma moderated-t on M-values for arrays/continuous; Welch quick-look only when explicitly continuous.
5. Build the design matrix, adding cell-fraction and other covariates for bulk-tissue EWAS.
6. Apply BH-FDR (or an EWAS genome-wide threshold for array consortium work).
7. Compute delta-beta from the raw beta values and gate hits on the intersection of FDR and |delta-beta|.
8. Optionally run a differential-variability scan (DiffVar/iEVORA on M-values) alongside the mean test.
9. Output a results table with per-CpG statistics, adjusted p-values, delta-beta, and significance calls.

## Tips

- The data object decides the model: counts -> a count model (DSS/methylKit) that uses coverage as precision; continuous beta -> a Gaussian model on M-values. A bare-beta t-test on sequencing counts discards coverage and is the central anti-pattern.
- methylKit's `calculateDiffMeth` defaults to `overdispersion="none"` (no correction, over-calls with replicates) and `adjust="SLIM"` (not BH). Set `overdispersion="MN", test="F", adjust="BH"`.
- Test on M-values, report effect on beta. The M-scale logFC (limma) and methylKit's `meth.diff` are NOT delta-beta - recompute delta-beta from raw betas.
- In Python, `scipy.stats.ttest_ind` defaults to Student's; pass `equal_var=False` for Welch. `statsmodels` `multipletests` defaults to Holm-Sidak; pass `method='fdr_bh'`.
- The limma adjusted-p column is `adj.P.Val`, not `padj` (DESeq2) or `FDR`.
- Fisher's exact is for unreplicated designs only (n=1 vs n=1); never pool biological replicates into one super-sample.
- Gate hits on BOTH FDR and |delta-beta|; at millions of CpGs a tiny delta-beta within noise can clear FDR.
- Discovery effect sizes are upward-biased (winner's curse) and attenuate on replication - treat them as upper bounds.
- Run differential-variability tests (DiffVar/iEVORA) on M-values, not beta, and co-report the mean to rule out a mean-at-boundary artifact.
- Per-site FDR is conservative when neighboring CpGs are correlated; for regional inference move to region-level methods.

## Related Skills

- methylation-calling - Produces the (M, coverage) counts tested here
- methylkit-analysis - methylKit object model and calculateDiffMeth mechanics
- dmr-detection - Region-level aggregation downstream of per-site testing
- cell-type-deconvolution - Cell-fraction covariates (the dominant bulk-tissue confounder)
- ewas-design - Covariate strategy, genomic inflation, and genome-wide thresholds
- experimental-design/multiple-testing - FDR/FWER theory behind the corrections applied here
- long-read-sequencing/nanopore-methylation - Long-read MM/ML calling; pipe per-site counts here for count-based statistics
- differential-expression/deseq2-basics - Analogous dispersion-shrinkage / empirical-Bayes machinery
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
