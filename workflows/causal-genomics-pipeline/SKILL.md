---
name: bio-workflows-causal-genomics-pipeline
description: End-to-end causal inference pipeline from GWAS summary statistics through Mendelian randomization, colocalization, fine-mapping, and mediation analysis. Use when performing post-GWAS causal inference to identify causal exposures, shared causal variants, and mediating mechanisms.
tool_type: r
primary_tool: TwoSampleMR
workflow: true
depends_on:
  - causal-genomics/mendelian-randomization
  - causal-genomics/colocalization-analysis
  - causal-genomics/fine-mapping
  - causal-genomics/pleiotropy-detection
  - causal-genomics/mediation-analysis
qc_checkpoints:
  - after_instrument_selection: "F-statistic > 10 for all instruments"
  - after_mr: "Consistent direction across MR methods"
  - after_sensitivity: "MR-PRESSO global test p > 0.05, Egger intercept p > 0.05"
  - after_coloc: "PP.H4 > 0.8 for strong colocalization"
---

# Causal Genomics Pipeline

Complete post-GWAS causal inference workflow from summary statistics to causal mechanisms.

## Pipeline Overview

```
GWAS Summary Statistics (exposure + outcome)
    |
    v
[1. Instrument Selection] ----> LD clumping, F-stat filtering
    |
    v
[2. Mendelian Randomization] -> IVW, MR-Egger, Weighted Median
    |
    +--> [3. Sensitivity Analysis] ---> MR-PRESSO, Egger intercept,
    |                                    leave-one-out, Steiger
    |
    v                           (can run in parallel with MR)
[4. Colocalization] ----------> coloc, eCAVIAR, HyPrColoc
    |
    v
[5. Fine-Mapping] ------------> SuSiE, FINEMAP (credible sets)
    |
    v
[6. Mediation Analysis] ------> Network MR, MVMR
    |
    v
Causal evidence summary (triangulation across methods)
```

## Step 1: Instrument Selection

```r
library(TwoSampleMR)

exposure_dat <- read_exposure_data(
    filename = 'exposure_gwas.tsv',
    sep = '\t',
    snp_col = 'SNP',
    beta_col = 'BETA',
    se_col = 'SE',
    effect_allele_col = 'A1',
    other_allele_col = 'A2',
    eaf_col = 'EAF',
    pval_col = 'P'
)

# Genome-wide significance threshold for instrument selection
# p < 5e-8: Standard GWAS threshold; use 5e-6 for underpowered exposures
exposure_dat <- subset(exposure_dat, pval.exposure < 5e-8)

# LD clumping to ensure independent instruments
# clump_r2 = 0.001: Strict independence; 0.01 for more instruments
# clump_kb = 10000: 10 Mb window, standard for MR
exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.001, clump_kb = 10000)
```

### QC Checkpoint: Instrument Strength

```r
check_instrument_strength <- function(exposure_dat) {
    # F-statistic > 10 to avoid weak instrument bias
    # F = beta^2 / se^2 for single-SNP instruments
    exposure_dat$F_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
    weak <- sum(exposure_dat$F_stat < 10)

    cat(sprintf('Instruments: %d SNPs\n', nrow(exposure_dat)))
    cat(sprintf('Mean F-statistic: %.1f\n', mean(exposure_dat$F_stat)))
    cat(sprintf('Weak instruments (F < 10): %d\n', weak))

    if (weak > 0) {
        cat('WARNING: Remove weak instruments before proceeding\n')
        exposure_dat <- subset(exposure_dat, F_stat >= 10)
    }

    if (nrow(exposure_dat) < 3) {
        cat('WARNING: < 3 instruments. Cannot run MR-Egger or weighted median.\n')
    }

    return(exposure_dat)
}

exposure_dat <- check_instrument_strength(exposure_dat)
```

## Step 2: Mendelian Randomization

```r
outcome_dat <- read_outcome_data(
    filename = 'outcome_gwas.tsv',
    sep = '\t',
    snp_col = 'SNP',
    beta_col = 'BETA',
    se_col = 'SE',
    effect_allele_col = 'A1',
    other_allele_col = 'A2',
    eaf_col = 'EAF',
    pval_col = 'P'
)

# Harmonize exposure and outcome data (align alleles)
dat <- harmonise_data(exposure_dat, outcome_dat)

# Run multiple MR methods for triangulation
mr_results <- mr(dat, method_list = c(
    'mr_ivw',              # Primary method (assumes no pleiotropy)
    'mr_egger_regression', # Allows directional pleiotropy (needs >= 3 SNPs)
    'mr_weighted_median',  # Robust to up to 50% invalid instruments
    'mr_weighted_mode'     # Most robust, least powerful
))

print(mr_results[, c('method', 'nsnp', 'b', 'se', 'pval')])
```

### QC Checkpoint: MR Consistency

```r
check_mr_consistency <- function(mr_results) {
    # All methods should agree on direction of effect
    directions <- sign(mr_results$b)
    consistent <- length(unique(directions)) == 1

    if (!consistent) {
        cat('WARNING: Inconsistent effect directions across methods.\n')
        cat('This suggests potential pleiotropy or weak instruments.\n')
    } else {
        cat('OK: All methods agree on direction of effect.\n')
    }

    # IVW should be significant (primary test)
    ivw <- mr_results[mr_results$method == 'Inverse variance weighted', ]
    cat(sprintf('IVW: beta = %.3f, p = %.2e\n', ivw$b, ivw$pval))

    return(consistent)
}

check_mr_consistency(mr_results)
```

## Step 3: Sensitivity Analysis

```r
# MR-PRESSO: detect and correct for outlier instruments
library(MRPRESSO)

presso_results <- mr_presso(
    BetaOutcome = 'beta.outcome',
    BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome',
    SdExposure = 'se.exposure',
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = dat,
    NbDistribution = 3000,
    SignifThreshold = 0.05
)

# Global test: p > 0.05 means no significant outliers
cat(sprintf('MR-PRESSO global test p: %.4f\n', presso_results$`MR-PRESSO results`$`Global Test`$Pvalue))

# MR-Egger intercept test for directional pleiotropy
egger_intercept <- mr_pleiotropy_test(dat)
cat(sprintf('Egger intercept: %.4f, p: %.4f\n', egger_intercept$egger_intercept, egger_intercept$pval))
# p > 0.05: No evidence of directional pleiotropy (good)

# Cochran Q test for heterogeneity
het <- mr_heterogeneity(dat)
cat(sprintf('Cochran Q p: %.4f\n', het$Q_pval[het$method == 'Inverse variance weighted']))

# Leave-one-out analysis
loo <- mr_leaveoneout(dat)

# Steiger directionality test
steiger <- directionality_test(dat)
cat(sprintf('Steiger: correct direction = %s, p = %.2e\n', steiger$correct_causal_direction, steiger$steiger_pval))
```

### QC Checkpoint: Sensitivity

```r
check_sensitivity <- function(presso_pval, egger_intercept_pval, steiger_correct) {
    issues <- c()

    if (presso_pval < 0.05) {
        issues <- c(issues, 'MR-PRESSO detects outliers - use corrected estimate')
    }
    if (egger_intercept_pval < 0.05) {
        issues <- c(issues, 'Egger intercept significant - directional pleiotropy present')
    }
    if (!steiger_correct) {
        issues <- c(issues, 'Steiger test suggests reverse causation')
    }

    if (length(issues) == 0) {
        cat('OK: All sensitivity tests passed.\n')
    } else {
        cat('WARNINGS:\n')
        for (issue in issues) cat(sprintf('  - %s\n', issue))
    }

    return(length(issues) == 0)
}
```

## Step 4: Colocalization

```r
library(coloc)

# Prepare datasets (need full summary statistics for the locus)
dataset1 <- list(
    beta = exposure_locus$BETA,
    varbeta = exposure_locus$SE^2,
    snp = exposure_locus$SNP,
    position = exposure_locus$BP,
    type = 'quant',
    N = exposure_n,
    MAF = exposure_locus$EAF
)

dataset2 <- list(
    beta = outcome_locus$BETA,
    varbeta = outcome_locus$SE^2,
    snp = outcome_locus$SNP,
    position = outcome_locus$BP,
    type = 'cc',
    N = outcome_n,
    s = case_fraction,
    MAF = outcome_locus$EAF
)

# Run colocalization
# Priors: p1=1e-4, p2=1e-4 (standard), p12=1e-5 (conservative)
# Use p12=5e-6 for stricter colocalization
result <- coloc.abf(dataset1, dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

cat(sprintf('PP.H0 (neither): %.3f\n', result$summary['PP.H0.abf']))
cat(sprintf('PP.H1 (exposure only): %.3f\n', result$summary['PP.H1.abf']))
cat(sprintf('PP.H2 (outcome only): %.3f\n', result$summary['PP.H2.abf']))
cat(sprintf('PP.H3 (both, different variants): %.3f\n', result$summary['PP.H3.abf']))
cat(sprintf('PP.H4 (shared causal variant): %.3f\n', result$summary['PP.H4.abf']))

# PP.H4 > 0.8: Strong evidence for shared causal variant
# PP.H4 > 0.5: Suggestive colocalization
# PP.H3 > 0.5: Different causal variants at same locus
```

## Step 5: Fine-Mapping with SuSiE

```r
library(susieR)

# Prepare LD matrix for the locus (from reference panel or in-sample)
# R: correlation matrix from genotype data
R <- as.matrix(read.csv('ld_matrix.csv', row.names = 1))

fitted <- susie_rss(
    bhat = locus_stats$BETA,
    shat = locus_stats$SE,
    R = R,
    n = sample_size,
    L = 10,  # Maximum number of causal variants (10 is standard)
    coverage = 0.95  # 95% credible set coverage
)

# Extract credible sets
cs <- fitted$sets$cs
for (i in seq_along(cs)) {
    snps_in_cs <- locus_stats$SNP[cs[[i]]]
    pip <- fitted$pip[cs[[i]]]
    cat(sprintf('Credible set %d: %d SNPs, min PIP = %.3f\n', i, length(snps_in_cs), min(pip)))
}

# SNPs with PIP > 0.5 are strong causal candidates
high_pip <- locus_stats$SNP[fitted$pip > 0.5]
cat(sprintf('High-confidence causal SNPs (PIP > 0.5): %s\n', paste(high_pip, collapse = ', ')))
```

## Step 6: Mediation Analysis with MVMR

```r
library(MendelianRandomization)

# Multivariable MR: test if mediator explains exposure-outcome relationship
# Example: BMI -> CRP -> CHD (is CRP a mediator?)

# Get instruments for all exposures
exposure_bmi <- mv_extract_exposures(c('ieu-a-2', 'ieu-a-1089'))  # BMI and CRP from IEU GWAS DB

outcome_chd <- mv_extract_outcome(exposure_bmi$SNP, 'ieu-a-7')

mvdat <- mv_harmonise_data(exposure_bmi, outcome_chd)
mvmr_result <- mv_multiple(mvdat)

cat('MVMR results (direct effects):\n')
print(mvmr_result$result[, c('exposure', 'b', 'se', 'pval')])

# If direct effect of exposure attenuates after adjusting for mediator,
# the mediator partially explains the causal path
```

## Complete Pipeline Script

```r
library(TwoSampleMR)
library(MRPRESSO)
library(coloc)
library(susieR)

run_causal_pipeline <- function(exposure_file, outcome_file, locus_chr, locus_start, locus_end,
                                 ld_matrix_file = NULL, exposure_n = NULL, outcome_n = NULL,
                                 case_fraction = NULL) {

    cat('=== Step 1: Instrument Selection ===\n')
    exposure_dat <- read_exposure_data(
        filename = exposure_file, sep = '\t',
        snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
        effect_allele_col = 'A1', other_allele_col = 'A2',
        eaf_col = 'EAF', pval_col = 'P'
    )
    exposure_dat <- subset(exposure_dat, pval.exposure < 5e-8)
    exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.001, clump_kb = 10000)
    exposure_dat$F_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
    exposure_dat <- subset(exposure_dat, F_stat >= 10)
    cat(sprintf('Selected %d instruments (mean F = %.1f)\n', nrow(exposure_dat), mean(exposure_dat$F_stat)))

    cat('\n=== Step 2: Mendelian Randomization ===\n')
    outcome_dat <- read_outcome_data(
        filename = outcome_file, sep = '\t',
        snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
        effect_allele_col = 'A1', other_allele_col = 'A2',
        eaf_col = 'EAF', pval_col = 'P'
    )
    dat <- harmonise_data(exposure_dat, outcome_dat)
    mr_results <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression',
                                           'mr_weighted_median', 'mr_weighted_mode'))
    print(mr_results[, c('method', 'nsnp', 'b', 'se', 'pval')])

    cat('\n=== Step 3: Sensitivity Analysis ===\n')
    presso <- mr_presso(
        BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure',
        SdOutcome = 'se.outcome', SdExposure = 'se.exposure',
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = dat, NbDistribution = 3000, SignifThreshold = 0.05
    )
    egger_int <- mr_pleiotropy_test(dat)
    steiger <- directionality_test(dat)
    cat(sprintf('MR-PRESSO global p: %.4f\n', presso$`MR-PRESSO results`$`Global Test`$Pvalue))
    cat(sprintf('Egger intercept p: %.4f\n', egger_int$pval))
    cat(sprintf('Steiger correct direction: %s\n', steiger$correct_causal_direction))

    cat('\n=== Step 4: Colocalization ===\n')
    if (!is.null(exposure_n) && !is.null(outcome_n)) {
        exposure_locus <- read.delim(exposure_file)
        exposure_locus <- subset(exposure_locus, CHR == locus_chr & BP >= locus_start & BP <= locus_end)
        outcome_locus <- read.delim(outcome_file)
        outcome_locus <- subset(outcome_locus, CHR == locus_chr & BP >= locus_start & BP <= locus_end)

        common_snps <- intersect(exposure_locus$SNP, outcome_locus$SNP)
        exposure_locus <- exposure_locus[match(common_snps, exposure_locus$SNP), ]
        outcome_locus <- outcome_locus[match(common_snps, outcome_locus$SNP), ]

        d1 <- list(beta = exposure_locus$BETA, varbeta = exposure_locus$SE^2,
                    snp = exposure_locus$SNP, position = exposure_locus$BP,
                    type = 'quant', N = exposure_n, MAF = exposure_locus$EAF)
        d2 <- list(beta = outcome_locus$BETA, varbeta = outcome_locus$SE^2,
                    snp = outcome_locus$SNP, position = outcome_locus$BP,
                    type = 'cc', N = outcome_n, s = case_fraction, MAF = outcome_locus$EAF)

        coloc_result <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
        cat(sprintf('PP.H4 (shared causal variant): %.3f\n', coloc_result$summary['PP.H4.abf']))
    }

    cat('\n=== Step 5: Fine-Mapping ===\n')
    if (!is.null(ld_matrix_file)) {
        R <- as.matrix(read.csv(ld_matrix_file, row.names = 1))
        fitted <- susie_rss(
            bhat = exposure_locus$BETA, shat = exposure_locus$SE,
            R = R, n = exposure_n, L = 10, coverage = 0.95
        )
        high_pip <- exposure_locus$SNP[fitted$pip > 0.5]
        cat(sprintf('Causal SNPs (PIP > 0.5): %s\n', paste(high_pip, collapse = ', ')))
    }

    cat('\n=== Pipeline Complete ===\n')
    return(list(mr = mr_results, harmonised = dat, sensitivity = list(
        presso = presso, egger_intercept = egger_int, steiger = steiger
    )))
}
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Instruments | p-value threshold | 5e-8 (standard), 5e-6 for underpowered |
| Instruments | clump_r2 | 0.001 (strict), 0.01 for more instruments |
| Instruments | F-statistic | > 10 required, > 30 preferred |
| Coloc | p12 prior | 1e-5 (standard), 5e-6 (conservative) |
| Coloc | PP.H4 threshold | > 0.8 strong, > 0.5 suggestive |
| SuSiE | L (max causal) | 10 (default), reduce for simple loci |
| SuSiE | coverage | 0.95 (95% credible set) |
| MR-PRESSO | NbDistribution | 3000 (standard), 10000 for publication |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| No instruments | Underpowered GWAS | Relax p-value to 5e-6 (with caution) |
| Weak instruments (F < 10) | Small effect SNPs | Remove weak instruments; use more powered GWAS |
| Inconsistent MR methods | Pleiotropy | Check MR-PRESSO outliers; use weighted median |
| Egger intercept p < 0.05 | Directional pleiotropy | Report Egger estimate as primary |
| PP.H3 > PP.H4 | Different causal variants | Loci associate through different mechanisms |
| No credible sets | LD matrix issues | Verify LD matrix matches summary stats population |
| Steiger reverse direction | Reverse causation | Consider bidirectional MR |

## Related Skills

- causal-genomics/mendelian-randomization - MR method details
- causal-genomics/colocalization-analysis - Coloc implementation
- causal-genomics/fine-mapping - SuSiE and FINEMAP details
- causal-genomics/pleiotropy-detection - MR-PRESSO and sensitivity tests
- causal-genomics/mediation-analysis - MVMR and network MR
- population-genetics/association-testing - Upstream GWAS methods
