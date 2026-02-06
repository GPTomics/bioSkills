# Complete post-GWAS causal inference pipeline
# Mendelian randomization -> sensitivity analysis -> colocalization -> fine-mapping

library(TwoSampleMR)
library(MRPRESSO)
library(coloc)
library(susieR)

# Configuration
# Example: LDL cholesterol (exposure) -> coronary heart disease (outcome)
# Use GWAS summary statistics from public repositories:
# - IEU OpenGWAS: https://gwas.mrcieu.ac.uk/
# - GWAS Catalog: https://www.ebi.ac.uk/gwas/
EXPOSURE_FILE <- 'exposure_gwas.tsv'
OUTCOME_FILE <- 'outcome_gwas.tsv'
EXPOSURE_N <- 100000
OUTCOME_N <- 50000
CASE_FRACTION <- 0.3

# ============================================================
# Step 1: Instrument Selection
# ============================================================

cat('=== Step 1: Instrument Selection ===\n')

exposure_dat <- read_exposure_data(
    filename = EXPOSURE_FILE,
    sep = '\t',
    snp_col = 'SNP',
    beta_col = 'BETA',
    se_col = 'SE',
    effect_allele_col = 'A1',
    other_allele_col = 'A2',
    eaf_col = 'EAF',
    pval_col = 'P'
)

# p < 5e-8: Standard GWAS significance; use 5e-6 for underpowered exposures
exposure_dat <- subset(exposure_dat, pval.exposure < 5e-8)
cat(sprintf('Genome-wide significant SNPs: %d\n', nrow(exposure_dat)))

# LD clumping for independent instruments
# clump_r2 = 0.001: Strict independence (standard for MR)
# clump_kb = 10000: 10 Mb window
exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.001, clump_kb = 10000)
cat(sprintf('After LD clumping: %d independent instruments\n', nrow(exposure_dat)))

# F-statistic > 10 to avoid weak instrument bias
exposure_dat$F_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
cat(sprintf('F-statistic range: %.1f - %.1f (mean: %.1f)\n',
            min(exposure_dat$F_stat), max(exposure_dat$F_stat), mean(exposure_dat$F_stat)))

weak <- sum(exposure_dat$F_stat < 10)
if (weak > 0) {
    cat(sprintf('WARNING: Removing %d weak instruments (F < 10)\n', weak))
    exposure_dat <- subset(exposure_dat, F_stat >= 10)
}

cat(sprintf('Final instrument count: %d\n', nrow(exposure_dat)))

# ============================================================
# Step 2: Mendelian Randomization
# ============================================================

cat('\n=== Step 2: Mendelian Randomization ===\n')

outcome_dat <- read_outcome_data(
    filename = OUTCOME_FILE,
    sep = '\t',
    snp_col = 'SNP',
    beta_col = 'BETA',
    se_col = 'SE',
    effect_allele_col = 'A1',
    other_allele_col = 'A2',
    eaf_col = 'EAF',
    pval_col = 'P'
)

dat <- harmonise_data(exposure_dat, outcome_dat)
cat(sprintf('Harmonised SNPs: %d\n', nrow(subset(dat, mr_keep == TRUE))))

# Run multiple MR methods for triangulation
# IVW: Primary (assumes balanced pleiotropy)
# Egger: Allows directional pleiotropy (needs >= 3 SNPs)
# Weighted median: Robust to up to 50% invalid instruments
# Weighted mode: Most robust, least power
mr_results <- mr(dat, method_list = c(
    'mr_ivw',
    'mr_egger_regression',
    'mr_weighted_median',
    'mr_weighted_mode'
))

cat('\nMR Results:\n')
print(mr_results[, c('method', 'nsnp', 'b', 'se', 'pval')])

# QC: Check consistency across methods
directions <- sign(mr_results$b)
if (length(unique(directions)) == 1) {
    cat('OK: All methods agree on direction of effect\n')
} else {
    cat('WARNING: Inconsistent effect directions - investigate pleiotropy\n')
}

# ============================================================
# Step 3: Sensitivity Analysis
# ============================================================

cat('\n=== Step 3: Sensitivity Analysis ===\n')

# MR-PRESSO: detect and correct for outlier instruments
# NbDistribution=3000: Standard; use 10000 for publication
presso <- mr_presso(
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

presso_global_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
cat(sprintf('MR-PRESSO global test p: %.4f', presso_global_p))
if (presso_global_p < 0.05) {
    cat(' (outliers detected)\n')
    n_outliers <- sum(presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, na.rm = TRUE)
    cat(sprintf('Outlier instruments: %d\n', n_outliers))
} else {
    cat(' (no significant outliers)\n')
}

# MR-Egger intercept for directional pleiotropy
# p > 0.05: No evidence of directional pleiotropy (good)
egger_int <- mr_pleiotropy_test(dat)
cat(sprintf('Egger intercept: %.4f, p: %.4f', egger_int$egger_intercept, egger_int$pval))
if (egger_int$pval < 0.05) {
    cat(' (directional pleiotropy present!)\n')
} else {
    cat(' (no directional pleiotropy)\n')
}

# Cochran Q for heterogeneity
het <- mr_heterogeneity(dat)
ivw_het <- het[het$method == 'Inverse variance weighted', ]
cat(sprintf('Cochran Q p: %.4f\n', ivw_het$Q_pval))

# Steiger directionality test
steiger <- directionality_test(dat)
cat(sprintf('Steiger correct direction: %s (p: %.2e)\n',
            steiger$correct_causal_direction, steiger$steiger_pval))

if (!steiger$correct_causal_direction) {
    cat('WARNING: Steiger suggests reverse causation. Consider bidirectional MR.\n')
}

# Leave-one-out analysis
loo <- mr_leaveoneout(dat)
loo_range <- range(loo$b[loo$SNP != 'All'])
cat(sprintf('Leave-one-out estimate range: [%.3f, %.3f]\n', loo_range[1], loo_range[2]))

# ============================================================
# Step 4: Colocalization
# ============================================================

cat('\n=== Step 4: Colocalization ===\n')

# Read full summary statistics for the top locus
# Identify top MR locus from strongest instrument
top_snp <- exposure_dat[which.max(exposure_dat$F_stat), ]
cat(sprintf('Top locus: %s (chr%s)\n', top_snp$SNP, top_snp$chr.exposure))

# Extract +/- 500 kb around top SNP for colocalization
exposure_full <- read.delim(EXPOSURE_FILE)
outcome_full <- read.delim(OUTCOME_FILE)

locus_chr <- top_snp$chr.exposure
locus_pos <- top_snp$pos.exposure
# 500 kb window: Standard for colocalization
window <- 500000

exposure_locus <- subset(exposure_full, CHR == locus_chr & BP >= (locus_pos - window) & BP <= (locus_pos + window))
outcome_locus <- subset(outcome_full, CHR == locus_chr & BP >= (locus_pos - window) & BP <= (locus_pos + window))

common_snps <- intersect(exposure_locus$SNP, outcome_locus$SNP)
exposure_locus <- exposure_locus[match(common_snps, exposure_locus$SNP), ]
outcome_locus <- outcome_locus[match(common_snps, outcome_locus$SNP), ]
cat(sprintf('SNPs in locus: %d\n', length(common_snps)))

d1 <- list(
    beta = exposure_locus$BETA, varbeta = exposure_locus$SE^2,
    snp = exposure_locus$SNP, position = exposure_locus$BP,
    type = 'quant', N = EXPOSURE_N, MAF = exposure_locus$EAF
)
d2 <- list(
    beta = outcome_locus$BETA, varbeta = outcome_locus$SE^2,
    snp = outcome_locus$SNP, position = outcome_locus$BP,
    type = 'cc', N = OUTCOME_N, s = CASE_FRACTION, MAF = outcome_locus$EAF
)

# Priors: p1=1e-4, p2=1e-4, p12=1e-5 (standard)
# Use p12=5e-6 for more conservative colocalization
coloc_result <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

cat('\nColocalization posterior probabilities:\n')
cat(sprintf('  PP.H0 (neither): %.3f\n', coloc_result$summary['PP.H0.abf']))
cat(sprintf('  PP.H1 (exposure only): %.3f\n', coloc_result$summary['PP.H1.abf']))
cat(sprintf('  PP.H2 (outcome only): %.3f\n', coloc_result$summary['PP.H2.abf']))
cat(sprintf('  PP.H3 (both, different variants): %.3f\n', coloc_result$summary['PP.H3.abf']))
cat(sprintf('  PP.H4 (shared causal variant): %.3f\n', coloc_result$summary['PP.H4.abf']))

# PP.H4 > 0.8: Strong evidence for shared causal variant
# PP.H4 > 0.5: Suggestive
pp4 <- coloc_result$summary['PP.H4.abf']
if (pp4 > 0.8) {
    cat('Strong colocalization evidence\n')
} else if (pp4 > 0.5) {
    cat('Suggestive colocalization\n')
} else {
    cat('Weak colocalization. Traits may have different causal variants at this locus.\n')
}

# ============================================================
# Step 5: Fine-Mapping with SuSiE
# ============================================================

cat('\n=== Step 5: Fine-Mapping ===\n')

# Fine-mapping requires an LD matrix from a matched reference panel
# In practice: compute from 1000 Genomes or UK Biobank
# LD_MATRIX_FILE <- 'ld_matrix.csv'
# R <- as.matrix(read.csv(LD_MATRIX_FILE, row.names = 1))

# If LD matrix is available:
# L=10: Maximum number of causal variants (standard)
# coverage=0.95: 95% credible set
# fitted <- susie_rss(
#     bhat = exposure_locus$BETA,
#     shat = exposure_locus$SE,
#     R = R,
#     n = EXPOSURE_N,
#     L = 10,
#     coverage = 0.95
# )
#
# cs <- fitted$sets$cs
# for (i in seq_along(cs)) {
#     snps_in_cs <- exposure_locus$SNP[cs[[i]]]
#     pip <- fitted$pip[cs[[i]]]
#     cat(sprintf('Credible set %d: %d SNPs, min PIP = %.3f\n', i, length(snps_in_cs), min(pip)))
# }
#
# high_pip <- exposure_locus$SNP[fitted$pip > 0.5]
# cat(sprintf('Causal SNPs (PIP > 0.5): %s\n', paste(high_pip, collapse = ', ')))

cat('Fine-mapping requires LD matrix from matched reference panel.\n')
cat('Compute from 1000 Genomes or in-sample genotypes, then run susie_rss().\n')

# ============================================================
# Summary
# ============================================================

cat('\n=== Causal Evidence Summary ===\n')
ivw <- mr_results[mr_results$method == 'Inverse variance weighted', ]
cat(sprintf('IVW estimate: beta = %.3f (95%% CI: %.3f to %.3f), p = %.2e\n',
            ivw$b, ivw$b - 1.96 * ivw$se, ivw$b + 1.96 * ivw$se, ivw$pval))
cat(sprintf('Direction consistency: %s\n', ifelse(length(unique(sign(mr_results$b))) == 1, 'YES', 'NO')))
cat(sprintf('Pleiotropy (Egger intercept p): %.4f\n', egger_int$pval))
cat(sprintf('Outliers (MR-PRESSO p): %.4f\n', presso_global_p))
cat(sprintf('Correct direction (Steiger): %s\n', steiger$correct_causal_direction))
cat(sprintf('Colocalization (PP.H4): %.3f\n', pp4))

cat('\nPipeline complete.\n')
