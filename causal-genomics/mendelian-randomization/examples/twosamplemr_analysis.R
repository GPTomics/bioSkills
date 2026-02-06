## Two-sample Mendelian randomization with TwoSampleMR
##
## Demonstrates: instrument selection, clumping, harmonization,
## MR analysis with multiple methods, and sensitivity checks.
## Uses local GWAS summary statistics (no OpenGWAS auth needed).

library(TwoSampleMR)

# --- Simulate GWAS summary statistics for demonstration ---
set.seed(42)
n_snps <- 50

exposure_gwas <- data.frame(
  SNP = paste0('rs', 1:n_snps),
  BETA = rnorm(n_snps, 0.05, 0.02),
  SE = runif(n_snps, 0.005, 0.015),
  A1 = sample(c('A', 'C', 'G', 'T'), n_snps, replace = TRUE),
  A2 = sample(c('A', 'C', 'G', 'T'), n_snps, replace = TRUE),
  EAF = runif(n_snps, 0.1, 0.9),
  stringsAsFactors = FALSE
)
exposure_gwas$P <- 2 * pnorm(-abs(exposure_gwas$BETA / exposure_gwas$SE))
# Keep genome-wide significant
exposure_gwas <- exposure_gwas[exposure_gwas$P < 5e-08, ]

outcome_gwas <- data.frame(
  SNP = exposure_gwas$SNP,
  BETA = exposure_gwas$BETA * 0.3 + rnorm(nrow(exposure_gwas), 0, 0.01),
  SE = runif(nrow(exposure_gwas), 0.008, 0.02),
  A1 = exposure_gwas$A1,
  A2 = exposure_gwas$A2,
  EAF = exposure_gwas$EAF + rnorm(nrow(exposure_gwas), 0, 0.02),
  stringsAsFactors = FALSE
)
outcome_gwas$P <- 2 * pnorm(-abs(outcome_gwas$BETA / outcome_gwas$SE))

# Write temporary files
write.table(exposure_gwas, 'exposure_gwas.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(outcome_gwas, 'outcome_gwas.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# --- Step 1: Read exposure data ---
exposure_dat <- read_exposure_data(
  filename = 'exposure_gwas.txt', sep = '\t',
  snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
  effect_allele_col = 'A1', other_allele_col = 'A2',
  pval_col = 'P', eaf_col = 'EAF'
)

# --- Step 2: Instrument strength ---
# F = (beta / se)^2; F > 10 required for valid instruments
exposure_dat$f_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
cat('Mean F-statistic:', mean(exposure_dat$f_stat), '\n')
cat('Weak instruments (F < 10):', sum(exposure_dat$f_stat < 10), '\n')
exposure_dat <- exposure_dat[exposure_dat$f_stat >= 10, ]

# --- Step 3: Read outcome data ---
outcome_dat <- read_outcome_data(
  filename = 'outcome_gwas.txt', sep = '\t',
  snps = exposure_dat$SNP,
  snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
  effect_allele_col = 'A1', other_allele_col = 'A2',
  pval_col = 'P', eaf_col = 'EAF'
)

# --- Step 4: Harmonize ---
# action = 2: Infer forward strand alleles (recommended)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
cat('SNPs after harmonization:', nrow(dat), '\n')

# --- Step 5: MR analysis ---
methods <- c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_weighted_mode')
results <- mr(dat, method_list = methods)
print(results[, c('method', 'nsnp', 'b', 'se', 'pval')])

# --- Step 6: Sensitivity analyses ---
het <- mr_heterogeneity(dat)
cat('\nHeterogeneity (Cochran Q):\n')
print(het[, c('method', 'Q', 'Q_pval')])

pleio <- mr_pleiotropy_test(dat)
cat('\nEgger intercept test:\n')
cat('  Intercept:', pleio$egger_intercept, '\n')
cat('  P-value:', pleio$pval, '\n')

# --- Step 7: Steiger directionality ---
steiger <- directionality_test(dat)
cat('\nSteiger test:\n')
cat('  Correct direction:', steiger$correct_causal_direction, '\n')
cat('  Steiger p-value:', steiger$steiger_pval, '\n')

# Clean up
file.remove('exposure_gwas.txt', 'outcome_gwas.txt')
