## Comprehensive MR sensitivity analyses
##
## Runs IVW, MR-Egger, weighted median, weighted mode, and MR-PRESSO
## on the same dataset and produces a comparison table plus diagnostics.

library(TwoSampleMR)
library(MRPRESSO)

# --- Simulate harmonized TwoSampleMR data ---
set.seed(42)
n <- 25

dat <- data.frame(
  SNP = paste0('rs', 1:n),
  beta.exposure = abs(rnorm(n, 0.06, 0.02)),
  se.exposure = rep(0.012, n),
  beta.outcome = NA,
  se.outcome = rep(0.015, n),
  effect_allele.exposure = rep('A', n),
  other_allele.exposure = rep('G', n),
  effect_allele.outcome = rep('A', n),
  other_allele.outcome = rep('G', n),
  eaf.exposure = runif(n, 0.15, 0.85),
  eaf.outcome = runif(n, 0.15, 0.85),
  id.exposure = rep('exposure', n),
  id.outcome = rep('outcome', n),
  exposure = rep('BMI', n),
  outcome = rep('T2D', n),
  mr_keep = rep(TRUE, n),
  stringsAsFactors = FALSE
)

# Outcome = causal effect + noise + some pleiotropy
dat$beta.outcome <- dat$beta.exposure * 0.35 + rnorm(n, 0, 0.006)
dat$beta.outcome[c(3, 12)] <- dat$beta.outcome[c(3, 12)] + c(0.04, -0.03)

# --- F-statistics ---
dat$f_stat <- (dat$beta.exposure / dat$se.exposure)^2
cat('Mean F-statistic:', round(mean(dat$f_stat), 1), '\n')
cat('Weak instruments (F < 10):', sum(dat$f_stat < 10), '\n\n')

# --- 1. Multiple MR methods ---
methods <- c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_weighted_mode')
mr_results <- mr(dat, method_list = methods)
cat('=== MR Method Comparison ===\n')
print(mr_results[, c('method', 'nsnp', 'b', 'se', 'pval')])

# --- 2. Heterogeneity ---
het <- mr_heterogeneity(dat)
cat('\n=== Heterogeneity ===\n')
print(het[, c('method', 'Q', 'Q_df', 'Q_pval')])

# --- 3. Egger intercept ---
pleio <- mr_pleiotropy_test(dat)
cat('\n=== Egger Intercept ===\n')
cat('Intercept:', round(pleio$egger_intercept, 5), '\n')
cat('SE:', round(pleio$se, 5), '\n')
cat('P-value:', format.pval(pleio$pval), '\n')

# --- 4. I-squared (NOME assumption for Egger) ---
# I^2 > 0.9: MR-Egger estimate is reliable
# I^2 < 0.6: Egger has low power, NOME assumption violated
isq <- Isq(dat$beta.exposure, dat$se.exposure)
cat('\nI-squared:', round(isq, 3), '\n')
if (isq < 0.9) {
  cat('Warning: I-squared < 0.9; MR-Egger may be unreliable\n')
} else {
  cat('I-squared adequate for MR-Egger\n')
}

# --- 5. MR-PRESSO ---
presso_input <- data.frame(
  bx = dat$beta.exposure, by = dat$beta.outcome,
  bxse = dat$se.exposure, byse = dat$se.outcome
)
presso <- mr_presso(
  BetaOutcome = 'by', BetaExposure = 'bx',
  SdOutcome = 'byse', SdExposure = 'bxse',
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = presso_input, NbDistribution = 5000, SignifThreshold = 0.05
)

cat('\n=== MR-PRESSO ===\n')
cat('Global test p-value:', presso$`MR-PRESSO results`$`Global Test`$Pvalue, '\n')
outlier_p <- presso$`MR-PRESSO results`$`Outlier Test`$Pvalue
cat('Outliers detected:', sum(outlier_p < 0.05, na.rm = TRUE), '\n')

# --- 6. Steiger directionality ---
steiger <- directionality_test(dat)
cat('\n=== Steiger Directionality ===\n')
cat('Correct direction:', steiger$correct_causal_direction, '\n')
cat('P-value:', format.pval(steiger$steiger_pval), '\n')

# --- 7. Leave-one-out ---
loo <- mr_leaveoneout(dat)
loo_range <- range(loo$b[!is.na(loo$b)])
cat('\n=== Leave-One-Out ===\n')
cat('Estimate range:', round(loo_range[1], 4), 'to', round(loo_range[2], 4), '\n')
cat('Full IVW estimate:', round(mr_results$b[mr_results$method == 'Inverse variance weighted'], 4), '\n')

# --- Summary ---
cat('\n=== Summary ===\n')
cat('All methods agree on direction: ', all(mr_results$b > 0) || all(mr_results$b < 0), '\n')
cat('Significant heterogeneity: ', any(het$Q_pval < 0.05), '\n')
cat('Significant Egger intercept: ', pleio$pval < 0.05, '\n')
cat('MR-PRESSO outliers: ', sum(outlier_p < 0.05, na.rm = TRUE), '\n')
cat('Steiger correct direction: ', steiger$correct_causal_direction, '\n')
