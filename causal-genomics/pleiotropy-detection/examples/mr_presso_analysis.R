## MR-PRESSO outlier detection and correction
##
## Detects horizontal pleiotropy among MR instruments, identifies
## outlier SNPs, and compares raw vs corrected causal estimates.

library(MRPRESSO)

# --- Simulate harmonized MR data ---
set.seed(42)
n_instruments <- 30

# True causal effect = 0.3
beta_exposure <- abs(rnorm(n_instruments, 0.05, 0.02))
se_exposure <- rep(0.01, n_instruments)

# Outcome effects: causal + noise + some pleiotropic outliers
beta_outcome <- beta_exposure * 0.3 + rnorm(n_instruments, 0, 0.005)
se_outcome <- rep(0.012, n_instruments)

# Introduce 3 pleiotropic outliers (large direct effects on outcome)
outlier_idx <- c(5, 15, 25)
beta_outcome[outlier_idx] <- beta_outcome[outlier_idx] + c(0.04, -0.03, 0.05)

dat <- data.frame(
  bx = beta_exposure,
  by = beta_outcome,
  bxse = se_exposure,
  byse = se_outcome,
  snp = paste0('rs', 1:n_instruments)
)

# --- Run MR-PRESSO ---
# NbDistribution: Number of random samplings for null distribution
# 5000: Recommended for publication (1000 minimum)
# SignifThreshold: P-value cutoff for outlier detection
presso <- mr_presso(
  BetaOutcome = 'by', BetaExposure = 'bx',
  SdOutcome = 'byse', SdExposure = 'bxse',
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = dat,
  NbDistribution = 5000,
  SignifThreshold = 0.05
)

# --- Global test ---
global_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
cat('=== MR-PRESSO Results ===\n\n')
cat('Global test p-value:', global_p, '\n')
if (global_p < 0.05) {
  cat('Interpretation: Significant pleiotropy detected among instruments\n\n')
} else {
  cat('Interpretation: No significant evidence of pleiotropy\n\n')
}

# --- Outlier test ---
outliers <- presso$`MR-PRESSO results`$`Outlier Test`
cat('Outlier test:\n')
detected <- which(outliers$Pvalue < 0.05)
cat('  Detected outliers:', length(detected), '\n')
if (length(detected) > 0) {
  cat('  Outlier SNPs:', paste(dat$snp[detected], collapse = ', '), '\n')
  cat('  True outliers were:', paste(dat$snp[outlier_idx], collapse = ', '), '\n')
}

# --- Distortion test ---
distortion_p <- presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
cat('\nDistortion test p-value:', distortion_p, '\n')
if (!is.na(distortion_p) && distortion_p < 0.05) {
  cat('Outliers significantly distorted the causal estimate\n')
}

# --- Compare raw vs corrected ---
main <- presso$`Main MR results`
cat('\nCausal estimates:\n')
cat('  Raw IVW:', round(main$`Causal Estimate`[1], 4), '\n')
cat('  Corrected IVW:', round(main$`Causal Estimate`[2], 4), '\n')
cat('  True causal effect: 0.3\n')
cat('  Raw p-value:', format.pval(main$`P-value`[1]), '\n')
cat('  Corrected p-value:', format.pval(main$`P-value`[2]), '\n')
