# Reference: TwoSampleMR 0.5+ | Verify API if version differs
## MR visualization with TwoSampleMR
##
## Demonstrates scatter, forest, leave-one-out, and funnel plots.
## Requires a harmonized dat object from twosamplemr_analysis.R.

library(TwoSampleMR)
library(ggplot2)

# --- Simulate harmonized data for standalone demo ---
set.seed(42)
n <- 20

exposure_dat <- data.frame(
  SNP = paste0('rs', 1:n),
  beta.exposure = rnorm(n, 0.05, 0.02),
  se.exposure = runif(n, 0.005, 0.012),
  effect_allele.exposure = rep('A', n),
  other_allele.exposure = rep('G', n),
  eaf.exposure = runif(n, 0.15, 0.85),
  pval.exposure = rep(1e-10, n),
  exposure = rep('BMI', n),
  id.exposure = rep('exposure', n),
  mr_keep.exposure = rep(TRUE, n),
  stringsAsFactors = FALSE
)

outcome_dat <- data.frame(
  SNP = paste0('rs', 1:n),
  beta.outcome = exposure_dat$beta.exposure * 0.4 + rnorm(n, 0, 0.008),
  se.outcome = runif(n, 0.008, 0.018),
  effect_allele.outcome = rep('A', n),
  other_allele.outcome = rep('G', n),
  eaf.outcome = exposure_dat$eaf.exposure,
  pval.outcome = rep(0.01, n),
  outcome = rep('T2D', n),
  id.outcome = rep('outcome', n),
  mr_keep.outcome = rep(TRUE, n),
  stringsAsFactors = FALSE
)

dat <- merge(exposure_dat, outcome_dat, by = 'SNP')
dat$mr_keep <- TRUE

# --- Run MR ---
results <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median'))
single <- mr_singlesnp(dat)
loo <- mr_leaveoneout(dat)

# --- Scatter plot ---
# Shows SNP-exposure effects (x) vs SNP-outcome effects (y)
# Each line represents a different MR method slope
p_scatter <- mr_scatter_plot(results, dat)
ggsave('mr_scatter.pdf', p_scatter[[1]], width = 8, height = 6)

# --- Forest plot ---
# Individual SNP Wald ratios with combined IVW/Egger estimates
p_forest <- mr_forest_plot(single)
ggsave('mr_forest.pdf', p_forest[[1]], width = 8, height = 8)

# --- Leave-one-out plot ---
# Causal estimate when each SNP is removed; checks for influential outliers
p_loo <- mr_leaveoneout_plot(loo)
ggsave('mr_leaveoneout.pdf', p_loo[[1]], width = 8, height = 8)

# --- Funnel plot ---
# Precision (1/SE) vs individual causal estimates
# Asymmetry suggests directional pleiotropy
p_funnel <- mr_funnel_plot(single)
ggsave('mr_funnel.pdf', p_funnel[[1]], width = 8, height = 6)

cat('Plots saved: mr_scatter.pdf, mr_forest.pdf, mr_leaveoneout.pdf, mr_funnel.pdf\n')
