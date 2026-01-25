# RNA-seq Power Analysis Examples
# Demonstrates power calculations for sequencing experiments

library(RNASeqPower)

# =============================================================================
# Basic Power Calculation
# =============================================================================

# Parameters explained:
# - depth: average read count per gene (20 is typical for well-expressed genes)
# - n: samples per group
# - cv: coefficient of variation (biological variability)
#   * Cell lines: 0.1-0.2
#   * Inbred mice: 0.2-0.3
#   * Human samples: 0.3-0.5
# - effect: fold change to detect (2 = 2-fold change)
# - alpha: significance level (0.05 standard)

# Calculate power with 3 replicates to detect 2-fold change
power_3rep <- rnapower(depth = 20, n = 3, cv = 0.4, effect = 2, alpha = 0.05)
cat('Power with n=3:', round(power_3rep, 3), '\n')

# Calculate power with 6 replicates
power_6rep <- rnapower(depth = 20, n = 6, cv = 0.4, effect = 2, alpha = 0.05)
cat('Power with n=6:', round(power_6rep, 3), '\n')

# =============================================================================
# Sample Size Estimation
# =============================================================================

# How many samples needed for 80% power to detect 2-fold change?
# Returns NA for n, need to solve iteratively
find_n_for_power <- function(target_power = 0.8, depth = 20, cv = 0.4,
                              effect = 2, alpha = 0.05) {
  for (n in 2:50) {
    p <- rnapower(depth = depth, n = n, cv = cv, effect = effect, alpha = alpha)
    if (p >= target_power) {
      return(n)
    }
  }
  return(NA)
}

n_needed <- find_n_for_power(target_power = 0.8, cv = 0.4, effect = 2)
cat('Samples needed for 80% power (2-fold, CV=0.4):', n_needed, '\n')

# =============================================================================
# Effect Size Sensitivity Analysis
# =============================================================================

# What effect sizes can we detect with different sample sizes?
effect_sizes <- c(1.25, 1.5, 2, 3, 4)
sample_sizes <- c(3, 5, 8, 12)

cat('\nPower table (CV=0.4, alpha=0.05, depth=20):\n')
cat('Effect\t', paste(paste0('n=', sample_sizes), collapse = '\t'), '\n')

for (effect in effect_sizes) {
  powers <- sapply(sample_sizes, function(n) {
    round(rnapower(depth = 20, n = n, cv = 0.4, effect = effect, alpha = 0.05), 2)
  })
  cat(paste0(effect, 'x\t'), paste(powers, collapse = '\t'), '\n')
}

# =============================================================================
# CV Impact Analysis
# =============================================================================

# How does biological variability affect required samples?
cv_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)

cat('\nSamples needed for 80% power to detect 2-fold change:\n')
for (cv in cv_values) {
  n <- find_n_for_power(target_power = 0.8, cv = cv, effect = 2)
  cat(sprintf('CV = %.1f: n = %d per group\n', cv, n))
}

# =============================================================================
# Sequencing Depth vs Sample Size
# =============================================================================

# More samples or deeper sequencing?
# Generally: more samples > deeper sequencing after ~20M reads

depths <- c(10, 20, 50, 100)
cat('\nPower with n=4, CV=0.4, effect=2 at different depths:\n')
for (d in depths) {
  p <- rnapower(depth = d, n = 4, cv = 0.4, effect = 2, alpha = 0.05)
  cat(sprintf('Depth %d: power = %.3f\n', d, p))
}

cat('\n# Key insight: Doubling depth from 20 to 40 gives less benefit\n')
cat('# than adding one more replicate (from n=4 to n=5)\n')
