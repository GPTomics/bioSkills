#!/usr/bin/env Rscript
# Reference: methylclock 1.8+ | Verify API if version differs
# Demonstrates the load-bearing clock mechanics on synthetic data (runs in base R):
#   1. a clock is a frozen weighted sum over a fixed CpG set,
#   2. missing clock CpGs imputed to the training mean bias age acceleration toward zero,
#   3. the endpoint is age ACCELERATION (residual of DNAm age on chronological age), not raw age.
# The methylclock / dnaMethyAge / DunedinPACE calls that do this on real arrays are shown in
# run_real_clocks() at the bottom; that path needs Bioconductor/GitHub packages installed.

set.seed(1)

n_samples <- 200
n_clock_cpgs <- 80

chronological_age <- runif(n_samples, 20, 80)

intercept <- 5
weights <- rnorm(n_clock_cpgs, 0, 0.4)

age_effect <- outer(chronological_age, rnorm(n_clock_cpgs, 0, 0.01))
beta <- plogis(matrix(rnorm(n_samples * n_clock_cpgs, 0, 0.3), n_samples) + age_effect)
colnames(beta) <- paste0('cg', sprintf('%05d', seq_len(n_clock_cpgs)))

apply_clock <- function(betas, w, b) as.numeric(b + betas %*% w)

raw_score <- scale(apply_clock(beta, weights, intercept))[, 1]
dnam_age <- mean(chronological_age) + 0.9 * (chronological_age - mean(chronological_age)) + 6 * raw_score

# Age ACCELERATION is the residual of DNAm age on chronological age (uncorrelated with age by construction).
eaa <- resid(lm(dnam_age ~ chronological_age))

cat('cor(DNAm age, chronological age):', round(cor(dnam_age, chronological_age), 3),
    '-> raw age recapitulates chronological age; report the residual instead\n')
cat('cor(EAA, chronological age):     ', round(cor(eaa, chronological_age), 3),
    '-> EAA is orthogonal to chronological age, as intended\n\n')

# Missing-CpG mean-imputation biases EAA toward zero: drop a third of the clock CpGs and
# impute them to the training (column) mean, then recompute.
training_mean <- colMeans(beta)
missing_frac <- 1 / 3
n_missing <- round(n_clock_cpgs * missing_frac)
missing_idx <- sample(seq_len(n_clock_cpgs), n_missing)

beta_imputed <- beta
beta_imputed[, missing_idx] <- matrix(training_mean[missing_idx], n_samples, n_missing, byrow = TRUE)

dnam_age_imputed <- apply_clock(beta_imputed, weights, intercept)
dnam_age_imputed <- dnam_age_imputed - mean(dnam_age_imputed) + mean(chronological_age)
eaa_imputed <- resid(lm(dnam_age_imputed ~ chronological_age))

cpg_coverage <- 1 - missing_frac      # always report the fraction of clock CpGs actually present
sd_ratio <- sd(eaa_imputed) / sd(eaa)

cat('clock CpG coverage:', round(cpg_coverage * 100), '%  (', n_clock_cpgs - n_missing, 'of', n_clock_cpgs, 'present )\n')
cat('SD(EAA) after mean-imputing missing CpGs / SD(EAA) full:', round(sd_ratio, 3),
    '-> imputation shrinks variance, biasing acceleration toward zero\n')
if (cpg_coverage < 0.8) {
  cat('VERDICT: coverage below the 0.8 floor (min.perc) -> flag/refuse these clock estimates\n')
} else {
  cat('VERDICT: coverage above the 0.8 floor -> clock estimates acceptable\n')
}

# The real path on an array beta matrix (CpGs in rows, samples in columns):
run_real_clocks <- function(beta_matrix, pheno) {
  library(methylclock)
  checkClocks(beta_matrix)                         # report missing clock CpGs BEFORE estimating
  ages <- DNAmAge(beta_matrix, clocks = c('Horvath', 'Hannum', 'Levine', 'skinHorvath'),
                  age = pheno$age, cell.count = FALSE, min.perc = 0.8)
  # ages$ageAcc2 is the EAA endpoint (residual of DNAm age on chronological age)

  library(dnaMethyAge)
  availableClock()                                 # confirm installed clock-name strings
  phenoage <- methyAge(beta_matrix, clock = 'LevineM2018', age_info = pheno, fit_method = 'Linear')

  library(DunedinPACE)
  pace <- PACEProjector(beta_matrix)               # returns the DunedinPACE pace values; a RATE (~1.0 = normal), never residualize like an age
  list(ages = ages, phenoage = phenoage, pace = pace)
}
