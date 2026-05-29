# Randomization, blocking, and the experimental unit
# Reference: designit 0.5+, lme4 1.1-35+, lmerTest 3.1+ | Verify API if version differs
#
# Demonstrates: identifying the experimental unit and aggregating observational units,
# restricted (block) randomization with run-order randomization, and the mixed-model
# random-effects structure that matches a nested / split-plot design.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

# ---------------------------------------------------------------------------
# 1. The experimental unit defines n: aggregate cells (observational units) to donors
# ---------------------------------------------------------------------------
set.seed(20260528)
n_donors_per_group <- 4
cells_per_donor <- 200

cells <- expand.grid(cell = seq_len(cells_per_donor),
                     donor = paste0('D', seq_len(2 * n_donors_per_group)))
cells$condition <- ifelse(as.integer(sub('D', '', cells$donor)) <= n_donors_per_group,
                         'ctrl', 'treat')
# donor-level mean shift + cell-level noise
donor_effect <- setNames(rnorm(length(unique(cells$donor)), 0, 0.5), unique(cells$donor))
cells$measurement <- donor_effect[cells$donor] +
  ifelse(cells$condition == 'treat', 0.3, 0) + rnorm(nrow(cells), 0, 1)

# Correct unit of inference: one value per donor (the experimental unit), not per cell.
eu_level <- cells |>
  group_by(donor, condition) |>
  summarise(value = mean(measurement), .groups = 'drop')
cat('Experimental units (n) per group:',
    paste(table(eu_level$condition), collapse = ' / '), '\n')   # 4 / 4, NOT 800 / 800

# Test on experimental-unit values (n = 8 donors), not on 1600 pseudoreplicated cells.
t.test(value ~ condition, data = eu_level)

# ---------------------------------------------------------------------------
# 2. Restricted (block) randomization + run-order randomization
# ---------------------------------------------------------------------------
# 24 samples processed 8/day over 3 days; randomize condition WITHIN day (block),
# and randomize the processing order so position is not confounded with condition.
units <- data.frame(id = sprintf('S%02d', 1:24),
                    day = rep(c('day1', 'day2', 'day3'), each = 8))
units$condition <- ave(units$id, units$day, FUN = function(ids)
  sample(rep(c('ctrl', 'treat'), length.out = length(ids))))
units$run_order <- sample(nrow(units))
# Confirm balance: each day holds equal ctrl/treat -> day is orthogonal to condition
print(table(units$day, units$condition))

# ---------------------------------------------------------------------------
# 3. Mixed model matching a nested design (cells within donor)
# ---------------------------------------------------------------------------
# Random intercept for donor accounts for the within-donor correlation that makes
# cells pseudoreplicates; condition stays a fixed effect.
fit_nested <- lmer(measurement ~ condition + (1 | donor), data = cells)
print(anova(fit_nested))   # Satterthwaite df via lmerTest; effective n driven by donors

# Split-plot pattern (whole plot = processing day/run, sub-plot = condition):
#   lmer(response ~ condition + (1 | day/sample), data = df)
# The whole-plot factor must be tested against whole-plot error, never sub-plot error.
