# Batch design: balanced assignment, hidden-batch detection, correction caveats
# Reference: sva 3.50+, designit 0.5+ | Verify API if version differs
#
# Demonstrates: the confounded-vs-balanced contrast, constrained sample-to-batch
# assignment, surrogate variable estimation for hidden batches, and the rule that a
# batch-"cleaned" matrix is for visualization only -- inference keeps batch in the model.

suppressPackageStartupMessages(library(sva))
set.seed(20260528)

# ---------------------------------------------------------------------------
# 1. Confounded vs balanced design
# ---------------------------------------------------------------------------
# BAD: condition aliased with batch -> non-identifiable, no correction can rescue it.
bad <- data.frame(condition = rep(c('ctrl', 'treat'), each = 12),
                  batch = rep(c('B1', 'B2'), each = 12))
cat('Confounded design (batch aliased with condition):\n'); print(table(bad$condition, bad$batch))

# GOOD: balance condition (and sex) across batches -> batch orthogonal, estimable.
samples <- data.frame(id = sprintf('S%02d', 1:24),
                     condition = rep(c('ctrl', 'treat'), each = 12),
                     sex = rep(c('M', 'F'), 12))
samples$batch <- NA_character_
for (cond in unique(samples$condition)) {
  idx <- which(samples$condition == cond)
  samples$batch[idx] <- sample(rep(paste0('B', 1:3), length.out = length(idx)))
}
cat('\nBalanced design (condition orthogonal to batch):\n'); print(table(samples$condition, samples$batch))

# Constrained assignment with designit (verify API vs installed vignette):
#   library(designit)
#   bc <- BatchContainer$new(dimensions = list(batch = 3, position = 8))
#   bc <- assign_in_order(bc, samples = samples)
#   bc <- optimize_design(bc, scoring = osat_score_generator(
#           batch_vars = 'batch', feature_vars = c('condition', 'sex')))

# ---------------------------------------------------------------------------
# 2. Simulate a batch effect and detect hidden structure with SVA
# ---------------------------------------------------------------------------
n_genes <- 1000; n <- nrow(samples)
counts <- matrix(rnbinom(n_genes * n, mu = 100, size = 10), nrow = n_genes,
                 dimnames = list(paste0('Gene', 1:n_genes), samples$id))
batch_mult <- c(B1 = 1.0, B2 = 1.5, B3 = 0.7)              # multiplicative batch effect
counts <- sweep(counts, 2, batch_mult[samples$batch], '*')
counts[1:50, samples$condition == 'treat'] <- counts[1:50, samples$condition == 'treat'] * 2  # true DE
expr <- log2(counts + 1)

mod  <- model.matrix(~ condition, data = samples)
mod0 <- model.matrix(~ 1, data = samples)
n_sv <- num.sv(expr, mod, method = 'leek')
cat('\nEstimated hidden factors (surrogate variables):', n_sv, '\n')
svobj <- sva(expr, mod, mod0, n.sv = max(n_sv, 1))
# Correct use: add svobj$sv as covariates to the DE model (in differential-expression),
# NOT subtract them from `expr` before testing.

# ---------------------------------------------------------------------------
# 3. Cleaned matrix for VISUALIZATION ONLY (never feed into the hypothesis test)
# ---------------------------------------------------------------------------
# ComBat on log-normalized data; for integer counts use sva::ComBat_seq() instead.
expr_viz <- ComBat(dat = expr, batch = samples$batch, mod = mod)
pc_before <- cor(prcomp(t(expr))$x[, 1], as.numeric(factor(samples$batch)))
pc_after  <- cor(prcomp(t(expr_viz))$x[, 1], as.numeric(factor(samples$batch)))
cat(sprintf('PC1-vs-batch correlation: %.2f (raw) -> %.2f (ComBat, for plots only)\n',
            pc_before, pc_after))
# Inference rule (Nygaard 2016): keep batch in the model -- e.g. ~ condition + batch (+ SVs) --
# rather than testing on the ComBat-cleaned matrix, which understates residual variance.
