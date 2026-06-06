# Reference: ropls 1.34+ | Verify API if version differs
# Permutation-validated OPLS-DA on synthetic metabolomics data, demonstrating that
# (1) only pQ2 -- not R2Y or the score plot -- licenses a discriminant claim, and
# (2) the scaling choice (Pareto vs unit-variance) is a hypothesis that changes the VIP list.
library(ropls)

set.seed(1)
n_per_group <- 20
n_features <- 300
n_true <- 15                                  # only the first 15 features actually differ

# Synthetic data is already homoscedastic and roughly normal; real MS intensities are
# right-skewed and heteroscedastic, so log/glog-transform them BEFORE opls() in practice.
group <- factor(rep(c('control', 'case'), each = n_per_group))
intensities <- matrix(rnorm(2 * n_per_group * n_features, mean = 10, sd = 1),
                      nrow = 2 * n_per_group, ncol = n_features)
colnames(intensities) <- paste0('M', seq_len(n_features))
case_rows <- group == 'case'
intensities[case_rows, seq_len(n_true)] <- intensities[case_rows, seq_len(n_true)] + 1.2

# A NULL model: same data, labels permuted. OPLS-DA will still separate it in p>>n,
# so a clean score plot here is the geometry, not signal -- pQ2 is the honest check.
null_group <- factor(sample(as.character(group)))

run_model <- function(y, scaleC) {
    opls(intensities, y, predI = 1, orthoI = NA, scaleC = scaleC,
         permI = 1000, crossvalI = 7, fig.pdfC = 'none', info.txtC = 'none')
}

cat('=== Real labels, Pareto scaling ===\n')
real_pareto <- run_model(group, 'pareto')
print(getSummaryDF(real_pareto)[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

cat('\n=== Permuted (null) labels, Pareto scaling -- expect high R2Y, failed pQ2 ===\n')
null_pareto <- run_model(null_group, 'pareto')
print(getSummaryDF(null_pareto)[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

cat('\n=== Real labels, unit-variance scaling -- compare the VIP ranking ===\n')
real_uv <- run_model(group, 'standard')
print(getSummaryDF(real_uv)[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

vip_pareto <- getVipVn(real_pareto)
vip_uv <- getVipVn(real_uv)
top10_pareto <- names(sort(vip_pareto, decreasing = TRUE))[1:10]
top10_uv <- names(sort(vip_uv, decreasing = TRUE))[1:10]
overlap <- length(intersect(top10_pareto, top10_uv))
cat('\nTop-10 VIP overlap between Pareto and UV scaling:', overlap, 'of 10\n')
cat('Scaling-fragile result if the two lists diverge.\n')

# pQ2 (the permutation p-value for Q2) is the licensing gate -- it, not R2Y or the score
# plot, exposes a model built on noise. Q2's magnitude is the effect size and should also
# clear the Triba >0.5 heuristic to be worth reporting.
summ <- getSummaryDF(real_pareto)
licensed <- summ[['pQ2']] < 0.05
cat(sprintf('\nReal model: Q2(cum)=%.2f, pQ2=%.3f -> permutation-licensed: %s\n',
            summ[['Q2(cum)']], summ[['pQ2']], licensed))
