# Reference: cause 1.2.0+, TwoSampleMR 0.5.11+ | Verify API if version differs
## CAUSE (Morrison 2020) for CHP-aware causal estimation
##
## CAUSE explicitly models correlated horizontal pleiotropy (CHP) via a shared-factor
## mixture component, distinguishing causal effect from shared genetic architecture.
## Required when LDSC rg(exposure, outcome) >= 0.3 or biology suggests a shared
## upstream factor. UHP-only methods (IVW, Egger, PRESSO) are blind to CHP.
##
## Prerequisites:
##   remotes::install_github('jean997/cause')
##   GWAS sumstats with matched effect-allele coding for both exposure and outcome
##   At least ~100 genome-wide significant SNPs (p < 5e-8) after LD pruning

library(cause)
library(dplyr)

set.seed(2024)
n_genome <- 50000

snps <- paste0('rs', 1:n_genome)
beta_hat_1 <- rnorm(n_genome, 0, 0.005)
seb1 <- abs(rnorm(n_genome, 0.01, 0.002))

n_sig <- 150
sig_idx <- sample(1:n_genome, n_sig)
beta_hat_1[sig_idx] <- rnorm(n_sig, 0, 0.05)
seb1[sig_idx] <- 0.008

true_gamma <- 0.3
true_eta <- 0.4
shared_factor_loadings <- rep(0, n_genome)
shared_factor_loadings[sig_idx[1:50]] <- rnorm(50, 0, 0.04)

beta_hat_2 <- true_gamma * beta_hat_1 + true_eta * shared_factor_loadings + rnorm(n_genome, 0, 0.005)
seb2 <- abs(rnorm(n_genome, 0.01, 0.002))

X <- new_cause_data(data.frame(
    snp = snps,
    beta_hat_1 = beta_hat_1,
    seb1 = seb1,
    beta_hat_2 = beta_hat_2,
    seb2 = seb2,
    A1 = 'A', A2 = 'G',
    stringsAsFactors = FALSE))

cat('=== Nuisance parameter estimation ===\n')
nuisance_snps <- sample(snps, min(1000000, n_genome))
params <- est_cause_params(X, nuisance_snps)
cat('rho (sample-overlap correction):', round(params$rho, 4), '\n')

cat('\n=== LD pruning to signature SNPs ===\n')
sig_filter <- X %>% filter(2 * pnorm(-abs(beta_hat_1 / seb1)) < 5e-8)
pruned_snps <- sig_filter$snp
cat('Sig SNPs after pruning:', length(pruned_snps), '\n')

if (length(pruned_snps) < 100) {
    cat('WARNING: <100 sig SNPs; CAUSE delta_ELPD CIs will be uninformative\n')
    cat('         Consider LHC-MR (genome-wide) or LCV (gcp) instead\n')
}

cat('\n=== CAUSE model fit ===\n')
res_cause <- cause(X = X, variants = pruned_snps, param_ests = params)

cat('\n=== Posterior summary ===\n')
summary_tab <- summary(res_cause)
print(summary_tab$tab)

cat('\n=== Model comparison (sharing vs causal) ===\n')
# res_cause$elpd is the ELPD comparison data.frame from `loo::loo_compare`-style output.
# Columns: model1, model2, delta_elpd, se_delta_elpd, z. Row 3 = sharing vs causal.
elpd_tab <- res_cause$elpd
print(elpd_tab)
row_sharing_vs_causal <- which(elpd_tab$model1 == 'sharing' & elpd_tab$model2 == 'causal')
delta_elpd <- elpd_tab$delta_elpd[row_sharing_vs_causal]
se_delta   <- elpd_tab$se_delta_elpd[row_sharing_vs_causal]
z_cause    <- elpd_tab$z[row_sharing_vs_causal]
p_one_sided <- pnorm(z_cause)
cat('delta_ELPD (sharing - causal; negative -> causal preferred):', round(delta_elpd, 3), '\n')
cat('SE:', round(se_delta, 3), '\n')
cat('z:', round(z_cause, 3), ' | one-sided p:', format.pval(p_one_sided), '\n')

if (p_one_sided < 0.05) {
    cat('Interpretation: Causal model preferred over CHP-only sharing model\n')
} else {
    cat('Interpretation: Cannot distinguish causal from shared-factor explanation\n')
    cat('               IVW / Egger / PRESSO estimates may be biased by CHP\n')
}

cat('\n=== Causal estimate (CAUSE causal model) ===\n')
# summary_tab$quants is a list with one entry per fitted model (sharing, causal).
# Each entry is a matrix of quantiles (rows = parameter, columns = lower/median/upper).
causal_quants <- summary_tab$quants[[2]]
gamma_est <- causal_quants[, 'gamma']
cat('gamma (causal effect) posterior:\n')
cat('  Lower 95%:', round(gamma_est[1], 4), '\n')
cat('  Median  :', round(gamma_est[2], 4), '\n')
cat('  Upper 95%:', round(gamma_est[3], 4), '\n')
cat('  True    :', true_gamma, '\n')

cat('\n=== Sharing parameters (q = CHP fraction; eta = shared-factor strength) ===\n')
q_est <- causal_quants[, 'q']
eta_est <- causal_quants[, 'eta']
cat('q (CHP fraction):\n')
cat('  Lower 95%:', round(q_est[1], 4), '\n')
cat('  Median  :', round(q_est[2], 4), '\n')
cat('  Upper 95%:', round(q_est[3], 4), '\n')
cat('eta (shared-factor effect):\n')
cat('  Lower 95%:', round(eta_est[1], 4), '\n')
cat('  Median  :', round(eta_est[2], 4), '\n')
cat('  Upper 95%:', round(eta_est[3], 4), '\n')
cat('  True    :', true_eta, '\n')

cat('\n=== Diagnostic: Pareto-k ===\n')
# Pareto-k diagnostics come from the loo objects in res_cause$loos[[2]] (sharing)
# and res_cause$loos[[3]] (causal); pareto_k_table() is loo's helper.
loo_causal <- res_cause$loos[[3]]
if (!is.null(loo_causal)) {
    pk <- loo::pareto_k_values(loo_causal)
    n_high_k <- sum(pk > 0.7, na.rm = TRUE)
    cat('Pareto-k > 0.7 (unstable):', n_high_k, '/', length(pk), '\n')
    if (n_high_k > 0.1 * length(pk)) {
        cat('WARNING: substantial unstable points; treat posterior with caution\n')
    }
}

cat('\n=== Reporting (STROBE-MR Item 17: CHP-aware sensitivity) ===\n')
cat('Method: CAUSE (Morrison 2020 Nat Genet 52:740)\n')
cat('Signature SNPs:', length(pruned_snps), '\n')
cat('rho (sample overlap):', round(params$rho, 4), '\n')
cat('gamma (causal effect): ', round(gamma_est[2], 4),
    ' (95% CI ', round(gamma_est[1], 4), ', ', round(gamma_est[3], 4), ')\n', sep = '')
cat('q (CHP fraction): ', round(q_est[2], 4),
    ' (95% CI ', round(q_est[1], 4), ', ', round(q_est[3], 4), ')\n', sep = '')
cat('delta_ELPD (sharing - causal): ', round(delta_elpd, 3), ' (one-sided p ',
    format.pval(p_one_sided), ')\n', sep = '')
