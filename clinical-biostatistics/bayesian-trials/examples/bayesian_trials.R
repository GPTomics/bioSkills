# Bayesian clinical trial design examples.
#
# Reference: RBesT 1.7+ | Verify API if version differs
# Reference: OncoBayes2 0.8+ | Verify API if version differs
# Reference: BOIN 2.7+ | Verify API if version differs
#
# Covers BOIN/CRM Phase 1, MAP + robust MAP, EXNEX basket, power prior,
# Berry-Berry hierarchical AE, and posterior probability stopping.


# ----------------------------------------------------------------------
# 1. BOIN Phase 1 dose-finding (FDA Fit-for-Purpose Dec 2021)
# ----------------------------------------------------------------------
library(BOIN)

# Generate the escalation table for the protocol
boundary_table <- get.boundary(
    target = 0.30,
    ncohort = 10,
    cohortsize = 3,
    n.earlystop = 12,
    p.saf = 0.6 * 0.30,
    p.tox = 1.4 * 0.30
)
print(boundary_table)
# Pre-print at investigator desk; transparent decision rule

# Operating characteristics over true DLT scenarios
oc <- get.oc(
    target = 0.30,
    p.true = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),
    ncohort = 10,
    cohortsize = 3,
    ntrial = 1000
)
print(oc)
# Reports: MTD selection accuracy, overdose risk, sample size distribution


# ----------------------------------------------------------------------
# 2. CRM with Lee-Cheung 2009 calibrated skeleton
# ----------------------------------------------------------------------
library(dfcrm)

# Calibrate prior skeleton via indifference-interval method
target_dlt <- 0.30
nlevel <- 6
prior_skeleton <- getprior(
    halfwidth = 0.05,
    target = target_dlt,
    nu = 3,                  # prior MTD level
    nlevel = nlevel
)
print(prior_skeleton)

# Simulate CRM OCs
crm_oc <- crmsim(
    PI = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),
    prior = prior_skeleton,
    target = target_dlt,
    n = 30,
    x0 = 1,
    nsim = 1000,
    method = 'bayes',
    model = 'logistic'
)
print(crm_oc)


# ----------------------------------------------------------------------
# 3. MAP prior and robust MAP via RBesT
# ----------------------------------------------------------------------
library(RBesT)

# Historical control data
historical <- data.frame(
    study = c('h1', 'h2', 'h3', 'h4'),
    n = c(40, 35, 50, 45),
    r = c(8, 6, 12, 9)
)

# Fit MAP via gMAP (Stan-based random-effects meta-analysis)
map_fit <- gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = historical,
    family = binomial,
    tau.dist = 'HalfNormal',
    tau.prior = 0.5,           # between-study SD prior
    beta.prior = cbind(0, 2)
)
print(map_fit)

# Approximate posterior with mixture for analytical use
map_mix <- automixfit(map_fit, Nc = 2)
print(map_mix)

# Effective sample size of the prior (Schmidli 2014)
ess(map_mix, method = 'morita')
ess(map_mix, method = 'elir')

# Robust MAP: add vague mixture component to guard against prior-data conflict
robust_map <- robustify(map_mix, weight = 0.2, mean = 0.5, n = 1)
print(robust_map)
ess(robust_map)


# ----------------------------------------------------------------------
# 4. Sensitivity over tau prior
# ----------------------------------------------------------------------
# Tau (between-study SD) drives borrowing strength
# Sensitivity over tau_prior = 0.1, 0.25, 0.5, 1.0
for (tau_prior in c(0.1, 0.25, 0.5, 1.0)) {
    map_sens <- gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = historical,
        family = binomial,
        tau.dist = 'HalfNormal',
        tau.prior = tau_prior,
        beta.prior = cbind(0, 2)
    )
    ess_sens <- ess(automixfit(map_sens, Nc = 2))
    cat(sprintf('tau_prior = %.2f: ESS = %.1f\n', tau_prior, ess_sens))
}


# ----------------------------------------------------------------------
# 5. EXNEX for basket trial across rare-disease strata
# ----------------------------------------------------------------------
# OncoBayes2 has the canonical EXNEX implementation for combination dose-finding
# For basket trials, bhmbasket is the simpler alternative

# library(bhmbasket)
# # Pseudocode -- consult bhmbasket documentation for current API
# basket_design <- create_basket_design(n_baskets = 5,
#                                          n_per_basket = c(20, 25, 15, 30, 20),
#                                          target_response = 0.40,
#                                          null_response = 0.20)
# # Default EXNEX weights 0.5 EX / 0.5 NEX (Neuenschwander 2016)


# ----------------------------------------------------------------------
# 6. Power prior for pediatric extrapolation
# ----------------------------------------------------------------------
# library(bayesDP)
# # Adult historical data: n_adult = 500, responders_adult = 200
# # Pediatric current: n_ped = 40, responders_ped = 20
# # Apply discount alpha (power) in [0, 1]; FDA Jan 2026 draft: 0.3-0.6 typical
# pp_fit <- bdpbinomial(
#     y_t = 20, N_t = 40,           # current pediatric treatment
#     y0_t = 200, N0_t = 500,        # historical adult treatment
#     y_c = 15, N_c = 40,           # current pediatric control
#     y0_c = 180, N0_c = 500,        # historical adult control
#     discount_function = 'identity',
#     alpha_max = 1.0,
#     fix_alpha = TRUE,
#     alpha_value = 0.50             # 50% borrowing
# )
# print(pp_fit)


# ----------------------------------------------------------------------
# 7. Berry-Berry 3-level hierarchical model for AE multiplicity
# ----------------------------------------------------------------------
# library(c212)
# # Model: AE within MedDRA PT within SOC
# # Spike-and-slab on log OR
# # Borrows within SOC; shrinks toward 0 if no evidence

# Conceptual data format:
# data.frame(USUBJID, treatment, AE_PT, AE_SOC)
# Berry-Berry returns posterior probability of treatment effect for each PT
# JMP Clinical also implements this for industry use


# ----------------------------------------------------------------------
# 8. Posterior probability stopping (Berry et al 2010)
# ----------------------------------------------------------------------
# Stop for efficacy at interim if P(theta > 0 | data) > p_efficacy_threshold
# Stop for futility if P(theta > delta_clin | data) < p_futility_threshold
# CALIBRATE thresholds via simulation under null to control frequentist Type-I

# Conceptual simulation loop:
n_sim <- 5000
type_I_count <- 0
for (sim in 1:n_sim) {
    # Simulate trial data under null (theta = 0)
    # Apply Bayesian stopping rule at interim looks
    # If trial stops for efficacy or final p_post > threshold, count as Type-I event
}
type_I_rate <- type_I_count / n_sim
# Adjust threshold until type_I_rate = nominal alpha (e.g., 0.025)


# ----------------------------------------------------------------------
# 9. I-SPY 2 style graduation criterion (conceptual)
# ----------------------------------------------------------------------
# Graduation: PP(success in 300-patient Phase 3 trial | current data) >= 0.85

# Implementation:
# 1. Fit hierarchical model to current platform data:
#    response ~ arm + biomarker_subtype + arm:subtype
# 2. Draw posterior treatment effects by subtype
# 3. For each draw, simulate Phase 3 trial of size 300:
#    treatment vs control, compute test statistic
# 4. Compute proportion of draws meeting Phase 3 success criterion (e.g., p < 0.025)
# 5. If proportion >= 0.85, arm graduates

# Pseudocode (requires custom Stan or FACTS commercial software):
# graduation_pp <- function(model_fit, target_subtype, phase3_n = 300, n_draws = 5000) {
#     posterior_effects <- as.matrix(model_fit)[, paste0('treatment_effect_', target_subtype)]
#     phase3_successes <- 0
#     for (i in 1:n_draws) {
#         theta <- posterior_effects[i]
#         # Simulate Phase 3 trial
#         control_responses <- rbinom(phase3_n/2, 1, baseline_rate)
#         treatment_responses <- rbinom(phase3_n/2, 1, baseline_rate + theta)
#         # Compute Phase 3 test
#         p_value <- prop.test(c(sum(treatment_responses), sum(control_responses)),
#                              c(phase3_n/2, phase3_n/2))$p.value
#         if (p_value < 0.025) phase3_successes <- phase3_successes + 1
#     }
#     return(phase3_successes / n_draws)
# }


# ----------------------------------------------------------------------
# 10. Stan convergence diagnostics for regulatory submissions
# ----------------------------------------------------------------------
# library(rstan)
# library(bayesplot)

# After fitting any Stan model:
# print(fit, pars = c('theta', 'tau'))
# # Check: R-hat < 1.01 for all parameters
# # Check: ESS (n_eff) > 1000-2000 per chain per parameter
# # Trace plots: mcmc_trace(fit)
# # Posterior predictive: ppc_dens_overlay(y_observed, y_rep)

# For FDA submission, document:
# - Number of chains, iterations, warmup
# - Initial values strategy
# - Seed
# - R/Stan version
# - Convergence diagnostics in submission appendix
