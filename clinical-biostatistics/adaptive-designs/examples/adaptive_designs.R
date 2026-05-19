# Adaptive clinical trial design examples.
#
# Reference: rpact 4.2+ | Verify API if version differs
# Reference: gsDesign 3.6+ | Verify API if version differs
#
# Covers group-sequential, blinded/unblinded SSR, Mehta-Pocock promising zone,
# BOIN phase 1 dose-finding, and adaptive enrichment.

# ----------------------------------------------------------------------
# 1. O'Brien-Fleming group-sequential design for survival
# ----------------------------------------------------------------------
library(gsDesign)

design_obf <- gsDesign(
    k = 3,                  # 2 interim + 1 final
    test.type = 1,          # 1-sided efficacy
    alpha = 0.025,
    beta = 0.10,            # power = 0.90
    sfu = sfLDOF,           # Lan-DeMets approximation of OBF
    timing = c(0.33, 0.67, 1.0)
)
print(design_obf)
# Boundaries: very conservative early (~0.0001 nominal at 33% info),
# near-nominal at final (~0.0234 of 0.025)

plot(design_obf)


# ----------------------------------------------------------------------
# 2. Sample size for survival group-sequential
# ----------------------------------------------------------------------
n_gs <- gsSurv(
    k = 3,
    test.type = 2,           # 2-sided
    alpha = 0.025,           # one-sided alpha for each side
    beta = 0.10,             # 90% power
    sfu = sfLDOF,
    lambdaC = 0.04,          # control hazard per month
    hr = 0.65,               # treatment HR
    eta = 0.005,             # dropout hazard per month
    T = 30,                  # total study duration in months
    minfup = 18,             # minimum follow-up
    ratio = 1                # 1:1 allocation
)
print(n_gs)


# ----------------------------------------------------------------------
# 3. Blinded SSR for continuous endpoint (Friede-Kieser 2006)
# ----------------------------------------------------------------------
library(rpact)

# Pre-specify design with uncertain variance
design_initial <- getDesignGroupSequential(
    kMax = 1,                # final-analysis only initially
    alpha = 0.025,
    beta = 0.20,             # 80% power
    sided = 1,
    typeOfDesign = 'asUser'
)

# Initial sample size assuming SD = 12
ss_initial <- getSampleSizeMeans(
    design = design_initial,
    alternative = 5,         # detect mean diff of 5
    stDev = 12,
    groups = 2
)
print(ss_initial)

# At internal pilot, re-estimate SD from blinded data
# (manual: pool all observations, compute SD ignoring arm)
# If observed SD = 14 (higher than assumed), recompute n
ss_recomputed <- getSampleSizeMeans(
    design = design_initial,
    alternative = 5,
    stDev = 14,
    groups = 2
)
print(ss_recomputed)
# Type-I error NOT affected because test ignores SSR


# ----------------------------------------------------------------------
# 4. Mehta-Pocock promising zone SSR with CHW weights
# ----------------------------------------------------------------------
# Inverse normal combination test with adaptive sample size
design_promising <- getDesignInverseNormal(
    kMax = 2,
    alpha = 0.025,
    beta = 0.20,
    sided = 1,
    informationRates = c(0.5, 1.0),
    typeOfDesign = 'asUser',
    userAlphaSpending = c(0.0125, 0.025)  # spending function
)

# Compute initial sample size
ss_promising <- getSampleSizeMeans(
    design = design_promising,
    alternative = 5,
    stDev = 12,
    groups = 2
)
print(ss_promising)

# At interim, if conditional power in (0.3, 0.8), increase to n_max
# Use rpact's getConditionalPower() and getSampleSizeMeans() for re-estimation
# Pre-specify CHW weights from original design


# ----------------------------------------------------------------------
# 5. BOIN Phase 1 dose-finding (FDA Fit-for-Purpose 2021)
# ----------------------------------------------------------------------
library(BOIN)

# Generate pre-tabulated escalation decisions
boin_table <- get.boundary(
    target = 0.30,           # target DLT rate
    ncohort = 10,            # 10 cohorts of size 3 -> max 30 patients
    cohortsize = 3,
    n.earlystop = 12,        # stop early at lowest dose if 12 patients show futility/safety
    p.saf = 0.6 * 0.30,      # boundary for "safe" (escalate)
    p.tox = 1.4 * 0.30       # boundary for "toxic" (de-escalate)
)
print(boin_table)
# This table is printed in the protocol; investigator looks up at bedside
# No real-time Bayesian software needed (Why FDA prefers BOIN)

# Simulate operating characteristics
boin_oc <- get.oc(
    target = 0.30,
    p.true = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),  # true DLT rates per dose
    ncohort = 10,
    cohortsize = 3,
    ntrial = 1000,
    n.earlystop = 12
)
print(boin_oc)
# Reports correct MTD selection rate, overdose risk, average sample size


# ----------------------------------------------------------------------
# 6. Continual Reassessment Method (CRM) -- alternative to BOIN
# ----------------------------------------------------------------------
library(dfcrm)

# CRM with logistic skeleton
prior_skeleton <- c(0.05, 0.10, 0.20, 0.30, 0.50)  # prior DLT rates per dose
target <- 0.30

# Simulate CRM operating characteristics
crm_sim <- crmsim(
    PI = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),  # true DLT rates
    prior = prior_skeleton,
    target = target,
    n = 30,
    x0 = 1,                  # starting dose
    nsim = 1000,
    mcohort = 1,             # cohort size
    method = 'bayes',
    model = 'logistic'
)
print(crm_sim)
# Compare to BOIN OCs; CRM more efficient under correct skeleton, sensitive to skeleton mis-spec


# ----------------------------------------------------------------------
# 7. Adaptive enrichment (population selection)
# ----------------------------------------------------------------------
# rpact has limited native enrichment support; use adaptr or custom implementation

# Conceptual: at interim, evaluate conditional power in:
# - Full population
# - Biomarker-positive subgroup
# If full CP < threshold and subgroup CP > threshold, drop biomarker-negative
# Closed-testing across full and subgroup preserves FWER

# library(adaptr)
# adapt_enrich_design <- create_adaptr_design(...)


# ----------------------------------------------------------------------
# 8. EXNEX for basket trial across rare-disease strata
# ----------------------------------------------------------------------
library(RBesT)

# Mixture of exchangeable + non-exchangeable per stratum
# Default weights 0.5 EX / 0.5 NEX (Neuenschwander 2016)
# Stratified borrowing -- one stratum can detach if truly different

# Simplified MAP prior via RBesT
historical_data <- data.frame(
    study = c('study1', 'study2', 'study3'),
    n = c(40, 35, 50),
    r = c(8, 6, 12)
)

map_prior <- gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = historical_data,
    family = binomial,
    tau.dist = 'HalfNormal',
    tau.prior = 0.5,
    beta.prior = cbind(0, 2)
)
print(map_prior)
# Effective sample size of MAP prior (Schmidli 2014)
print(ess(map_prior))


# ----------------------------------------------------------------------
# 9. I-SPY 2 style graduation criterion (conceptual)
# ----------------------------------------------------------------------
# Posterior predictive probability of success in Phase 3
# graduation = PP(Phase 3 trial of size 300 succeeds | current data) >= 0.85
# Implementation requires Stan/JAGS or FACTS (Berry Consultants commercial)

# Conceptual scaffold (requires custom Bayesian implementation):
# 1. Fit hierarchical model to current platform data
# 2. Simulate forward: draw treatment effect from posterior
# 3. For each draw, simulate Phase 3 trial of size 300
# 4. Compute proportion of simulations achieving Phase 3 success
# 5. If proportion >= 0.85, arm graduates


# ----------------------------------------------------------------------
# 10. RAR with time-trend adjustment (Robertson 2023 consensus)
# ----------------------------------------------------------------------
# RAR appropriate for multi-arm (>=3 arms) and rare disease
# Pre-specify time-trend covariate (e.g., enrollment quarter) in primary analysis
# Use proper analysis weights (e.g., inverse-probability weighting)

# Conceptual: in 4-arm trial, update allocation probabilities based on posterior
# P(superior | data) every quarter; floor at minimum allocation (e.g., 10% per arm)
# Primary analysis includes enrollment quarter as covariate in Cox or logistic
