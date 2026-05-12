# Reference: mediation 4.5+, CMAverse 0.1.0+, EValue 4.1+ | Verify API if version differs
## 4-way mediation decomposition (VanderWeele 2014) with exposure-mediator interaction
##
## Decomposes total effect into:
##   CDE    = controlled direct effect
##   PIE    = pure indirect effect (mediation without interaction)
##   INTref = interaction-reference effect (interaction without mediation)
##   INTmed = mediated interaction (interaction AND mediation)
## Total = CDE + INTref + INTmed + PIE
##
## CMAverse install:
##   remotes::install_github('BS1125/CMAverse')

library(CMAverse)
library(mediation)
library(EValue)

set.seed(42)
n <- 1500

## Simulate exposure, mediator with exposure-dependent variance,
## and binary outcome with exposure-mediator interaction
genotype <- rbinom(n, 2, 0.3)
age <- scale(rnorm(n, 55, 10))[, 1]
sex <- rbinom(n, 1, 0.5)
pc1 <- rnorm(n)
pc2 <- rnorm(n)

## Mediator depends on exposure (alpha = 0.5)
expression <- 0.5 * genotype + 0.05 * age - 0.1 * sex + 0.2 * pc1 + rnorm(n, 0, 0.8)

## Outcome: exposure direct effect 0.15, mediator effect 0.3, interaction 0.25
## Rare-ish disease (~12% prevalence) so OR-based 4-way is approximately valid
logit_p <- -2.5 + 0.15 * genotype + 0.3 * expression +
           0.25 * genotype * expression + 0.04 * age + 0.1 * sex + 0.1 * pc1
disease <- rbinom(n, 1, plogis(logit_p))

dat <- data.frame(genotype, expression, disease, age, sex, pc1, pc2)
cat('Outcome prevalence:', round(mean(dat$disease), 3), '\n\n')

## --- CMAverse regression-based 4-way decomposition ---
## EMint=TRUE adds the exposure-mediator interaction term to the outcome model
## astar=0 = reference exposure level, a=1 = "increment by 1 allele" contrast
## yreg='logistic' for binary outcome; mreg=list('linear') for continuous mediator
## inference='bootstrap' with nboot=1000 (use 5000 for publication)
result_4way <- cmest(
  data=dat, model='rb',
  outcome='disease', exposure='genotype', mediator='expression',
  basec=c('age','sex','pc1','pc2'),
  EMint=TRUE,
  mreg=list('linear'), yreg='logistic',
  astar=0, a=1, mval=list(0),
  estimation='paramfunc', inference='bootstrap', nboot=1000,
  boot.ci.type='per'
)

cat('--- 4-Way Decomposition ---\n')
summary(result_4way)

## Component meaning under rare-disease logistic 4-way decomposition (Valeri 2013).
## CMAverse output column names (case-sensitive; verify with summary(result_4way)$results):
##   Rcde   = controlled direct effect on OR scale (mediator fixed at mval)
##   Rpnie  = pure natural indirect effect (mediator without interaction)
##   intref = reference interaction term (interaction WITHOUT mediation)
##   intmed = mediated interaction term (interaction AND mediation)
##   Rte    = total OR
##   pm     = proportion mediated; int, pe = proportion attributable to interaction
## Component names differ between continuous outcome (cde/pnde/tnde/pnie/tnie) and
## ratio scale (Rcde/Rpnde/Rtnde/Rpnie/Rtnie); inspect names(summary(result)$results).

## --- Compare to standard (no-interaction) mediation ---
med_simple <- lm(expression ~ genotype + age + sex + pc1 + pc2, data=dat)
out_simple <- glm(disease ~ genotype + expression + age + sex + pc1 + pc2,
                  data=dat, family=binomial)

med_imai <- mediate(med_simple, out_simple, treat='genotype', mediator='expression',
                    boot=TRUE, sims=1000)
cat('\n--- Standard counterfactual mediation (no interaction term) ---\n')
cat('  ACME:', round(med_imai$d0, 4),
    ' [', round(med_imai$d0.ci[1], 4), ',', round(med_imai$d0.ci[2], 4), ']\n')
cat('  ADE :', round(med_imai$z0, 4),
    ' [', round(med_imai$z0.ci[1], 4), ',', round(med_imai$z0.ci[2], 4), ']\n')
cat('  Total:', round(med_imai$tau.coef, 4), '\n')
cat('  PM   :', round(med_imai$n0, 3), '\n')

## Without the interaction term, ACME and ADE silently absorb INTref + INTmed
## into the wrong components. 4-way is the diagnostic.

## --- Mediational E-value (Smith & VanderWeele 2019) ---
## Risk-ratio-scale strength an unmeasured M-Y confounder would need to nullify ACME
acme_rr <- exp(med_imai$d0)
acme_lower_rr <- exp(med_imai$d0.ci[1])
acme_upper_rr <- exp(med_imai$d0.ci[2])
ev <- evalues.RR(est=acme_rr, lo=acme_lower_rr, hi=acme_upper_rr)
cat('\n--- Mediational E-value ---\n')
print(ev)
## E-value > 2.0 suggests robustness; < 1.5 suggests strong sensitivity to confounding

## --- Imai rho-based sensitivity (probit-link refit required for medsens) ---
out_probit <- glm(disease ~ genotype + expression + age + sex + pc1 + pc2,
                  data=dat, family=binomial(link='probit'))
med_probit <- mediate(med_simple, out_probit, treat='genotype', mediator='expression',
                      boot=TRUE, sims=500)
sens <- medsens(med_probit, rho.by=0.05, effect.type='indirect', sims=500)
cat('\n--- Imai rho sensitivity ---\n')
summary(sens)
## rho_crit (critical rho where ACME crosses 0):
##   |rho_crit| > 0.3 -> reasonably robust
##   |rho_crit| < 0.1 -> highly sensitive to unmeasured M-Y confounding
