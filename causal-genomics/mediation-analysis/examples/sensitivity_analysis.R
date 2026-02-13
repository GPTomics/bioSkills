# Reference: R stats (base), ggplot2 3.5+ | Verify API if version differs
## Mediation sensitivity analysis
##
## Tests robustness of mediation results to unmeasured confounding.
## The sequential ignorability assumption is untestable, but sensitivity
## analysis quantifies how much confounding would be needed to nullify
## the mediation effect.

library(mediation)

# --- Simulate data ---
set.seed(42)
n <- 500

genotype <- rbinom(n, 2, 0.3)
age <- rnorm(n, 55, 10)
sex <- rbinom(n, 1, 0.5)

expression <- 0.5 * genotype + 0.01 * age + rnorm(n, 0, 0.8)
logit_p <- -2 + 0.3 * expression + 0.15 * genotype + 0.02 * age
disease <- rbinom(n, 1, plogis(logit_p))

dat <- data.frame(genotype, expression, disease, age, sex)

# --- Fit models and run mediation ---
med_model <- lm(expression ~ genotype + age + sex, data = dat)
out_model <- glm(disease ~ genotype + expression + age + sex, data = dat, family = binomial)

med_result <- mediate(med_model, out_model,
                      treat = 'genotype', mediator = 'expression',
                      boot = TRUE, sims = 1000)

cat('Mediation result:\n')
cat('  ACME:', round(med_result$d0, 4), '\n')
cat('  ACME p-value:', format.pval(med_result$d0.p), '\n')
cat('  Proportion mediated:', round(med_result$n0, 3), '\n\n')

# --- Sensitivity analysis ---
# rho.by: Step size for correlation between residuals
# rho ranges from -1 to 1; at rho = 0, assumes no unmeasured confounding
# The critical rho where ACME crosses 0 indicates robustness
# |critical rho| > 0.3: Reasonably robust
# |critical rho| < 0.1: Sensitive to even small confounding
sens <- medsens(med_result, rho.by = 0.05, effect.type = 'indirect', sims = 1000)

cat('Sensitivity analysis:\n')
summary(sens)

# --- Plot sensitivity ---
pdf('sensitivity_plot.pdf', width = 8, height = 6)
plot(sens, main = 'Sensitivity of ACME to Unmeasured Confounding',
     xlab = 'Sensitivity Parameter (rho)',
     ylab = 'Average Causal Mediation Effect (ACME)')
dev.off()

# --- Alternative: R-squared sensitivity ---
# Interpret in terms of R^2 rather than rho
pdf('sensitivity_rsq.pdf', width = 8, height = 6)
plot(sens, sens.par = 'R2', r.type = 'total',
     main = 'Sensitivity: R-squared Framework')
dev.off()

cat('\nPlots saved: sensitivity_plot.pdf, sensitivity_rsq.pdf\n')
cat('\nInterpretation guide:\n')
cat('  - The x-axis shows hypothetical correlation (rho) between residuals\n')
cat('  - At rho = 0, no unmeasured confounding (standard mediation estimate)\n')
cat('  - The rho where ACME crosses 0 is the critical sensitivity parameter\n')
cat('  - Larger |critical rho| = more robust to confounding\n')
