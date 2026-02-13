# Reference: R stats (base), ggplot2 3.5+ | Verify API if version differs
## eQTL mediation analysis
##
## Tests whether gene expression mediates the effect of a SNP on disease.
## Demonstrates the mediation package workflow with bootstrap CIs.

library(mediation)

# --- Simulate individual-level data ---
set.seed(42)
n <- 500

genotype <- rbinom(n, 2, 0.3)
age <- rnorm(n, 55, 10)
sex <- rbinom(n, 1, 0.5)
pc1 <- rnorm(n)
pc2 <- rnorm(n)

# SNP affects expression (a path)
expression <- 0.5 * genotype + 0.01 * age - 0.1 * sex + 0.2 * pc1 + rnorm(n, 0, 0.8)

# Expression affects disease risk (b path), with residual direct effect (c' path)
logit_p <- -2 + 0.3 * expression + 0.15 * genotype + 0.02 * age + 0.1 * sex
disease <- rbinom(n, 1, plogis(logit_p))

dat <- data.frame(genotype, expression, disease, age, sex, pc1, pc2)

# --- Step 1: Mediator model (genotype -> expression) ---
med_model <- lm(expression ~ genotype + age + sex + pc1 + pc2, data = dat)
cat('Mediator model (genotype -> expression):\n')
cat('  Genotype coefficient:', round(coef(med_model)['genotype'], 4), '\n')
cat('  P-value:', format.pval(summary(med_model)$coefficients['genotype', 4]), '\n\n')

# --- Step 2: Outcome model (genotype + expression -> disease) ---
out_model <- glm(disease ~ genotype + expression + age + sex + pc1 + pc2,
                 data = dat, family = binomial)
cat('Outcome model (expression -> disease):\n')
cat('  Expression coefficient:', round(coef(out_model)['expression'], 4), '\n')
cat('  Direct genotype coefficient:', round(coef(out_model)['genotype'], 4), '\n\n')

# --- Step 3: Mediation analysis ---
# sims = 1000: Bootstrap iterations for CIs (use 5000 for publication)
med_result <- mediate(med_model, out_model,
                      treat = 'genotype', mediator = 'expression',
                      boot = TRUE, sims = 1000)

cat('Mediation Results:\n')
cat('  ACME (indirect):', round(med_result$d0, 4),
    ' [', round(med_result$d0.ci[1], 4), ',', round(med_result$d0.ci[2], 4), ']\n')
cat('  ADE (direct):', round(med_result$z0, 4),
    ' [', round(med_result$z0.ci[1], 4), ',', round(med_result$z0.ci[2], 4), ']\n')
cat('  Total:', round(med_result$tau.coef, 4), '\n')
cat('  Proportion mediated:', round(med_result$n0, 3), '\n')
cat('  ACME p-value:', format.pval(med_result$d0.p), '\n')

# --- Multiple genes ---
cat('\n--- Testing mediation for multiple genes ---\n')

run_mediation <- function(dat, gene_col, covars = c('age', 'sex', 'pc1', 'pc2')) {
  covar_str <- paste(covars, collapse = ' + ')
  med_form <- as.formula(paste(gene_col, '~ genotype +', covar_str))
  out_form <- as.formula(paste('disease ~ genotype +', gene_col, '+', covar_str))

  med_mod <- lm(med_form, data = dat)
  out_mod <- glm(out_form, data = dat, family = binomial)

  result <- mediate(med_mod, out_mod, treat = 'genotype', mediator = gene_col,
                    boot = TRUE, sims = 500)

  data.frame(gene = gene_col, acme = result$d0, acme_p = result$d0.p,
             ade = result$z0, prop_mediated = result$n0)
}

# Simulate a second gene (no mediation effect)
dat$expression2 <- 0.05 * genotype + rnorm(n, 0, 1.2)

results <- rbind(
  run_mediation(dat, 'expression'),
  run_mediation(dat, 'expression2')
)
results$acme_fdr <- p.adjust(results$acme_p, method = 'BH')
print(results)
