# Reference: TwoSampleMR 0.6+, MVMR 0.4+ | Verify API if version differs
## MR-mediation via multivariable Mendelian randomization
##
## Estimates total effect (univariable MR for E -> Y) minus direct effect
## (MVMR with E and M as joint exposures) to derive the indirect effect
## of E on Y through mediator M.
##
## Sanderson 2021 IJE 50:1651 (MVMR-mediation tutorial)
## Sanderson 2019 IJE 48:713 (conditional F-statistic; > 10 per exposure required)
##
## Install:
##   install.packages('TwoSampleMR',
##                    repos=c('https://mrcieu.r-universe.dev','https://cloud.r-project.org'))
##   remotes::install_github('WSpiller/MVMR')

library(TwoSampleMR)
library(MVMR)

set.seed(42)

## --- Simulate independent SNP-level summary statistics for E and M ---
## Real workflow: extract genome-wide-significant SNPs for E and M from
## the IEU OpenGWAS catalog (extract_instruments + extract_outcome_data),
## harmonise, and run mr_ivw for total and ivw_mvmr for direct.

n_snps <- 80
beta_E <- rnorm(n_snps, 0, 0.05)      ## SNP -> E effect
beta_M <- 0.3 * beta_E + rnorm(n_snps, 0, 0.03)  ## SNP -> M effect (some via E)
se_E <- abs(rnorm(n_snps, 0.01, 0.002))
se_M <- abs(rnorm(n_snps, 0.01, 0.002))

## True causal model: E -> Y has total effect 0.4, with 60% going through M
## Direct E -> Y: 0.16; M -> Y: 0.8
beta_Y <- 0.16 * beta_E + 0.8 * beta_M + rnorm(n_snps, 0, 0.005)
se_Y <- abs(rnorm(n_snps, 0.005, 0.001))

snps <- paste0('rs', sprintf('%06d', sample(1e6, n_snps)))

## --- Univariable MR for total effect (E -> Y) ---
## Construct a TwoSampleMR-style data frame
harm_total <- data.frame(
  SNP=snps,
  beta.exposure=beta_E, se.exposure=se_E,
  beta.outcome=beta_Y, se.outcome=se_Y,
  effect_allele.exposure='A', other_allele.exposure='G',
  effect_allele.outcome='A', other_allele.outcome='G',
  eaf.exposure=runif(n_snps, 0.1, 0.5),
  eaf.outcome=runif(n_snps, 0.1, 0.5),
  exposure='E', outcome='Y', id.exposure='E', id.outcome='Y',
  mr_keep=TRUE,
  pval.exposure=runif(n_snps, 1e-12, 1e-7),
  pval.outcome=runif(n_snps, 1e-3, 0.5)
)

total <- mr_ivw(b_exp=harm_total$beta.exposure, b_out=harm_total$beta.outcome,
                se_exp=harm_total$se.exposure, se_out=harm_total$se.outcome)
cat('--- Univariable MR (Total Effect E -> Y) ---\n')
cat('  Total beta:', round(total$b, 4), ' SE:', round(total$se, 4),
    ' p:', format.pval(total$pval), '\n')
cat('  95% CI: [', round(total$b - 1.96*total$se, 4), ',',
    round(total$b + 1.96*total$se, 4), ']\n\n')

## --- MVMR for direct effect (E and M joint exposures -> Y) ---
mvmr_dat <- format_mvmr(
  BXGs=cbind(beta_E, beta_M),
  BYG=beta_Y,
  seBXGs=cbind(se_E, se_M),
  seBYG=se_Y,
  RSID=snps
)

## Conditional F-statistic per exposure (Sanderson 2019)
## gencov=0 assumes no genetic covariance between E and M instruments;
## for non-overlapping samples this is reasonable. For overlapping samples,
## estimate gencov from LDSC cross-trait genetic covariance.
fstat <- strength_mvmr(mvmr_dat, gencov=0)
cat('--- Conditional F-statistics (require > 10 per exposure) ---\n')
print(fstat)

## Q-statistic for pleiotropy
qstat <- pleiotropy_mvmr(mvmr_dat, gencov=0)
cat('\n--- MVMR Q-statistic for instrument heterogeneity ---\n')
print(qstat)

## Fit MVMR-IVW
mvmr_fit <- ivw_mvmr(mvmr_dat)
cat('\n--- MVMR (Direct Effects) ---\n')
print(mvmr_fit)

direct_E <- mvmr_fit[1, 'Estimate']
direct_E_se <- mvmr_fit[1, 'Std. Error']

## --- Indirect effect via delta method ---
## indirect = total - direct; assume independence -> SE^2 = total_se^2 + direct_se^2
## (slightly conservative because total and direct are positively correlated;
## bootstrap CIs from individual-level data are preferred when available)
indirect <- total$b - direct_E
indirect_se <- sqrt(total$se^2 + direct_E_se^2)
indirect_ci <- indirect + c(-1.96, 1.96) * indirect_se
prop_med <- indirect / total$b

cat('\n--- Indirect Effect (Total - Direct) ---\n')
cat('  Indirect:', round(indirect, 4), ' SE:', round(indirect_se, 4), '\n')
cat('  95% CI: [', round(indirect_ci[1], 4), ',', round(indirect_ci[2], 4), ']\n')
cat('  Proportion mediated:', round(prop_med, 3), '\n\n')

## --- Decision rules ---
## If fstat conditional F < 10 for either E or M: report weak-instrument bias warning
## If qstat Q-pval < 0.05: substantial pleiotropy -> use ivw_mvmr robust alternatives
##   (qhet_mvmr) or report indirect effect as exploratory
## If indirect CI crosses 0: insufficient evidence for mediation through M
## If indirect / total > 1 or < 0: suspect Steiger reversal (M may cause E)
##   or pleiotropic instruments; investigate before reporting

if (any(fstat < 10)) {
  cat('WARNING: conditional F < 10 detected. Weak-instrument bias possible.\n')
  cat('         Consider Q-statistic-adjusted methods (qhet_mvmr).\n')
}
