# Reference: TwoSampleMR 0.5.11+, MR-PRESSO 1.0+, MendelianRandomization 0.9.0+, mr.raps 0.4+ (GitHub qingyuanzhao/mr.raps), simex 1.8+ | Verify API if version differs
## Postdoc-grade UHP-focused sensitivity battery for two-sample MR
##
## Runs IVW + Egger (with NOME / I^2_GX check and optional SIMEX) + weighted median +
## weighted mode + Cochran Q + Egger intercept + MR-PRESSO (5000 distributions) +
## Steiger + leave-one-out. Produces a STROBE-MR compatible reporting table.
##
## Scope: addresses UHP. For CHP-aware estimation (CAUSE / LHC-MR / LCV),
## see cause_analysis.R or run those tools separately.

library(TwoSampleMR)
library(MRPRESSO)
library(MendelianRandomization)

set.seed(42)
n <- 30

dat <- data.frame(
    SNP = paste0('rs', 1:n),
    beta.exposure = abs(rnorm(n, 0.06, 0.02)),
    se.exposure = rep(0.012, n),
    se.outcome = rep(0.015, n),
    effect_allele.exposure = rep('A', n),
    other_allele.exposure = rep('G', n),
    effect_allele.outcome = rep('A', n),
    other_allele.outcome = rep('G', n),
    eaf.exposure = runif(n, 0.15, 0.85),
    eaf.outcome = runif(n, 0.15, 0.85),
    id.exposure = rep('exposure', n),
    id.outcome = rep('outcome', n),
    exposure = rep('BMI', n),
    outcome = rep('T2D', n),
    samplesize.exposure = rep(500000, n),
    samplesize.outcome = rep(500000, n),
    mr_keep = rep(TRUE, n),
    stringsAsFactors = FALSE)

true_beta <- 0.35
dat$beta.outcome <- dat$beta.exposure * true_beta + rnorm(n, 0, 0.006)
dat$beta.outcome[c(3, 12, 21)] <- dat$beta.outcome[c(3, 12, 21)] + c(0.04, -0.03, 0.05)

dat$f_stat <- (dat$beta.exposure / dat$se.exposure)^2
mean_f <- mean(dat$f_stat)
weak_n <- sum(dat$f_stat < 10)
cat('=== Instrument strength ===\n')
cat('Mean F-statistic:', round(mean_f, 1), '\n')
cat('Weak (F < 10):', weak_n, '/', n, '\n\n')

cat('=== Core MR methods ===\n')
methods <- c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_weighted_mode')
res_mr <- mr(dat, method_list = methods)
print(res_mr[, c('method', 'nsnp', 'b', 'se', 'pval')])

cat('\n=== Heterogeneity (Cochran Q) ===\n')
het <- mr_heterogeneity(dat)
print(het[, c('method', 'Q', 'Q_df', 'Q_pval')])

cat('\n=== Egger intercept (directional UHP) ===\n')
pleio <- mr_pleiotropy_test(dat)
cat('Intercept:', round(pleio$egger_intercept, 5), '\n')
cat('SE:', round(pleio$se, 5), '\n')
cat('P-value:', format.pval(pleio$pval), '\n')

isq <- Isq(dat$beta.exposure, dat$se.exposure)
cat('\nI^2_GX (NOME):', round(isq, 3), '\n')
if (isq >= 0.9) {
    cat('NOME assumption holds; Egger slope reliable\n')
} else if (isq >= 0.6) {
    cat('NOME partial; consider SIMEX correction\n')
} else {
    cat('NOME violated; Egger slope unreliable; prefer MR-RAPS\n')
}

cat('\n=== Steiger directionality ===\n')
steiger <- directionality_test(dat)
cat('Correct direction:', steiger$correct_causal_direction, '\n')
cat('Steiger p-value:', format.pval(steiger$steiger_pval), '\n')

dat_steiger_filt <- steiger_filtering(dat)
n_pass <- sum(dat_steiger_filt$steiger_dir == TRUE)
cat('Per-SNP pass:', n_pass, '/', n, '\n')

cat('\n=== MR-PRESSO ===\n')
presso_input <- data.frame(
    bx = dat$beta.exposure, by = dat$beta.outcome,
    bxse = dat$se.exposure, byse = dat$se.outcome)

presso <- mr_presso(
    BetaOutcome = 'by', BetaExposure = 'bx',
    SdOutcome = 'byse', SdExposure = 'bxse',
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = presso_input, NbDistribution = 5000, SignifThreshold = 0.05)

global_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
outlier_p <- presso$`MR-PRESSO results`$`Outlier Test`$Pvalue
n_outliers <- sum(outlier_p < 0.05, na.rm = TRUE)
distortion_p <- presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
main <- presso$`Main MR results`

cat('Global p:', global_p, '\n')
cat('Outliers detected:', n_outliers, '/', n, '\n')
cat('Distortion p:', distortion_p, '\n')
cat('Raw IVW:', round(main$`Causal Estimate`[1], 4), '\n')
cat('Corrected IVW:', round(main$`Causal Estimate`[2], 4), '\n')

if (n_outliers / n > 0.5) {
    cat('WARNING: >50% pleiotropic; PRESSO majority-valid assumption fails\n')
    cat('         Switch to weighted-mode or CAUSE / LHC-MR\n')
}

cat('\n=== Contamination mixture ===\n')
mr_obj <- mr_input(
    bx = dat$beta.exposure, bxse = dat$se.exposure,
    by = dat$beta.outcome, byse = dat$se.outcome,
    snps = dat$SNP)
conmix <- mr_conmix(mr_obj)
# MendelianRandomization returns S4 objects; use @ to access slots
cat('Estimate:', round(conmix@Estimate, 4), '\n')
cat('95% CI:', round(conmix@CILower, 4), 'to', round(conmix@CIUpper, 4), '\n')

cat('\n=== MR-RAPS ===\n')
# TwoSampleMR::mr_raps wraps mr.raps::mr.raps (GitHub qingyuanzhao/mr.raps; CRAN-archived 2025-03-01)
# Returns a list with $b, $se, $pval, $nsnp; MendelianRandomization does NOT export mr_raps
raps <- TwoSampleMR::mr_raps(b_exp = dat$beta.exposure, b_out = dat$beta.outcome,
                              se_exp = dat$se.exposure, se_out = dat$se.outcome)
cat('Estimate:', round(raps$b, 4), '\n')
cat('SE:', round(raps$se, 4), '\n')
cat('P-value:', format.pval(raps$pval), '\n')

cat('\n=== Leave-one-out (IVW) ===\n')
loo <- mr_leaveoneout(dat)
loo_b <- loo$b[!is.na(loo$b)]
cat('LOO range:', round(min(loo_b), 4), 'to', round(max(loo_b), 4), '\n')

cat('\n=== STROBE-MR reporting summary ===\n')
report_tab <- data.frame(
    Method = c('IVW', 'MR-Egger', 'Weighted median', 'Weighted mode',
               'MR-PRESSO (corrected)', 'Contamination mixture', 'MR-RAPS'),
    Estimate = c(res_mr$b[1], res_mr$b[2], res_mr$b[3], res_mr$b[4],
                 main$`Causal Estimate`[2], conmix@Estimate, raps$b),
    SE = c(res_mr$se[1], res_mr$se[2], res_mr$se[3], res_mr$se[4],
           main$Sd[2], (conmix@CIUpper - conmix@CILower) / (2 * 1.96), raps$se),
    P = c(res_mr$pval[1], res_mr$pval[2], res_mr$pval[3], res_mr$pval[4],
          main$`P-value`[2], NA, raps$pval))
report_tab$Estimate <- round(report_tab$Estimate, 4)
report_tab$SE <- round(report_tab$SE, 4)
print(report_tab)

cat('\nDiagnostics:\n')
cat('  Mean F-statistic:', round(mean_f, 1), '(>=10 strong)\n')
cat('  I^2_GX (NOME):', round(isq, 3), '(>=0.9 Egger reliable)\n')
cat('  Egger intercept p:', format.pval(pleio$pval), '\n')
cat('  Cochran Q p (IVW):', format.pval(het$Q_pval[het$method == 'Inverse variance weighted']), '\n')
cat('  PRESSO global p:', global_p, '\n')
cat('  PRESSO outliers:', n_outliers, '\n')
cat('  Steiger correct direction:', steiger$correct_causal_direction, '\n')

cat('\n=== Interpretation rule ===\n')
all_same_sign <- all(report_tab$Estimate > 0) || all(report_tab$Estimate < 0)
cat('All methods agree on direction:', all_same_sign, '\n')
cat('If LDSC rg(exposure, outcome) >= 0.3 OR biology suggests shared factor:\n')
cat('  Additionally run CAUSE or LHC-MR (see cause_analysis.R)\n')
cat('  UHP-method agreement alone is INSUFFICIENT under CHP\n')
