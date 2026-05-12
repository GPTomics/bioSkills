# Reference: TwoSampleMR 0.5.11+, MendelianRandomization 0.10+, MR-PRESSO 1.0+, mr.raps 0.4+ (GitHub qingyuanzhao/mr.raps), ieugwasr 1.0+ | Verify API if version differs
## Two-sample MR with full sensitivity battery.
##
## Demonstrates: instrument selection from local GWAS, F-statistic filtering (from exposure),
## local plink clumping, harmonization with action=2, primary IVW + Egger + median + mode,
## MR-RAPS, MR-PRESSO outlier/distortion, Steiger directionality, and STROBE-MR-style summary.

library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)

set.seed(42)
n_snps <- 60

exposure_gwas <- data.frame(
    SNP = paste0('rs', 1:n_snps),
    BETA = rnorm(n_snps, 0.05, 0.02),
    SE = runif(n_snps, 0.005, 0.012),
    A1 = sample(c('A', 'C', 'G', 'T'), n_snps, replace = TRUE),
    A2 = sample(c('A', 'C', 'G', 'T'), n_snps, replace = TRUE),
    EAF = runif(n_snps, 0.15, 0.85),
    N = rep(50000, n_snps),
    stringsAsFactors = FALSE
)
exposure_gwas$P <- 2 * pnorm(-abs(exposure_gwas$BETA / exposure_gwas$SE))
exposure_gwas <- subset(exposure_gwas, P < 5e-08)

outcome_gwas <- data.frame(
    SNP = exposure_gwas$SNP,
    BETA = exposure_gwas$BETA * 0.4 + rnorm(nrow(exposure_gwas), 0, 0.01),
    SE = runif(nrow(exposure_gwas), 0.008, 0.018),
    A1 = exposure_gwas$A1, A2 = exposure_gwas$A2,
    EAF = exposure_gwas$EAF + rnorm(nrow(exposure_gwas), 0, 0.02),
    N = rep(80000, nrow(exposure_gwas)),
    stringsAsFactors = FALSE
)
outcome_gwas$P <- 2 * pnorm(-abs(outcome_gwas$BETA / outcome_gwas$SE))

write.table(exposure_gwas, 'exposure.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(outcome_gwas, 'outcome.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

exposure_dat <- read_exposure_data(
    filename = 'exposure.tsv', sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P', samplesize_col = 'N'
)

exposure_dat$f_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2  # F = (beta/se)^2, computed from EXPOSURE
exposure_dat <- subset(exposure_dat, f_stat >= 10)  # Staiger-Stock 1997 weak-IV heuristic; debated for one-sample (Zhao 2020)
cat('Mean F:', round(mean(exposure_dat$f_stat), 1), '| n SNPs after F filter:', nrow(exposure_dat), '\n')

# Skipping live clumping in demo. Production code:
# clumped <- ieugwasr::ld_clump(
#     dplyr::tibble(rsid = exposure_dat$SNP, pval = exposure_dat$pval.exposure),
#     clump_r2 = 0.001, clump_kb = 10000,
#     plink_bin = genetics.binaRies::get_plink_binary(),
#     bfile = '1kg_EUR/EUR'
# )
# exposure_dat <- subset(exposure_dat, SNP %in% clumped$rsid)

outcome_dat <- read_outcome_data(
    filename = 'outcome.tsv', snps = exposure_dat$SNP, sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P', samplesize_col = 'N'
)

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)  # infer from EAF; drops MAF~0.5 palindromes
cat('SNPs after harmonization:', nrow(dat), '\n')

primary <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression',
                                    'mr_weighted_median', 'mr_weighted_mode'))
cat('\n--- Primary MR ---\n')
print(primary[, c('method', 'nsnp', 'b', 'se', 'pval')])

het <- mr_heterogeneity(dat)
cat('\n--- Heterogeneity ---\n')
print(het[, c('method', 'Q', 'Q_df', 'Q_pval')])

pleio <- mr_pleiotropy_test(dat)
cat('\n--- Egger intercept ---\n')
cat('intercept:', signif(pleio$egger_intercept, 3),
    '| se:', signif(pleio$se, 3), '| p:', signif(pleio$pval, 3), '\n')

# NOME / I^2_GX check; below 0.9 invalidates Egger without SIMEX correction (Bowden 2016 Stat Med 35:1880)
mr_obj <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                   by = dat$beta.outcome, byse = dat$se.outcome)
isq <- TwoSampleMR::Isq(dat$beta.exposure, dat$se.exposure)
cat('I^2_GX:', round(isq, 3), if (isq < 0.9) ' -- NOME VIOLATED; SIMEX-correct Egger\n' else ' -- NOME OK\n')

# MR-RAPS via TwoSampleMR wrapper (mr.raps CRAN-archived 2025-03-01; install from github qingyuanzhao/mr.raps)
# TwoSampleMR::mr_raps takes four numeric vectors and returns a list with $b, $se, $pval, $nsnp
raps <- TwoSampleMR::mr_raps(b_exp = dat$beta.exposure, b_out = dat$beta.outcome,
                              se_exp = dat$se.exposure, se_out = dat$se.outcome)
cat('\n--- MR-RAPS (Huber) ---\n')
cat('beta:', signif(raps$b, 3), '| se:', signif(raps$se, 3),
    '| p:', signif(raps$pval, 3), '\n')

# MR-PRESSO: NbDistribution >= 10000 for publication-grade p-value precision (default 1000 is exploratory)
presso <- mr_presso(
    BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome', SdExposure = 'se.exposure',
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = dat, NbDistribution = 10000, SignifThreshold = 0.05
)
cat('\n--- MR-PRESSO ---\n')
cat('Global RSSobs:', signif(presso$`MR-PRESSO results`$`Global Test`$RSSobs, 3),
    '| p:', signif(presso$`MR-PRESSO results`$`Global Test`$Pvalue, 3), '\n')
if (!is.null(presso$`MR-PRESSO results`$`Distortion Test`)) {
    cat('Distortion p:', signif(presso$`MR-PRESSO results`$`Distortion Test`$Pvalue, 3), '\n')
}

# Steiger directionality with Hemani Tilling 2022 confounder caveat: heuristic only, not definitive
steiger <- directionality_test(dat)
cat('\n--- Steiger ---\n')
cat('correct direction:', steiger$correct_causal_direction,
    '| p:', signif(steiger$steiger_pval, 3), '\n')

loo <- mr_leaveoneout(dat)
cat('\n--- Leave-one-out range of IVW estimate ---\n')
cat('min:', signif(min(loo$b), 3), '| max:', signif(max(loo$b), 3), '\n')

cat('\n--- STROBE-MR summary ---\n')
cat('Exposure SNPs:', nrow(dat), '| Mean F:', round(mean(exposure_dat$f_stat), 1), '\n')
cat('Primary IVW beta:', signif(primary$b[primary$method == 'Inverse variance weighted'], 3),
    '| p:', signif(primary$pval[primary$method == 'Inverse variance weighted'], 3), '\n')
cat('Egger intercept p:', signif(pleio$pval, 3),
    '| MR-PRESSO global p:', signif(presso$`MR-PRESSO results`$`Global Test`$Pvalue, 3), '\n')

file.remove('exposure.tsv', 'outcome.tsv')
