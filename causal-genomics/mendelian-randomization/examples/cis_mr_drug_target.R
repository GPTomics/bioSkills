# Reference: TwoSampleMR 0.5.11+, MendelianRandomization 0.10+, coloc 5.2+, ieugwasr 1.0+ | Verify API if version differs
## Cis-MR drug-target analysis: cis-pQTL exposure -> outcome with colocalization triangulation.
##
## Demonstrates the Schmidt 2020 Nat Commun 11:3255 framework:
## restrict instruments to +/- 500 kb of the target gene, clump at r2 < 0.1 within window
## (looser than polygenic MR's r2 < 0.001 to retain power), run IVW + Wald ratio,
## then coloc.abf to confirm shared causal variant (PP.H4 >= 0.7 required for drug-target claim).

library(TwoSampleMR)
library(MendelianRandomization)
library(coloc)

set.seed(7)
gene_start <- 50000000
gene_end <-   50050000
cis_window <- 500000  # 500 kb either side per Schmidt 2020 drug-target convention
n_snps <- 40

cis_pqtl <- data.frame(
    SNP = paste0('rs', 1:n_snps),
    CHR = rep(1, n_snps),
    POS = sample(seq(gene_start - cis_window, gene_end + cis_window, by = 100), n_snps),
    BETA = rnorm(n_snps, 0.15, 0.05),
    SE = runif(n_snps, 0.02, 0.04),
    A1 = sample(c('A', 'C', 'G', 'T'), n_snps, replace = TRUE),
    A2 = sample(c('A', 'C', 'G', 'T'), n_snps, replace = TRUE),
    EAF = runif(n_snps, 0.10, 0.90),
    N = rep(54000, n_snps),
    stringsAsFactors = FALSE
)
cis_pqtl$P <- 2 * pnorm(-abs(cis_pqtl$BETA / cis_pqtl$SE))

outcome_assoc <- data.frame(
    SNP = cis_pqtl$SNP, CHR = cis_pqtl$CHR, POS = cis_pqtl$POS,
    BETA = cis_pqtl$BETA * 0.5 + rnorm(n_snps, 0, 0.03),
    SE = runif(n_snps, 0.025, 0.05),
    A1 = cis_pqtl$A1, A2 = cis_pqtl$A2,
    EAF = cis_pqtl$EAF + rnorm(n_snps, 0, 0.02),
    N = rep(250000, n_snps),
    stringsAsFactors = FALSE
)
outcome_assoc$P <- 2 * pnorm(-abs(outcome_assoc$BETA / outcome_assoc$SE))

# Window restriction: instruments must lie within +/- 500 kb of gene boundary
in_window <- with(cis_pqtl, POS >= gene_start - cis_window & POS <= gene_end + cis_window)
cis_pqtl <- cis_pqtl[in_window, ]
cat('cis-pQTLs in window:', nrow(cis_pqtl), '\n')

instruments_raw <- subset(cis_pqtl, P < 5e-08)
instruments_raw$f_stat <- (instruments_raw$BETA / instruments_raw$SE)^2
instruments_raw <- subset(instruments_raw, f_stat >= 10)  # weak-IV filter
cat('Genome-wide-sig cis instruments after F filter:', nrow(instruments_raw), '\n')

# Production clumping (r2 < 0.1 within the 1 Mb window; cis-MR looser than polygenic r2 < 0.001):
# clumped <- ieugwasr::ld_clump(
#     dplyr::tibble(rsid = instruments_raw$SNP, pval = instruments_raw$P),
#     clump_r2 = 0.1, clump_kb = 250,
#     plink_bin = genetics.binaRies::get_plink_binary(),
#     bfile = '1kg_EUR/EUR'
# )
# instruments_raw <- subset(instruments_raw, SNP %in% clumped$rsid)

write.table(instruments_raw, 'cis_pqtl.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(outcome_assoc,   'outcome_locus.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

exposure_dat <- read_exposure_data(
    'cis_pqtl.tsv', sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P', samplesize_col = 'N'
)
outcome_dat <- read_outcome_data(
    'outcome_locus.tsv', snps = exposure_dat$SNP, sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P', samplesize_col = 'N'
)

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
cat('SNPs after harmonization:', nrow(dat), '\n')

# Wald ratio per SNP if only one instrument survives clumping; IVW otherwise
if (nrow(dat) == 1) {
    wr <- mr(dat, method_list = 'mr_wald_ratio')
    cat('\n--- Wald ratio (single instrument) ---\n')
    print(wr[, c('method', 'b', 'se', 'pval')])
} else {
    res <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median'))
    cat('\n--- cis-MR primary ---\n')
    print(res[, c('method', 'nsnp', 'b', 'se', 'pval')])
}

# Colocalization triangulation -- required for drug-target claim
# PP.H4 >= 0.7 indicates shared causal variant; <0.5 indicates likely distinct signals
all_snps <- merge(cis_pqtl, outcome_assoc, by = 'SNP', suffixes = c('.exp', '.out'))
coloc_input_exp <- list(
    beta = all_snps$BETA.exp, varbeta = all_snps$SE.exp^2,
    snp = all_snps$SNP, position = all_snps$POS.exp,
    type = 'quant', N = all_snps$N.exp[1], MAF = pmin(all_snps$EAF.exp, 1 - all_snps$EAF.exp)
)
coloc_input_out <- list(
    beta = all_snps$BETA.out, varbeta = all_snps$SE.out^2,
    snp = all_snps$SNP, position = all_snps$POS.out,
    type = 'quant', N = all_snps$N.out[1], MAF = pmin(all_snps$EAF.out, 1 - all_snps$EAF.out)
)

coloc_res <- coloc.abf(coloc_input_exp, coloc_input_out,
                       p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)  # coloc default priors
cat('\n--- Coloc triangulation ---\n')
cat('PP.H0 (none):', signif(coloc_res$summary['PP.H0.abf'], 3), '\n')
cat('PP.H1 (exp only):', signif(coloc_res$summary['PP.H1.abf'], 3), '\n')
cat('PP.H2 (out only):', signif(coloc_res$summary['PP.H2.abf'], 3), '\n')
cat('PP.H3 (distinct causals):', signif(coloc_res$summary['PP.H3.abf'], 3), '\n')
cat('PP.H4 (shared causal):', signif(coloc_res$summary['PP.H4.abf'], 3),
    if (coloc_res$summary['PP.H4.abf'] >= 0.7) ' -- DRUG-TARGET SUPPORTED\n' else ' -- NOT SUPPORTED\n')

cat('\nInterpretation:\n')
cat('cis-MR estimates the effect of genetically perturbed target protein on outcome.\n')
cat('Coloc PP.H4 >= 0.7 confirms the same causal variant drives both signals,\n')
cat('ruling out the alternative that distinct nearby variants drive each independently.\n')

file.remove('cis_pqtl.tsv', 'outcome_locus.tsv')
