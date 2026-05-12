# Reference: TwoSampleMR 0.5.11+, MendelianRandomization 0.10+, coloc 5.2+, ieugwasr 1.0+ | Verify API if version differs
# Cis-pQTL MR pipeline for a single drug target with coloc triangulation and PAV-aware sensitivity.
# Worked example: PCSK9 (UKB-PPP Olink) -> coronary artery disease (CARDIoGRAMplusC4D).

library(TwoSampleMR)
library(MendelianRandomization)
library(coloc)
library(ieugwasr)
library(dplyr)

target_gene <- 'PCSK9'
gene_chr <- 1
gene_start <- 55039548   # hg38 PCSK9 start
gene_end <- 55064852     # hg38 PCSK9 end
cis_window_bp <- 500000  # +/- 500 kb per Schmidt 2020 cis-MR convention

pqtl_file <- 'data/ukbppp_pcsk9_full_window.tsv'
outcome_file <- 'data/cad_gwas.tsv'
pav_annotation_file <- 'data/pcsk9_vep_pav.tsv'   # SNP -> consequence from Ensembl VEP
ld_bfile <- '1kg_EUR/EUR'
plink_bin <- genetics.binaRies::get_plink_binary()

pqtl_full <- read.table(pqtl_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
pqtl_cis <- subset(pqtl_full,
    CHR == gene_chr &
    POS >= (gene_start - cis_window_bp) &
    POS <= (gene_end + cis_window_bp))

pqtl_sig <- subset(pqtl_cis, P < 5e-8)
pqtl_sig$f_stat <- (pqtl_sig$BETA / pqtl_sig$SE)^2
pqtl_sig <- subset(pqtl_sig, f_stat >= 10)   # Staiger-Stock 1997 weak-IV floor

pav_consequences <- c('missense_variant', 'stop_gained', 'stop_lost', 'frameshift_variant',
                       'splice_acceptor_variant', 'splice_donor_variant', 'start_lost',
                       'protein_altering_variant')

vep <- read.table(pav_annotation_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
pqtl_sig$is_pav <- pqtl_sig$SNP %in% vep$SNP[vep$Consequence %in% pav_consequences]
cat('Total cis-pQTLs after F-filter:', nrow(pqtl_sig), '\n')
cat('PAV cis-pQTLs:', sum(pqtl_sig$is_pav), '\n')

exposure_dat <- format_data(pqtl_sig, type = 'exposure',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P')

clumped <- ld_clump(
    tibble(rsid = exposure_dat$SNP, pval = exposure_dat$pval.exposure),
    clump_r2 = 0.1, clump_kb = 500,   # cis-window clumping per Schmidt 2020
    plink_bin = plink_bin, bfile = ld_bfile)

exposure_clumped <- subset(exposure_dat, SNP %in% clumped$rsid)
cat('Clumped cis-pQTLs:', nrow(exposure_clumped), '\n')

outcome_dat <- read_outcome_data(outcome_file, snps = exposure_clumped$SNP, sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P')

dat <- harmonise_data(exposure_clumped, outcome_dat, action = 2)

run_cis_panel <- function(d) {
    methods <- if (nrow(d) == 1) 'mr_wald_ratio' else
               if (nrow(d) >= 10) c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median') else
               c('mr_ivw', 'mr_weighted_median')
    mr(d, method_list = methods)
}

primary <- run_cis_panel(dat)
print(primary)

dat$is_pav <- dat$SNP %in% vep$SNP[vep$Consequence %in% pav_consequences]
dat_no_pav <- subset(dat, !is_pav)
cat('Instruments after PAV exclusion:', nrow(dat_no_pav), '\n')
pav_excluded <- if (nrow(dat_no_pav) >= 1) run_cis_panel(dat_no_pav) else NULL
print(pav_excluded)

if (nrow(dat) >= 2) {
    ld <- ld_matrix(dat$SNP, bfile = ld_bfile, plink_bin = plink_bin, with_alleles = FALSE)
    common <- intersect(dat$SNP, rownames(ld))
    d_corr <- subset(dat, SNP %in% common)
    ld_sub <- ld[common, common]
    mr_obj <- mr_input(bx = d_corr$beta.exposure, bxse = d_corr$se.exposure,
                       by = d_corr$beta.outcome, byse = d_corr$se.outcome,
                       correlation = ld_sub)
    result_correl <- mr_ivw(mr_obj, model = 'default', correl = TRUE)
    print(result_correl)
}

gwas_window <- read.table(outcome_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
gwas_window <- subset(gwas_window,
    CHR == gene_chr &
    POS >= (gene_start - cis_window_bp) &
    POS <= (gene_end + cis_window_bp))

shared_snps <- intersect(pqtl_cis$SNP, gwas_window$SNP)
pqtl_coloc <- subset(pqtl_cis, SNP %in% shared_snps)
gwas_coloc <- subset(gwas_window, SNP %in% shared_snps)
gwas_coloc <- gwas_coloc[match(pqtl_coloc$SNP, gwas_coloc$SNP), ]

coloc_res <- coloc.abf(
    dataset1 = list(beta = pqtl_coloc$BETA, varbeta = pqtl_coloc$SE^2,
                    snp = pqtl_coloc$SNP, type = 'quant',
                    N = 54219, sdY = 1),                          # UKB-PPP N; inverse-rank-normal protein expression
    dataset2 = list(beta = gwas_coloc$BETA, varbeta = gwas_coloc$SE^2,
                    snp = gwas_coloc$SNP, type = 'cc',
                    N = 122733, s = 0.34),                        # CARDIoGRAMplusC4D 1000G CAD; ~34% case fraction
    p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)                             # conservative p12 for drug-target

cat('PP.H0:', coloc_res$summary['PP.H0.abf'], '\n')
cat('PP.H3:', coloc_res$summary['PP.H3.abf'], '\n')
cat('PP.H4:', coloc_res$summary['PP.H4.abf'], '\n')

sensitivity_res <- coloc::sensitivity(coloc_res, rule = 'H4 > 0.7')

decision <- list(
    cis_mr_p = primary$pval[primary$method == 'Inverse variance weighted'][1],
    cis_mr_beta = primary$b[primary$method == 'Inverse variance weighted'][1],
    coloc_pp_h4 = unname(coloc_res$summary['PP.H4.abf']),
    n_instruments = nrow(dat),
    n_after_pav_exclusion = nrow(dat_no_pav),
    pav_excluded_p = if (!is.null(pav_excluded)) pav_excluded$pval[1] else NA,
    triangulation_passed = !is.null(pav_excluded) &&
                            !is.na(primary$pval[1]) &&
                            primary$pval[1] < 0.05 / 2923 &&
                            unname(coloc_res$summary['PP.H4.abf']) >= 0.7
)
print(decision)
