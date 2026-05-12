# Reference: TwoSampleMR 0.5.11+, MendelianRandomization 0.10+, coloc 5.2+, ieugwasr 1.0+ | Verify API if version differs
# Phenome-wide cis-MR for a single drug target across OpenGWAS outcomes.
# Worked example: PCSK9 cis-pQTL instrument set scanned for on-target adverse effects.
# Recreates the discovery logic behind Schmidt 2017 Lancet Diabetes Endocrinol 5:97 (PCSK9 -> T2D risk).

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)

target_gene <- 'PCSK9'
gene_chr <- 1
gene_start <- 55039548
gene_end <- 55064852
cis_window_bp <- 500000
min_outcome_n <- 50000
min_outcome_population <- 'European'
ld_bfile <- '1kg_EUR/EUR'
plink_bin <- genetics.binaRies::get_plink_binary()

pqtl_full <- read.table('data/ukbppp_pcsk9_full_window.tsv',
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)
pqtl_cis <- subset(pqtl_full,
    CHR == gene_chr &
    POS >= (gene_start - cis_window_bp) &
    POS <= (gene_end + cis_window_bp) &
    P < 5e-8)
pqtl_cis$f_stat <- (pqtl_cis$BETA / pqtl_cis$SE)^2
pqtl_cis <- subset(pqtl_cis, f_stat >= 10)

exposure_dat <- format_data(pqtl_cis, type = 'exposure',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P')

clumped <- ld_clump(
    tibble(rsid = exposure_dat$SNP, pval = exposure_dat$pval.exposure),
    clump_r2 = 0.1, clump_kb = 500,
    plink_bin = plink_bin, bfile = ld_bfile)

exposure_clumped <- subset(exposure_dat, SNP %in% clumped$rsid)
cat('Instrument set:', nrow(exposure_clumped), 'cis-pQTLs\n')

all_outcomes <- available_outcomes()
outcomes_filt <- subset(all_outcomes,
    sample_size >= min_outcome_n &
    population == min_outcome_population &
    !grepl('Sun et al.*pQTL', author, ignore.case = TRUE))   # avoid scanning pQTLs as outcomes

cat('Outcomes to scan:', nrow(outcomes_filt), '\n')

scan_one_outcome <- function(outcome_id) {
    tryCatch({
        outcome_dat <- extract_outcome_data(snps = exposure_clumped$SNP, outcomes = outcome_id)
        if (is.null(outcome_dat) || nrow(outcome_dat) < 2) return(NULL)
        dat <- harmonise_data(exposure_clumped, outcome_dat, action = 2)
        if (nrow(dat) < 2) return(NULL)
        res <- mr(dat, method_list = 'mr_ivw')
        res$outcome_id <- outcome_id
        res$n_snp <- nrow(dat)
        res
    }, error = function(e) NULL)
}

# Sequential loop kept simple; for production parallelise via future.apply::future_lapply
results <- lapply(outcomes_filt$id[1:200], scan_one_outcome)
results_df <- do.call(rbind, Filter(Negate(is.null), results))

n_tests <- nrow(results_df)
results_df$p_bonf <- pmin(results_df$pval * n_tests, 1)
results_df$fdr <- p.adjust(results_df$pval, method = 'BH')

results_df <- merge(results_df,
                    outcomes_filt[, c('id', 'trait', 'author', 'sample_size')],
                    by.x = 'outcome_id', by.y = 'id')

top_hits <- subset(results_df, p_bonf < 0.05)
top_hits <- top_hits[order(top_hits$p_bonf), ]
cat('Bonferroni-significant outcomes:', nrow(top_hits), '\n')
print(top_hits[, c('trait', 'b', 'se', 'pval', 'p_bonf', 'n_snp')])

write.table(results_df, file = 'pcsk9_phewas_mr.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

# Forest-plot input: each row is a phenotype with effect estimate, 95% CI, p_bonf
# Plot externally with ggplot2 or forestplot package (see data-visualization/ggplot2-fundamentals)
forest_input <- top_hits %>%
    transmute(trait,
              beta = b,
              lci = b - 1.96 * se,
              uci = b + 1.96 * se,
              p_bonf,
              n_snp)
write.table(forest_input, file = 'pcsk9_phewas_forest.tsv',
            sep = '\t', row.names = FALSE, quote = FALSE)
