# Reference: GenomicSEM 0.0.5+, lavaan 0.6+, MTAG 1.0+ | Verify API if version differs
#
# Common-factor model + common-factor GWAS with Q_SNP heterogeneity (Grotzinger 2019).
# Produces:
#   - LDSC genetic covariance (S) and sampling covariance (V)
#   - Common-factor CFA with CFI / RMSEA / SRMR
#   - Common-factor GWAS classifying SNPs as factor-mediated vs Q_SNP-heterogeneous
#
# Reference paths assume HapMap3-munged sumstats and the alkesgroup EUR LD-score panel.

library(GenomicSEM)
library(lavaan)
library(Matrix)

trait_files <- c(
    'data/trait1.sumstats.gz',
    'data/trait2.sumstats.gz',
    'data/trait3.sumstats.gz'
)

trait_names <- c('trait1', 'trait2', 'trait3')

sample_prev <- c(0.5, 0.5, NA)
population_prev <- c(0.05, 0.05, NA)

ld_dir <- 'reference/eur_w_ld_chr/'
wld_dir <- 'reference/eur_w_ld_chr/'

# Step 1: LDSC produces S (genetic covariance) and V (sampling covariance with overlap correction)
ldsc_output <- ldsc(
    traits = trait_files,
    sample.prev = sample_prev,
    population.prev = population_prev,
    ld = ld_dir,
    wld = wld_dir,
    trait.names = trait_names,
    stand = FALSE
)

cat('Genetic covariance S:\n')
print(ldsc_output$S)
cat('\nSampling covariance V eigenvalues (must be > 0):\n')
print(eigen(ldsc_output$V)$values)

if (any(eigen(ldsc_output$V)$values < 0)) {
    stop('V matrix not positive definite; check per-trait mean chi-square > 1.02 and re-munge.')
}

# Step 2: Common-factor CFA
cf_fit <- commonfactor(covstruc = ldsc_output, estimation = 'DWLS')

cat('\n--- Common-factor model fit ---\n')
print(cf_fit$modelfit)

cat('\n--- Standardized loadings ---\n')
print(cf_fit$results)

cfi <- cf_fit$modelfit$CFI
rmsea <- cf_fit$modelfit$RMSEA

if (cfi < 0.9 || rmsea > 0.08) {
    warning('Poor model fit: CFI=', round(cfi, 3), ' RMSEA=', round(rmsea, 3),
            '. Consider ESEM or multi-factor model.')
}

# Step 3: Prepare per-SNP sumstats for common-factor GWAS
ref_file <- 'reference/reference.1000G.maf.0.005.txt'

ss_input <- sumstats(
    files = c('data/trait1.raw.txt', 'data/trait2.raw.txt', 'data/trait3.raw.txt'),
    ref = ref_file,
    trait.names = trait_names,
    se.logit = c(TRUE, TRUE, FALSE),   # case-control on logit scale; continuous trait FALSE
    OLS = c(FALSE, FALSE, TRUE),
    linprob = c(FALSE, FALSE, FALSE),
    N = c(150000, 200000, 175000),
    info.filter = 0.9,                 # standard imputation INFO threshold
    maf.filter = 0.01                  # standard MAF threshold
)

# Step 4: Common-factor GWAS with parallel execution
# commonfactorGWAS always reports Q + Q_pval (Q_SNP heterogeneity) in the output;
# no `Q_SNP=TRUE` argument is required (it is on by default in modern releases).
cfgwas <- commonfactorGWAS(
    covstruc = ldsc_output,
    SNPs = ss_input,
    estimation = 'DWLS',
    parallel = TRUE,
    cores = 8
)

# Step 5: Classify SNPs. Output columns include est, se_c, Z_Estimate, Pval_Estimate
# (factor effect) and Q, Q_df, Q_pval (per-SNP heterogeneity = Q_SNP in the literature).
factor_threshold <- 5e-08
n_factor_sig <- sum(cfgwas$Pval_Estimate < factor_threshold, na.rm = TRUE)
qsnp_threshold <- 0.05 / max(n_factor_sig, 1)  # Bonferroni within discovered factor SNPs

cfgwas$factor_sig <- cfgwas$Pval_Estimate < factor_threshold
cfgwas$qsnp_sig <- cfgwas$Q_pval < qsnp_threshold
cfgwas$class <- ifelse(cfgwas$factor_sig & !cfgwas$qsnp_sig, 'common_factor',
                ifelse(cfgwas$factor_sig & cfgwas$qsnp_sig, 'heterogeneous',
                       'non_significant'))

cat('\n--- SNP classification counts ---\n')
print(table(cfgwas$class))

write.table(
    subset(cfgwas, class == 'common_factor'),
    'results/common_factor_snps.tsv',
    sep = '\t', quote = FALSE, row.names = FALSE
)
write.table(
    subset(cfgwas, class == 'heterogeneous'),
    'results/qsnp_heterogeneous_snps.tsv',
    sep = '\t', quote = FALSE, row.names = FALSE
)

# Step 6 (optional): User model -- single-factor confirmed via lavaan syntax.
# With only 3 input traits a two-factor model is under-identified (each factor
# needs >= 3 indicators OR an anchor loading fixed to 1 with factor variance free).
# Below replicates the common factor via explicit syntax for parameterisation control;
# scale F by free loadings (NA*) plus unit factor variance (F~~1*F).
user_syntax <- '
    F =~ NA*trait1 + trait2 + trait3
    F ~~ 1*F
'

user_fit <- usermodel(covstruc = ldsc_output, model = user_syntax, estimation = 'DWLS')

cat('\n--- User model (explicit single-factor) fit ---\n')
print(user_fit$modelfit)
print(user_fit$results)

# AIC / BIC for nested-model comparison
aic_cf <- cf_fit$modelfit$AIC
aic_user <- user_fit$modelfit$AIC
cat('\nAIC: common-factor=', aic_cf, ' user-syntax=', aic_user, '\n')
cat('Lower AIC preferred.\n')
