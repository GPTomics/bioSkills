library(bigsnpr)
library(data.table)

# --- LDpred2 PRS Calculation ---
# LDpred2 uses Bayesian shrinkage for more accurate PRS
# Recommended: LDpred2-auto for automatic hyperparameter tuning

# --- ALTERNATIVE: Use real data ---
# Download example genotypes from UK Biobank or 1000 Genomes
# Download GWAS from GWAS Catalog

# Step 1: Convert plink to bigSNP format (one-time)
# snp_readBed('genotypes.bed')
# obj.bigsnp <- snp_attach('genotypes.rds')

# For demonstration, simulate data
set.seed(42)
n <- 1000   # samples
m <- 10000  # SNPs

message('Simulating genotype data for demonstration...')
G_sim <- matrix(sample(0:2, n * m, replace = TRUE, prob = c(0.25, 0.5, 0.25)), n, m)
map_sim <- data.frame(
    chr = rep(1:22, length.out = m),
    pos = cumsum(sample(1000:5000, m, replace = TRUE)),
    rsid = paste0('rs', 1:m),
    a0 = sample(c('A', 'C', 'G', 'T'), m, replace = TRUE),
    a1 = sample(c('A', 'C', 'G', 'T'), m, replace = TRUE)
)

# Simulate GWAS summary stats
# True betas for 1% of SNPs (causal)
true_beta <- rep(0, m)
causal <- sample(m, m * 0.01)
true_beta[causal] <- rnorm(length(causal), 0, 0.1)

# Simulated phenotype (h2 ~ 0.5)
pheno <- G_sim %*% true_beta + rnorm(n, 0, sqrt(0.5))
pheno_binary <- ifelse(pheno > median(pheno), 1, 0)

# GWAS (simplified)
gwas_beta <- sapply(1:m, function(j) {
    coef(lm(pheno ~ G_sim[, j]))[2]
})
gwas_se <- sapply(1:m, function(j) {
    summary(lm(pheno ~ G_sim[, j]))$coef[2, 2]
})
gwas_p <- 2 * pnorm(-abs(gwas_beta / gwas_se))

sumstats <- data.frame(
    rsid = map_sim$rsid,
    chr = map_sim$chr,
    pos = map_sim$pos,
    a0 = map_sim$a0,
    a1 = map_sim$a1,
    beta = gwas_beta,
    beta_se = gwas_se,
    p = gwas_p,
    n_eff = n
)

message('Simulated GWAS with ', sum(gwas_p < 5e-8), ' genome-wide significant SNPs')

# --- LDpred2-auto (Recommended) ---
# In practice, use real genotype data and compute LD from reference panel
message('\nFor real data, LDpred2 workflow:')
cat('
# 1. Load genotypes
obj.bigsnp <- snp_attach("genotypes.rds")
G <- obj.bigsnp$genotypes
map <- obj.bigsnp$map

# 2. Match GWAS to genotypes
df_beta <- snp_match(sumstats, map, strand_flip = TRUE)

# 3. Compute LD correlation matrix
corr <- snp_cor(G, ind.col = df_beta$`_NUM_ID_`)

# 4. Estimate heritability with LDSC
ldsc <- snp_ldsc2(corr, df_beta)
h2_est <- ldsc[["h2"]]

# 5. Run LDpred2-auto
multi_auto <- snp_ldpred2_auto(
    corr,
    df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.2, 30),
    ncores = 4
)

# 6. Get posterior betas
beta_auto <- sapply(multi_auto, function(x) x$beta_est)

# 7. Compute PRS
prs <- big_prodMat(G, beta_auto)

# 8. Evaluate (if phenotype available)
library(pROC)
auc(pheno_binary ~ prs[, 1])
')

# Simple PRS for demonstration
message('\nSimple PRS (weighted sum):')

# Clump by selecting top SNPs with low p-value
# In practice use proper LD clumping
sig_idx <- which(sumstats$p < 0.01)
message('Using ', length(sig_idx), ' SNPs with p < 0.01')

# Calculate PRS as weighted sum
prs_simple <- G_sim[, sig_idx] %*% sumstats$beta[sig_idx]

# Normalize to Z-scores
prs_z <- scale(prs_simple)[, 1]

# Correlation with true phenotype
cor_with_pheno <- cor(prs_z, pheno)
message('Correlation with phenotype: ', round(cor_with_pheno, 3))

# Risk stratification
# Standard thresholds: 16th percentile (low), 84th percentile (high)
risk_category <- cut(prs_z,
                     breaks = c(-Inf, -1, 1, 2, Inf),
                     labels = c('Low', 'Average', 'High', 'Very High'))
table(risk_category)

message('\nPRS analysis complete')
