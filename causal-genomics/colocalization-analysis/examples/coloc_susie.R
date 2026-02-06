## SuSiE-coloc for multiple causal variants
##
## When a locus has multiple independent signals, standard coloc.abf
## can give misleading results. SuSiE-coloc tests each pair of
## credible sets between datasets.

library(coloc)
library(susieR)

# --- Simulate data with two causal variants ---
set.seed(42)
n_snps <- 500
positions <- sort(sample(30000000:31000000, n_snps))

# Two causal SNPs
causal1 <- which.min(abs(positions - 30300000))
causal2 <- which.min(abs(positions - 30700000))

# Simulated LD matrix (block diagonal for simplicity)
# In practice, compute from a reference panel with plink --r square
ld_matrix <- diag(n_snps)
for (i in 1:n_snps) {
  for (j in max(1, i - 5):min(n_snps, i + 5)) {
    if (i != j) {
      ld_matrix[i, j] <- 0.8^abs(i - j)
    }
  }
}

gwas_beta <- rnorm(n_snps, 0, 0.01)
gwas_beta[causal1] <- 0.12
gwas_beta[causal2] <- 0.10
gwas_se <- rep(0.025, n_snps)

eqtl_beta <- rnorm(n_snps, 0, 0.02)
# Shared signal at causal1 only
eqtl_beta[causal1] <- 0.35
eqtl_se <- rep(0.04, n_snps)

snp_ids <- paste0('rs', 1:n_snps)

# --- Format for SuSiE ---
gwas_data <- list(
  beta = gwas_beta, varbeta = gwas_se^2,
  snp = snp_ids, position = positions,
  type = 'cc', s = 0.3, N = 50000, LD = ld_matrix
)

eqtl_data <- list(
  beta = eqtl_beta, varbeta = eqtl_se^2,
  snp = snp_ids, position = positions,
  type = 'quant', N = 500, sdY = 1, LD = ld_matrix
)

# --- Run SuSiE on each dataset ---
# L: Maximum number of causal variants to consider
# L = 10 is a reasonable default; increase if many signals expected
susie_gwas <- runsusie(gwas_data, L = 10)
susie_eqtl <- runsusie(eqtl_data, L = 10)

cat('GWAS credible sets found:', length(summary(susie_gwas)$cs), '\n')
cat('eQTL credible sets found:', length(summary(susie_eqtl)$cs), '\n')

# --- Run SuSiE-coloc ---
result <- coloc.susie(susie_gwas, susie_eqtl)

if (!is.null(result$summary)) {
  cat('\nSuSiE-coloc results:\n')
  print(result$summary)
  # Each row: one pair of credible sets tested
  # hit1, hit2: Credible set indices from GWAS and eQTL
  # PP.H4.abf: Posterior probability of shared causal variant for this pair
} else {
  cat('\nNo overlapping credible sets found\n')
}
