# Reference: ggplot2 3.5+ | Verify API if version differs
## Bayesian colocalization with coloc.abf
##
## Tests whether a GWAS signal and eQTL share a causal variant.
## Demonstrates input preparation, running coloc, and interpreting results.

library(coloc)

# --- Simulate locus data (1 Mb region, ~1000 SNPs) ---
set.seed(42)
n_snps <- 1000
positions <- sort(sample(30000000:31000000, n_snps))

# Shared causal variant at position ~30500000
causal_idx <- which.min(abs(positions - 30500000))

# GWAS summary stats (case-control)
gwas_beta <- rnorm(n_snps, 0, 0.02)
gwas_beta[causal_idx] <- 0.15
gwas_se <- rep(0.03, n_snps)
gwas_p <- 2 * pnorm(-abs(gwas_beta / gwas_se))

# eQTL summary stats (quantitative)
eqtl_beta <- rnorm(n_snps, 0, 0.03)
eqtl_beta[causal_idx] <- 0.4
eqtl_se <- rep(0.05, n_snps)
eqtl_p <- 2 * pnorm(-abs(eqtl_beta / eqtl_se))

snp_ids <- paste0('rs', 1:n_snps)

# --- Format coloc input ---
# type = 'cc' for case-control, 'quant' for quantitative
# s = proportion of cases (required for cc)
gwas_input <- list(
  beta = gwas_beta, varbeta = gwas_se^2,
  snp = snp_ids, position = positions,
  type = 'cc', s = 0.3, N = 50000
)

eqtl_input <- list(
  beta = eqtl_beta, varbeta = eqtl_se^2,
  snp = snp_ids, position = positions,
  type = 'quant', N = 500, sdY = 1
)

# --- Run coloc.abf ---
# Default priors: p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
result <- coloc.abf(dataset1 = gwas_input, dataset2 = eqtl_input)

cat('Posterior probabilities:\n')
print(round(result$summary, 4))

# --- Interpret ---
# PP.H4 > 0.8: Strong evidence for shared causal variant
# PP.H3 > 0.8: Distinct causal variants at locus
pp4 <- result$summary['PP.H4.abf']
cat('\nPP.H4 =', round(pp4, 3), '\n')
if (pp4 > 0.8) {
  cat('Strong colocalization: traits likely share a causal variant\n')
} else if (pp4 > 0.5) {
  cat('Suggestive colocalization: consider SuSiE-coloc or larger sample\n')
} else {
  cat('No strong evidence for colocalization\n')
}

# --- Per-SNP posterior probabilities ---
snp_pp <- result$results
top_snps <- snp_pp[order(-snp_pp$SNP.PP.H4), ][1:5, ]
cat('\nTop 5 SNPs by PP.H4:\n')
print(top_snps[, c('snp', 'SNP.PP.H4')])

# --- Prior sensitivity ---
# Check if conclusion depends on prior choice
sensitivity(result, 'H4 > 0.8')
