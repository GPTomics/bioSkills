# Reference: ggplot2 3.5+ | Verify API if version differs
## SuSiE fine-mapping from GWAS summary statistics
##
## Demonstrates susie_rss with a simulated locus: computes PIPs,
## extracts credible sets, and reports results.

library(susieR)

# --- Simulate a locus with 2 causal variants ---
set.seed(42)
n_snps <- 500
n_samples <- 50000

# Simulated LD matrix (block structure)
ld_matrix <- diag(n_snps)
for (i in 1:n_snps) {
  for (j in max(1, i - 10):min(n_snps, i + 10)) {
    if (i != j) {
      ld_matrix[i, j] <- 0.9^abs(i - j)
    }
  }
}

# Two causal variants
causal1 <- 100
causal2 <- 350

# Generate Z-scores: large Z at causal variants, small elsewhere
z_scores <- rnorm(n_snps, 0, 1)
z_scores[causal1] <- 6.5
z_scores[(causal1 - 3):(causal1 + 3)] <- z_scores[(causal1 - 3):(causal1 + 3)] + c(3, 4, 5.5, 6.5, 5.2, 3.8, 2.5)
z_scores[causal2] <- 5.0
z_scores[(causal2 - 2):(causal2 + 2)] <- z_scores[(causal2 - 2):(causal2 + 2)] + c(2.5, 4.0, 5.0, 3.5, 2.0)

snp_ids <- paste0('rs', 1:n_snps)

# --- Run SuSiE ---
# L = 10: Max causal variants to search for (SuSiE prunes unused effects)
fit <- susie_rss(z = z_scores, R = ld_matrix, n = n_samples, L = 10)

# --- Credible sets ---
cs <- fit$sets$cs
cat('Number of credible sets found:', length(cs), '\n\n')

for (i in seq_along(cs)) {
  cat(sprintf('Credible set %d:\n', i))
  cat('  Size:', length(cs[[i]]), 'variants\n')
  cat('  Purity (min |correlation|):', round(fit$sets$purity[i, 1], 3), '\n')
  cat('  Top variant:', snp_ids[cs[[i]][which.max(fit$pip[cs[[i]]])]], '\n')
  cat('  Max PIP in set:', round(max(fit$pip[cs[[i]]]), 4), '\n')

  # Check if true causal is captured
  has_causal1 <- causal1 %in% cs[[i]]
  has_causal2 <- causal2 %in% cs[[i]]
  if (has_causal1) cat('  ** Contains true causal 1 **\n')
  if (has_causal2) cat('  ** Contains true causal 2 **\n')
  cat('\n')
}

# --- Top PIPs ---
pip_df <- data.frame(SNP = snp_ids, PIP = fit$pip)
pip_df <- pip_df[order(-pip_df$PIP), ]

cat('Top 10 variants by PIP:\n')
print(head(pip_df, 10))

# --- Summary ---
cat('\nVariants with PIP > 0.5:', sum(fit$pip > 0.5), '\n')
cat('Variants with PIP > 0.95:', sum(fit$pip > 0.95), '\n')

cat('\nTrue causal 1 (variant', causal1, ') PIP:', round(fit$pip[causal1], 4), '\n')
cat('True causal 2 (variant', causal2, ') PIP:', round(fit$pip[causal2], 4), '\n')
