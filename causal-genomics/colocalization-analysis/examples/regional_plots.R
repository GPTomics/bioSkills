## Regional association and LocusCompare plots for colocalization
##
## Visualizes GWAS and eQTL signals at a locus to accompany coloc results.

library(ggplot2)
library(patchwork)

# --- Simulate locus data ---
set.seed(42)
n_snps <- 500
positions <- sort(sample(30000000:31000000, n_snps))
causal_idx <- which.min(abs(positions - 30500000))

gwas_p <- runif(n_snps, 0.001, 1)
gwas_p[causal_idx] <- 1e-12
gwas_p[(causal_idx - 3):(causal_idx + 3)] <- runif(7, 1e-8, 1e-5)

eqtl_p <- runif(n_snps, 0.01, 1)
eqtl_p[causal_idx] <- 5e-10
eqtl_p[(causal_idx - 3):(causal_idx + 3)] <- runif(7, 1e-6, 1e-3)

df <- data.frame(
  SNP = paste0('rs', 1:n_snps),
  POS = positions,
  GWAS_P = gwas_p,
  eQTL_P = eqtl_p
)

# --- Regional association plots (stacked) ---
p_gwas <- ggplot(df, aes(x = POS / 1e6, y = -log10(GWAS_P))) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = -log10(5e-8), linetype = 'dashed', color = 'red', alpha = 0.5) +
  labs(x = NULL, y = '-log10(P)', title = 'GWAS') +
  theme_minimal() +
  theme(axis.text.x = element_blank())

p_eqtl <- ggplot(df, aes(x = POS / 1e6, y = -log10(eQTL_P))) +
  geom_point(alpha = 0.5, size = 1, color = 'steelblue') +
  labs(x = 'Position (Mb)', y = '-log10(P)', title = 'eQTL') +
  theme_minimal()

p_regional <- p_gwas / p_eqtl
ggsave('regional_association.pdf', p_regional, width = 8, height = 6)

# --- LocusCompare plot ---
p_compare <- ggplot(df, aes(x = -log10(GWAS_P), y = -log10(eQTL_P))) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = 'lm', se = FALSE, linetype = 'dashed', color = 'grey50', linewidth = 0.5) +
  labs(x = '-log10(P) GWAS', y = '-log10(P) eQTL', title = 'LocusCompare') +
  theme_minimal()

ggsave('locuscompare.pdf', p_compare, width = 6, height = 6)

# --- Color by LD with lead SNP ---
# Simulated r2 values (in practice, compute from reference panel)
df$r2_lead <- exp(-abs(df$POS - positions[causal_idx]) / 50000)

p_ld <- ggplot(df, aes(x = POS / 1e6, y = -log10(GWAS_P), color = r2_lead)) +
  geom_point(size = 1.2) +
  scale_color_gradientn(
    colors = c('darkblue', 'lightblue', 'green', 'orange', 'red'),
    name = expression(r^2)
  ) +
  geom_hline(yintercept = -log10(5e-8), linetype = 'dashed', color = 'red', alpha = 0.3) +
  labs(x = 'Position (Mb)', y = '-log10(P)', title = 'GWAS (colored by LD with lead)') +
  theme_minimal()

ggsave('regional_ld.pdf', p_ld, width = 8, height = 4)

cat('Plots saved: regional_association.pdf, locuscompare.pdf, regional_ld.pdf\n')
