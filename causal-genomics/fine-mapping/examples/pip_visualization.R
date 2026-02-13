# Reference: ggplot2 3.5+ | Verify API if version differs
## PIP visualization for fine-mapping results
##
## Generates PIP plots, combined GWAS + PIP panels, and credible set
## membership plots.

library(ggplot2)
library(patchwork)

# --- Simulate fine-mapping results ---
set.seed(42)
n_snps <- 500
positions <- sort(sample(30000000:31000000, n_snps))

causal_idx <- which.min(abs(positions - 30500000))

gwas_p <- runif(n_snps, 0.01, 1)
gwas_p[causal_idx] <- 1e-15
gwas_p[(causal_idx - 5):(causal_idx + 5)] <- runif(11, 1e-10, 1e-4)

pip <- rbeta(n_snps, 0.5, 10)
pip[causal_idx] <- 0.98
pip[causal_idx - 1] <- 0.65
pip[causal_idx + 1] <- 0.45
pip[causal_idx - 2] <- 0.12

# Credible set: variants with cumulative PIP >= 0.95
cs_variants <- c(causal_idx, causal_idx - 1, causal_idx + 1)

df <- data.frame(
  SNP = paste0('rs', 1:n_snps),
  POS = positions,
  P = gwas_p,
  PIP = pip,
  in_cs = 1:n_snps %in% cs_variants
)

# --- PIP plot ---
p_pip <- ggplot(df, aes(x = POS / 1e6, y = PIP)) +
  geom_point(alpha = 0.4, size = 1.2, color = 'grey50') +
  geom_point(data = df[df$in_cs, ], aes(x = POS / 1e6, y = PIP),
             color = 'red', size = 2.5) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'orange', alpha = 0.5) +
  geom_hline(yintercept = 0.95, linetype = 'dashed', color = 'red', alpha = 0.5) +
  annotate('text', x = max(df$POS) / 1e6, y = 0.52, label = 'PIP = 0.5',
           hjust = 1, size = 3, color = 'orange') +
  annotate('text', x = max(df$POS) / 1e6, y = 0.97, label = 'PIP = 0.95',
           hjust = 1, size = 3, color = 'red') +
  labs(x = 'Position (Mb)', y = 'Posterior Inclusion Probability') +
  theme_minimal()

ggsave('pip_plot.pdf', p_pip, width = 8, height = 4)

# --- Combined GWAS + PIP ---
p_gwas <- ggplot(df, aes(x = POS / 1e6, y = -log10(P))) +
  geom_point(alpha = 0.4, size = 1, color = 'grey50') +
  geom_point(data = df[df$in_cs, ], color = 'red', size = 2) +
  geom_hline(yintercept = -log10(5e-8), linetype = 'dashed', color = 'red', alpha = 0.3) +
  labs(x = NULL, y = '-log10(P)', title = 'GWAS Association') +
  theme_minimal() + theme(axis.text.x = element_blank())

p_pip2 <- ggplot(df, aes(x = POS / 1e6, y = PIP)) +
  geom_point(alpha = 0.4, size = 1, color = 'grey50') +
  geom_point(data = df[df$in_cs, ], color = 'red', size = 2) +
  labs(x = 'Position (Mb)', y = 'PIP', title = 'Fine-Mapping (credible set in red)') +
  theme_minimal()

p_combined <- p_gwas / p_pip2
ggsave('gwas_pip_combined.pdf', p_combined, width = 8, height = 6)

# --- PIP vs -log10(P) scatter ---
p_scatter <- ggplot(df, aes(x = -log10(P), y = PIP)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_point(data = df[df$in_cs, ], color = 'red', size = 2.5) +
  geom_vline(xintercept = -log10(5e-8), linetype = 'dashed', alpha = 0.3) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.3) +
  labs(x = '-log10(P)', y = 'PIP',
       title = 'Fine-mapping refines GWAS signal',
       subtitle = 'Many significant variants but few high-PIP variants') +
  theme_minimal()

ggsave('pip_vs_pvalue.pdf', p_scatter, width = 6, height = 5)

cat('Plots saved: pip_plot.pdf, gwas_pip_combined.pdf, pip_vs_pvalue.pdf\n')
