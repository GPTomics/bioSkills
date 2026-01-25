# Multiple Testing Correction Examples
# Demonstrates different correction methods and when to use them

# =============================================================================
# Generate Example P-values
# =============================================================================

set.seed(42)
n_genes <- 10000
n_true_de <- 500  # 5% truly DE

# Simulate p-values
# True nulls: uniform(0,1)
# True positives: beta distribution skewed toward 0
pvalues <- c(
  rbeta(n_true_de, 0.3, 5),      # True DE genes (small p-values)
  runif(n_genes - n_true_de)     # True nulls (uniform)
)
is_true_de <- c(rep(TRUE, n_true_de), rep(FALSE, n_genes - n_true_de))

cat('Simulated', n_genes, 'genes with', n_true_de, 'truly DE\n')
cat('Significant at p < 0.05:', sum(pvalues < 0.05), '\n\n')

# =============================================================================
# Bonferroni Correction
# =============================================================================

# Most conservative - controls family-wise error rate (FWER)
# FWER = P(at least one false positive)
# Use for: small targeted studies, confirmatory analyses

p_bonf <- p.adjust(pvalues, method = 'bonferroni')
sig_bonf <- sum(p_bonf < 0.05)

cat('=== Bonferroni Correction ===\n')
cat('Threshold: p <', 0.05/n_genes, '\n')
cat('Significant genes:', sig_bonf, '\n')
cat('True positives:', sum(p_bonf < 0.05 & is_true_de), '\n')
cat('False positives:', sum(p_bonf < 0.05 & !is_true_de), '\n\n')

# =============================================================================
# Benjamini-Hochberg FDR
# =============================================================================

# Standard for genomics - controls false discovery rate
# FDR = E[false positives / total positives]
# More powerful than Bonferroni when many true positives exist

p_bh <- p.adjust(pvalues, method = 'BH')
sig_bh <- sum(p_bh < 0.05)

cat('=== Benjamini-Hochberg FDR ===\n')
cat('Significant at FDR < 0.05:', sig_bh, '\n')
cat('True positives:', sum(p_bh < 0.05 & is_true_de), '\n')
cat('False positives:', sum(p_bh < 0.05 & !is_true_de), '\n')
cat('Observed FDR:', round(sum(p_bh < 0.05 & !is_true_de) / sig_bh, 3), '\n\n')

# =============================================================================
# q-value Method
# =============================================================================

# More powerful than BH by estimating pi0 (proportion true nulls)
# Directly estimates FDR for each feature

library(qvalue)
qobj <- qvalue(pvalues)
qvalues <- qobj$qvalues
sig_qval <- sum(qvalues < 0.05)

cat('=== q-value Method ===\n')
cat('Estimated pi0 (true null proportion):', round(qobj$pi0, 3), '\n')
cat('Significant at q < 0.05:', sig_qval, '\n')
cat('True positives:', sum(qvalues < 0.05 & is_true_de), '\n')
cat('False positives:', sum(qvalues < 0.05 & !is_true_de), '\n\n')

# =============================================================================
# Method Comparison
# =============================================================================

cat('=== Method Comparison ===\n')
comparison <- data.frame(
  Method = c('None (p < 0.05)', 'Bonferroni', 'BH FDR', 'q-value'),
  Significant = c(sum(pvalues < 0.05), sig_bonf, sig_bh, sig_qval),
  TruePos = c(sum(pvalues < 0.05 & is_true_de),
              sum(p_bonf < 0.05 & is_true_de),
              sum(p_bh < 0.05 & is_true_de),
              sum(qvalues < 0.05 & is_true_de)),
  FalsePos = c(sum(pvalues < 0.05 & !is_true_de),
               sum(p_bonf < 0.05 & !is_true_de),
               sum(p_bh < 0.05 & !is_true_de),
               sum(qvalues < 0.05 & !is_true_de))
)
comparison$Sensitivity <- round(comparison$TruePos / n_true_de, 3)
comparison$FDR <- round(comparison$FalsePos / comparison$Significant, 3)

print(comparison, row.names = FALSE)

# =============================================================================
# Method Selection Guide
# =============================================================================

cat('\n=== When to Use Each Method ===\n')
cat('Bonferroni:\n')
cat('  - Small, targeted gene panels (<100 genes)\n')
cat('  - Confirmatory/validation studies\n')
cat('  - When ANY false positive is unacceptable\n\n')

cat('Benjamini-Hochberg:\n')
cat('  - Standard genome-wide DE analysis\n')
cat('  - Exploratory studies\n')
cat('  - When some false positives are acceptable\n\n')

cat('q-value:\n')
cat('  - Large-scale studies with many true positives\n')
cat('  - Maximum power needed\n')
cat('  - When pi0 < 0.9 (many true effects expected)\n\n')

cat('GWAS:\n')
cat('  - Use genome-wide threshold p < 5e-8\n')
cat('  - Bonferroni for ~1 million independent tests\n')
