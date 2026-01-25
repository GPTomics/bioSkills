# Sample Size Estimation Examples
# Demonstrates sample size calculations for various study types

library(ssizeRNA)

# =============================================================================
# RNA-seq Sample Size
# =============================================================================

# Parameters:
# - nGenes: total genes tested
# - pi0: proportion of true null hypotheses (non-DE genes)
# - m: target number of DE genes to detect
# - mu: average read count
# - disp: dispersion parameter (related to CV)
# - fc: fold change
# - fdr: target false discovery rate
# - power: target statistical power

# Typical RNA-seq study
result <- ssizeRNA_single(
  nGenes = 20000,     # Total genes
  pi0 = 0.95,         # 95% genes are not DE
  m = 200,            # Want to detect 200 DE genes
  mu = 10,            # Average count ~10
  disp = 0.1,         # Dispersion (CV^2 for NB)
  fc = 2,             # 2-fold change
  fdr = 0.05,         # 5% FDR
  power = 0.8         # 80% power
)

cat('Required sample size per group:', result$ssize, '\n')

# =============================================================================
# Sensitivity to Effect Size
# =============================================================================

cat('\nSample sizes for different fold changes (80% power, FDR 0.05):\n')
fold_changes <- c(1.5, 2, 2.5, 3, 4)

for (fc in fold_changes) {
  res <- tryCatch({
    ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200, mu = 10,
                    disp = 0.1, fc = fc, fdr = 0.05, power = 0.8)
  }, error = function(e) list(ssize = NA))
  cat(sprintf('  %.1f-fold: n = %s per group\n', fc,
              ifelse(is.na(res$ssize), '>50', res$ssize)))
}

# =============================================================================
# Dispersion Impact (CV equivalent)
# =============================================================================

# dispersion ~ CV^2 for negative binomial
# CV 0.3 -> disp ~ 0.09
# CV 0.4 -> disp ~ 0.16
# CV 0.5 -> disp ~ 0.25

cat('\nSample sizes for different dispersions (2-fold, 80% power):\n')
dispersions <- c(0.05, 0.1, 0.2, 0.3)

for (d in dispersions) {
  cv_equiv <- sqrt(d)
  res <- tryCatch({
    ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200, mu = 10,
                    disp = d, fc = 2, fdr = 0.05, power = 0.8)
  }, error = function(e) list(ssize = NA))
  cat(sprintf('  disp=%.2f (CV~%.2f): n = %s per group\n', d, cv_equiv,
              ifelse(is.na(res$ssize), '>50', res$ssize)))
}

# =============================================================================
# Practical Guidelines by Assay
# =============================================================================

cat('\n=== Recommended Minimum Sample Sizes ===\n')
cat('Based on typical effect sizes and variability:\n\n')

guidelines <- data.frame(
  Assay = c('Bulk RNA-seq (cell lines)', 'Bulk RNA-seq (human)',
            'ATAC-seq', 'ChIP-seq', 'Proteomics (TMT)',
            'Methylation (WGBS)', 'scRNA-seq'),
  MinReplicates = c(3, 5, 2, 2, 4, 4, 3),
  ForSmallEffects = c(5, 10, 4, 4, 8, 8, 6),
  Notes = c('Low CV', 'High donor variability',
            'High technical reproducibility', 'IP efficiency varies',
            'Missing values common', 'CpG site variability',
            '+ 1000 cells/sample minimum')
)

print(guidelines, row.names = FALSE)

# =============================================================================
# Budget Optimization
# =============================================================================

cat('\n=== Budget Optimization ===\n')
cat('With fixed budget, prioritize:\n')
cat('1. Biological replicates > technical replicates\n')
cat('2. More samples > deeper sequencing (after 20M reads for RNA-seq)\n')
cat('3. Balanced designs (equal n per group)\n')
cat('\nAdd 15-20% extra samples for:\n')
cat('- QC failures\n')
cat('- Sample degradation\n')
cat('- Technical issues\n')
