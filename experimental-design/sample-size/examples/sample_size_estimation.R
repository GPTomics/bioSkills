# Sample size for genomics: FDR-aware NB sizing, pilot dispersions, verifying a fixed n
# Reference: ssizeRNA 1.3+, DESeq2 1.42+ | Verify API if version differs
#
# Demonstrates: ssizeRNA_single/_vary (note: `m` is the pseudo sample size for SIMULATION,
# default 200 -- NOT the number of DE genes or n per group), estimating dispersion from a
# pilot, and check.power() to confirm a budget-fixed n meets the target FDR.

suppressPackageStartupMessages(library(ssizeRNA))
set.seed(20260528)

# ---------------------------------------------------------------------------
# 1. FDR-aware sample size (single mean/dispersion for all genes)
# ---------------------------------------------------------------------------
# nGenes: total genes; pi0: proportion NON-DE; m: pseudo sample size for simulation (default 200);
# mu: mean control count; disp: NB dispersion; fc: fold change; fdr/power: targets.
res <- ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200,
                       mu = 10, disp = 0.2, fc = 1.5, fdr = 0.05, power = 0.80, maxN = 30)
cat('Minimum n per group (1.5-fold, disp=0.2, FDR 0.05, 80% power):', res$ssize, '\n')

# Sensitivity to fold change and dispersion (the two dominant levers)
for (fc in c(1.5, 2, 3)) {
  r <- ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200, mu = 10, disp = 0.2,
                       fc = fc, fdr = 0.05, power = 0.80, maxN = 30)
  cat(sprintf('fc=%.1f -> n=%s per group\n', fc, r$ssize))
}

# ---------------------------------------------------------------------------
# 2. Estimate dispersion from a pilot (defensible input) and size with a vector of dispersions
# ---------------------------------------------------------------------------
if (requireNamespace('DESeq2', quietly = TRUE)) {
  library(DESeq2)
  # ... fit DESeq() on pilot counts, then:
  # disp_vec <- dispersions(dds); mu_vec <- rowMeans(counts(dds, normalized = TRUE))
  # ssizeRNA_vary(nGenes = length(disp_vec), pi0 = 0.95, mu = mu_vec, disp = disp_vec,
  #               fc = 1.5, fdr = 0.05, power = 0.80, maxN = 30)$ssize
  cat('\nUse DESeq2::dispersions(dds) from a pilot to feed ssizeRNA_vary (mu/disp as vectors).\n')
}

# ---------------------------------------------------------------------------
# 3. Verify a budget-fixed n: average power and TRUE realized FDR (here m = n per group)
# ---------------------------------------------------------------------------
cp <- check.power(nGenes = 20000, pi0 = 0.95, m = 6, mu = 10, disp = 0.2,
                  fc = 1.5, fdr = 0.05, sims = 50)
cat(sprintf('\nAt n=6/group: BH average power = %.2f, true FDR = %.3f\n',
            cp$pow_bh_ave, cp$fdr_bh_ave))

# ---------------------------------------------------------------------------
# 4. Assay floors (NOT targets) -- defensible only after pilot/literature dispersion supports them
# ---------------------------------------------------------------------------
floors <- data.frame(
  assay = c('Bulk RNA-seq', 'scRNA-seq (population DE)', 'ATAC-seq', 'ChIP-seq',
            'Proteomics (DIA/TMT)', 'Methylation (WGBS)'),
  min_replicates = c(3, 3, 2, 2, 3, 4),                 # floors under low dispersion + large effects
  for_small_effects = c('6-12', '6+ donors', '4-6', '3-4', '6-10', '8-12'),
  note = c('Schurch 2016: >=6 recovers most true DE', 'donors, not cells, set power (Squair 2021)',
           'library complexity floor', 'IDR reproducibility (ENCODE)',
           'higher missingness (MNAR)', 'high per-CpG variance'))
print(floors, row.names = FALSE)
cat('\n# Count biological replicates, not measurements; add 10-20% for sample failures.\n')
