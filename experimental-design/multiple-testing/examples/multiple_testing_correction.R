# Multiple testing: FDR vs FWER, dependence (BH vs BY), q-value/pi0, local FDR, IHW
# Reference: R stats (base), qvalue 2.34+, IHW 1.30+ | Verify API if version differs
#
# Demonstrates the error-rate choices that matter in genomics discovery and the levers
# (pi0, covariate weighting) that buy back power. Python statsmodels note: the DEFAULT
# method is 'hs' (Holm-Sidak), NOT BH -- always pass method= explicitly.

set.seed(20260528)
n_genes <- 10000; n_de <- 500
pvalues <- c(rbeta(n_de, 0.3, 5), runif(n_genes - n_de))   # 5% true DE (small p), rest null
is_de  <- c(rep(TRUE, n_de), rep(FALSE, n_genes - n_de))

# ---------------------------------------------------------------------------
# 1. FWER (Bonferroni/Holm) vs FDR (BH) vs FDR-under-dependence (BY)
# ---------------------------------------------------------------------------
methods <- c(none = NA, bonferroni = 'bonferroni', holm = 'holm', BH = 'BH', BY = 'BY')
tab <- data.frame(method = names(methods), significant = NA_integer_, true_pos = NA_integer_,
                  false_pos = NA_integer_)
for (i in seq_along(methods)) {
  padj <- if (is.na(methods[i])) pvalues else p.adjust(pvalues, method = methods[i])
  sig <- padj < 0.05
  tab[i, c('significant', 'true_pos', 'false_pos')] <- c(sum(sig), sum(sig & is_de), sum(sig & !is_de))
}
tab$realized_fdr <- round(tab$false_pos / pmax(tab$significant, 1), 3)
print(tab, row.names = FALSE)
# BH is the discovery default (valid under independence/PRDS); BY is conservative but valid
# under arbitrary/negative dependence.

# ---------------------------------------------------------------------------
# 2. q-value: estimate pi0 (true-null proportion) for more power; local FDR per feature
# ---------------------------------------------------------------------------
if (requireNamespace('qvalue', quietly = TRUE)) {
  library(qvalue)
  qobj <- qvalue(pvalues)
  cat(sprintf('\nq-value: pi0 = %.3f, discoveries at q<0.05 = %d\n',
              qobj$pi0, sum(qobj$qvalues < 0.05)))
  # qobj$lfdr is the local FDR: posterior P(null | statistic) for EACH feature.
}

# ---------------------------------------------------------------------------
# 3. IHW: weight hypotheses by an independent covariate (must be null-independent)
# ---------------------------------------------------------------------------
if (requireNamespace('IHW', quietly = TRUE)) {
  library(IHW)
  mean_expr <- rgamma(n_genes, shape = 2, rate = 0.5)        # covariate independent of null p
  res <- ihw(pvalues, mean_expr, alpha = 0.05)
  cat(sprintf('IHW discoveries at FDR 0.05 = %d (vs BH = %d)\n',
              rejections(res), sum(p.adjust(pvalues, 'BH') < 0.05)))
}

# ---------------------------------------------------------------------------
# 4. Python equivalent (mind the default method)
# ---------------------------------------------------------------------------
# from statsmodels.stats.multitest import multipletests
# rej, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')   # BH -- NOT the default
# rej, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_by')   # BY under dependence

# GWAS: genome-wide significance ~5e-8 (Dudbridge & Gusnanto 2008 derived ~7.2e-8);
# the GWAS test machinery lives in population-genetics/association-testing.
