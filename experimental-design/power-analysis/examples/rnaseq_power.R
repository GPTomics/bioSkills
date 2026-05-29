# Power analysis for RNA-seq: closed-form vs simulation, per-gene power
# Reference: RNASeqPower 1.42+, PROPER 1.34+ | Verify API if version differs
#
# Demonstrates: closed-form NB power (both directions), why a single CV is only an
# approximation, simulation-based marginal power at a target FDR, and the depth-vs-replicate
# tradeoff. Observed/post-hoc power is deliberately NOT computed -- it is uninformative.

suppressPackageStartupMessages(library(RNASeqPower))

# ---------------------------------------------------------------------------
# 1. Closed-form NB power -- solves for whichever of n / power is omitted
# ---------------------------------------------------------------------------
# depth = reads/gene; cv = biological coefficient of variation; effect = fold change
cat('Power, n=3, 2-fold, CV=0.4:', round(rnapower(depth = 20, n = 3, cv = 0.4, effect = 2, alpha = 0.05), 3), '\n')
cat('Power, n=6, 2-fold, CV=0.4:', round(rnapower(depth = 20, n = 6, cv = 0.4, effect = 2, alpha = 0.05), 3), '\n')
cat('n for 80% power, 2-fold, CV=0.4:',
    ceiling(rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.80)), 'per group\n')

# Effect / CV sensitivity: a SINGLE CV gives one curve, but real dispersion varies with mean.
for (cv in c(0.2, 0.4)) {
  n <- ceiling(rnapower(depth = 20, cv = cv, effect = 1.5, alpha = 0.05, power = 0.80))
  cat(sprintf('CV %.1f -> n=%d per group for 1.5-fold at 80%% power\n', cv, n))
}

# ---------------------------------------------------------------------------
# 2. Simulation-based MARGINAL power at a target FDR (the honest default)
# ---------------------------------------------------------------------------
# Closed-form uses one CV; PROPER simulates from an empirical mean-dispersion trend and
# reports average power across the expression distribution at a controlled FDR.
if (requireNamespace('PROPER', quietly = TRUE)) {
  library(PROPER)
  sim_opts <- RNAseq.SimOptions.2grp(ngenes = 20000, p.DE = 0.05,
                                     lOD = 'cheung', lBaselineExpr = 'cheung')
  sims <- runSims(Nreps = c(3, 5, 8, 12), sim.opts = sim_opts, nsims = 20, DEmethod = 'edgeR')
  powr <- comparePower(sims, alpha.type = 'fdr', alpha.nominal = 0.05,
                       stratify.by = 'expr', delta = log(1.5))   # delta is natural-log lfc in PROPER (not log2)
  print(summaryPower(powr))   # marginal power per replicate number at FDR 0.05
} else {
  cat('\nInstall PROPER for simulation-based marginal power (the reported figure).\n')
}

# ---------------------------------------------------------------------------
# 3. Depth vs replicates
# ---------------------------------------------------------------------------
# Past ~10-20M mapped reads, deeper sequencing adds little; replicates keep helping
# (Liu, Zhou & White 2014, Bioinformatics 30:301).
for (d in c(10, 20, 50, 100)) {
  cat(sprintf('Depth %3d, n=4: power = %.3f\n', d,
              rnapower(depth = d, n = 4, cv = 0.4, effect = 2, alpha = 0.05)))
}
cat('# Adding a replicate beats doubling depth once depth is adequate.\n')
