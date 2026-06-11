# Reference: MultiAssayExperiment 1.36+ | Verify API if version differs
# Pre-integration design diagnostic: assemble a MultiAssayExperiment, then read off the
# three decisions that determine whether and how to integrate - sample correspondence
# (paired/mosaic/horizontal), the n<<p regime, and the per-block variance imbalance.

library(MultiAssayExperiment)

set.seed(1)

n_total <- 60                                    # patients in the study
samples <- paste0('S', seq_len(n_total))

make_block <- function(ids, n_feat, prefix) {
    m <- matrix(rnorm(length(ids) * n_feat), nrow=n_feat, dimnames=list(paste0(prefix, seq_len(n_feat)), ids))
    m
}

rna    <- make_block(samples, 2000, 'gene')              # all 60 patients
prot   <- make_block(samples[1:48], 300, 'prot')         # mosaic: 12 patients lack proteomics
methyl <- make_block(samples[6:60], 5000, 'cg')          # 5000 CpGs, 55 patients; 5 lack methylation

clinical <- DataFrame(row.names=samples, group=factor(rep(c('responder', 'nonresponder'), length.out=n_total)))

mae <- MultiAssayExperiment(experiments=ExperimentList(rna=rna, prot=prot, methyl=methyl), colData=clinical)

cat('--- correspondence ---\n')
n_samp <- nrow(colData(mae))
per_omic_n <- vapply(experiments(mae), ncol, integer(1))
complete_n <- sum(complete.cases(mae))
mosaic_frac <- 1 - complete_n / n_samp
print(per_omic_n)
cat('samples with every omic:', complete_n, 'of', n_samp, '| mosaic fraction:', round(mosaic_frac, 2), '\n')

cat('\n--- n vs p ---\n')
for (a in names(experiments(mae))) {
    cat(sprintf('%-7s n=%d p=%d  n/p=%.4f\n', a, ncol(mae[[a]]), nrow(mae[[a]]), ncol(mae[[a]]) / nrow(mae[[a]])))
}

cat('\n--- variance imbalance (after per-feature scaling) ---\n')
block_var <- vapply(experiments(mae), function(x) sum(apply(t(scale(t(as.matrix(x)))), 1, var, na.rm=TRUE)), numeric(1))
share <- block_var / sum(block_var)                       # each block's share of stacked variance
print(round(share, 3))

dominance_flag <- 0.80                                     # one view > 80% of stacked variance hijacks shared factors (Tini 2019)
recommend <- if (mosaic_frac > 0.10) {
    'mosaic cohort -> MOFA2 (models missing-view samples) rather than complete-case intersection'
} else if (max(share) > dominance_flag) {
    'one omic dominates variance -> equalize blocks (MFA / per-block keepX / z-score) or use SNF before trusting shared factors'
} else {
    'paired and balanced -> MOFA2 for an unsupervised map; escalate to DIABLO (supervised) or SNF (subtypes) per question'
}
cat('\nrecommendation:', recommend, '\n')
