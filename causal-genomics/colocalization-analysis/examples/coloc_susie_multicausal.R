# Reference: coloc 5.2.3+, susieR 0.12.35+ | Verify API if version differs
## Multi-causal coloc.susie pipeline with LD-z-score consistency diagnostics
##
## Demonstrates:
##   1. estimate_s_rss diagnostic (LD reference must match z-scores)
##   2. SuSiE per-trait credible sets with L = 10 max signals
##   3. coloc.susie per (CS1, CS2) pair
##   4. Reporting per-pair PP.H4 with credible-set lead SNPs

library(coloc)
library(susieR)

set.seed(42)
n_snps <- 500
positions <- sort(sample(30000000:31000000, n_snps))

## Two independent causal SNPs in GWAS; only one is shared with eQTL
causal1 <- which.min(abs(positions - 30300000))
causal2 <- which.min(abs(positions - 30700000))

## Block-banded LD matrix proxy (in practice: plink --r-phased square on ancestry-matched reference)
ld_matrix <- diag(n_snps)
for (i in 1:n_snps) {
    j_range <- max(1, i - 5):min(n_snps, i + 5)
    for (j in j_range) {
        if (i != j) ld_matrix[i, j] <- 0.8^abs(i - j)
    }
}
ld_matrix <- (ld_matrix + t(ld_matrix)) / 2   # symmetrise

gwas_beta <- rnorm(n_snps, 0, 0.01)
gwas_beta[causal1] <- 0.12
gwas_beta[causal2] <- 0.10
gwas_se <- rep(0.025, n_snps)
gwas_N <- 50000

eqtl_beta <- rnorm(n_snps, 0, 0.02)
eqtl_beta[causal1] <- 0.35   # only causal1 shared with eQTL
eqtl_se <- rep(0.04, n_snps)
eqtl_N <- 500

snp_ids <- paste0('rs', 1:n_snps)

## coloc::runsusie matches the `snp` vector to dimnames(LD); unnamed LD errors out.
dimnames(ld_matrix) <- list(snp_ids, snp_ids)

## Critical diagnostic: lambda > 0.05 means LD does not match z-scores -> abort
z_gwas <- gwas_beta / gwas_se
lam_gwas <- susieR::estimate_s_rss(z=z_gwas, R=ld_matrix, n=gwas_N)
cat(sprintf('GWAS estimate_s_rss lambda = %.4f\n', lam_gwas))
if (lam_gwas > 0.05) stop('GWAS LD-z mismatch; abort or use coloc.abf')

z_eqtl <- eqtl_beta / eqtl_se
lam_eqtl <- susieR::estimate_s_rss(z=z_eqtl, R=ld_matrix, n=eqtl_N)
cat(sprintf('eQTL estimate_s_rss lambda = %.4f\n', lam_eqtl))
if (lam_eqtl > 0.05) stop('eQTL LD-z mismatch; abort or use coloc.abf')

gwas_data <- list(beta=gwas_beta, varbeta=gwas_se^2, snp=snp_ids, position=positions,
                   type='cc', s=0.3, N=gwas_N, LD=ld_matrix)
eqtl_data <- list(beta=eqtl_beta, varbeta=eqtl_se^2, snp=snp_ids, position=positions,
                   type='quant', sdY=1, N=eqtl_N, LD=ld_matrix)

## L = 10 = max number of credible sets to detect; SuSiE returns fewer when fewer are supported
s_gwas <- runsusie(gwas_data, L=10)
s_eqtl <- runsusie(eqtl_data, L=10)

n_cs_gwas <- if (!is.null(summary(s_gwas)$cs)) nrow(summary(s_gwas)$cs) else 0
n_cs_eqtl <- if (!is.null(summary(s_eqtl)$cs)) nrow(summary(s_eqtl)$cs) else 0
cat(sprintf('\nGWAS credible sets: %d | eQTL credible sets: %d\n', n_cs_gwas, n_cs_eqtl))

res_susie <- coloc.susie(s_gwas, s_eqtl)

if (is.null(res_susie$summary) || nrow(res_susie$summary) == 0) {
    cat('\nNo overlapping credible sets between traits -> no colocalization signal\n')
} else {
    cat('\nPer-(CS1, CS2) colocalization PP:\n')
    print(res_susie$summary[, c('idx1', 'idx2', 'hit1', 'hit2',
                                  'PP.H3.abf', 'PP.H4.abf')])
    best <- res_susie$summary[which.max(res_susie$summary$PP.H4.abf), ]
    cat(sprintf('\nBest CS pair: GWAS CS%d (%s) x eQTL CS%d (%s) | PP.H4 = %.3f\n',
                 best$idx1, best$hit1, best$idx2, best$hit2, best$PP.H4.abf))
}
