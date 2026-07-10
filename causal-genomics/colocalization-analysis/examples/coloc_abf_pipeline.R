# Reference: coloc 5.2.3+, susieR 0.12.35+ | Verify API if version differs
## Single-causal coloc.abf pipeline with harmonisation, p12 sensitivity, and reporting
##
## Demonstrates:
##   1. Allele harmonisation between two summary statistics
##   2. Locus extraction (1 Mb window)
##   3. coloc.abf with conservative p12 = 5e-6
##   4. p12 sensitivity analysis (Wallace 2020 recommendation)
##   5. Reporting PP.H4 alongside PP.H3 with the safe-range over p12

library(coloc)

harmonise_summary_stats <- function(df1, df2) {
    m <- merge(df1, df2, by='SNP', suffixes=c('.1','.2'))
    same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
    flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2
    palindromic <- (m$A1.1 %in% c('A','T') & m$A2.1 %in% c('A','T')) |
                   (m$A1.1 %in% c('C','G') & m$A2.1 %in% c('C','G'))
    high_maf_palin <- palindromic & pmin(m$MAF.1, 1 - m$MAF.1) > 0.42
    m$BETA.2[flip] <- -m$BETA.2[flip]
    keep <- (same | flip) & !high_maf_palin
    m[keep, ]
}

extract_locus <- function(sumstats, chr, lead_pos, window=500000) {
    locus <- sumstats[sumstats$CHR == chr &
                       sumstats$POS >= (lead_pos - window) &
                       sumstats$POS <= (lead_pos + window), ]
    locus[order(locus$POS), ]
}

set.seed(42)
n_snps <- 1000
positions <- sort(sample(30000000:31000000, n_snps))
causal_idx <- which.min(abs(positions - 30500000))

gwas_beta <- rnorm(n_snps, 0, 0.02); gwas_beta[causal_idx] <- 0.15
gwas_se <- rep(0.03, n_snps)
gwas_df <- data.frame(SNP=paste0('rs', 1:n_snps), CHR=6, POS=positions,
                       A1='A', A2='G', MAF=runif(n_snps, 0.05, 0.5),
                       BETA=gwas_beta, SE=gwas_se)

eqtl_beta <- rnorm(n_snps, 0, 0.03); eqtl_beta[causal_idx] <- 0.4
eqtl_se <- rep(0.05, n_snps)
eqtl_df <- data.frame(SNP=paste0('rs', 1:n_snps), CHR=6, POS=positions,
                       A1='A', A2='G', MAF=runif(n_snps, 0.05, 0.5),
                       BETA=eqtl_beta, SE=eqtl_se)

harm <- harmonise_summary_stats(gwas_df, eqtl_df)
cat('Harmonised SNPs:', nrow(harm), 'of', nrow(gwas_df), '\n')

gwas_input <- list(beta=harm$BETA.1, varbeta=harm$SE.1^2,
                    snp=harm$SNP, position=harm$POS.1,
                    type='cc', s=0.30, N=50000)

eqtl_input <- list(beta=harm$BETA.2, varbeta=harm$SE.2^2,
                    snp=harm$SNP, position=harm$POS.1,
                    type='quant', sdY=1, N=500)

## p12 = 5e-6 is the conservative cis-eQTL default (lower than coloc default 1e-5)
res <- coloc.abf(dataset1=gwas_input, dataset2=eqtl_input,
                  p1=1e-4, p2=1e-4, p12=5e-6)

cat('\nPosterior probabilities:\n')
print(round(res$summary, 4))

pp4 <- res$summary['PP.H4.abf']
pp3 <- res$summary['PP.H3.abf']
cat(sprintf('\nPP.H4 = %.3f | PP.H3 = %.3f\n', pp4, pp3))

if (pp4 >= 0.9) {
    cat('Stringent threshold passed (>= 0.90)\n')
} else if (pp4 >= 0.80) {
    cat('Published-grade threshold passed (>= 0.80)\n')
} else if (pp4 >= 0.75) {
    cat('Open Targets / screening threshold passed (>= 0.75)\n')
} else if (pp4 >= 0.50) {
    cat('Suggestive only; recommend coloc.susie or larger N\n')
} else {
    cat('Below threshold; report H0/H1/H2/H3 dominance instead\n')
}

## p12 sensitivity is non-negotiable per Wallace 2020
## Identifies the lowest p12 at which PP.H4 still passes 0.75
cat('\nRunning p12 sensitivity over [1e-8, 1e-4]\n')
sens <- coloc::sensitivity(res, rule='H4 > 0.75')

top_snps <- res$results[order(-res$results$SNP.PP.H4), ][1:5, ]
cat('\nTop 5 SNPs by per-SNP PP.H4:\n')
print(top_snps[, c('snp', 'SNP.PP.H4')])
