# Reference: susieR 0.12+, coloc 5.2+, PolyFun (head 2024) | Verify API if version differs
##
## Fine-mapping a GWAS locus with susie_rss including the LD diagnostic block,
## optional PolyFun functional priors, credible set extraction with purity
## filter, and reporting conventions.

library(susieR)
library(dplyr)

# ----- Inputs ---------------------------------------------------------------
# Replace these paths with real data; the simulated block below is for demo.
gwas_path <- 'locus_sumstats.tsv'
ld_path   <- 'locus_ld_matrix.tsv'
n_gwas    <- 500000
locus_chr <- 6
locus_start <- 30000000
locus_end   <- 31000000

# ----- Simulated locus for self-contained demo ------------------------------
# Comment this block out when feeding real data.
set.seed(42)
n_snps <- 500
n_causal <- 2
true_causal <- c(120, 360)

block_ld <- diag(n_snps)
for (i in seq_len(n_snps)) {
    for (j in max(1, i - 15):min(n_snps, i + 15)) {
        if (i != j) block_ld[i, j] <- 0.85 ^ abs(i - j)
    }
}

z_sim <- rnorm(n_snps, 0, 1)
for (c in true_causal) {
    window <- max(1, c - 6):min(n_snps, c + 6)
    bump <- 6.5 * exp(-0.18 * (window - c) ^ 2)
    z_sim[window] <- z_sim[window] + bump
}

gwas_df <- data.frame(
    SNP = sprintf('rs%07d', seq_len(n_snps)),
    CHR = locus_chr,
    POS = locus_start + seq_len(n_snps) * 2000,
    Z   = z_sim,
    P   = 2 * pnorm(-abs(z_sim))
)
ld_matrix <- block_ld

# ----- Real-data loader (uncomment when ready) -------------------------------
# gwas_df <- read.table(gwas_path, header = TRUE, sep = '\t')
# stopifnot(all(c('SNP', 'BETA', 'SE') %in% names(gwas_df)))
# gwas_df$Z <- gwas_df$BETA / gwas_df$SE
# ld_matrix <- as.matrix(read.table(ld_path))
# stopifnot(nrow(ld_matrix) == nrow(gwas_df))

# ----- LD positive semi-definite check -------------------------------------
# Negative eigenvalues from finite-precision storage break susie_rss internals.
eig <- eigen(ld_matrix, only.values = TRUE)$values
if (min(eig) < -1e-6) {
    cat(sprintf('LD matrix has min eigenvalue %.2e; adding ridge\n', min(eig)))
    ridge <- abs(min(eig)) + 1e-4
    ld_matrix <- ld_matrix + diag(ridge, nrow(ld_matrix))
}

# ----- Mandatory LD diagnostic block ---------------------------------------
# estimate_s_rss returns the inferred LD inconsistency scale.
# Source: susieR vignette "Diagnostic for summary statistic"; Zou 2022 PLoS Genet.
# Threshold convention: < 0.05 acceptable; 0.05-0.10 marginal; > 0.10 refit.
s_hat <- estimate_s_rss(z = gwas_df$Z, R = ld_matrix, n = n_gwas)
cat(sprintf('estimate_s_rss lambda = %.4f\n', s_hat))
if (s_hat > 0.10) {
    warning('Lambda > 0.10: LD reference likely mismatches the GWAS sample. ',
            'Consider in-sample LD or ancestry-stratified reference.')
}

# kriging_rss flags per-SNP inconsistency (typically strand flips / coding mismatches)
cond_z <- kriging_rss(z = gwas_df$Z, R = ld_matrix, n = n_gwas)
flag_idx <- which(abs(cond_z$conditional_dist$z_std_diff) > 3)
cat(sprintf('kriging_rss flagged %d SNPs with |z_obs - z_exp| > 3\n', length(flag_idx)))

# ----- Optional PolyFun functional priors ----------------------------------
# Skip this section if no PolyFun output exists; uniform priors are the default.
use_polyfun <- FALSE
if (use_polyfun) {
    priors <- read.table('polyfun_priors.6.snpvar.gz', header = TRUE)
    priors <- priors[match(gwas_df$SNP, priors$SNP), ]
    prior_w <- priors$SNPVAR / sum(priors$SNPVAR, na.rm = TRUE)
    # NOTE: per-SNP causal probability goes to prior_weights, NOT prior_variance.
} else {
    prior_w <- NULL
}

# ----- Fit susie_rss --------------------------------------------------------
# L=10 is the default; raise to 20-30 for HLA or loci expected to have >5 signals.
fit <- susie_rss(
    z = gwas_df$Z,
    R = ld_matrix,
    n = n_gwas,
    L = 10,
    prior_weights = prior_w,
    estimate_residual_variance = TRUE
)

if (!isTRUE(fit$converged)) {
    warning('SuSiE IBSS did not converge; increase max_iter or check LD')
}

cat(sprintf('Effective L used (credible sets returned): %d / requested %d\n',
            length(fit$sets$cs),
            10))

# ----- Extract credible sets with purity filter ----------------------------
# Purity = min absolute correlation among variants in the set.
# Convention: purity >= 0.5 (r2 >= 0.25) for a well-resolved set.
cs_list <- fit$sets$cs
purity_mat <- fit$sets$purity

gwas_df$PIP <- fit$pip
gwas_df$in_credible_set <- NA_integer_

valid_sets <- list()
for (i in seq_along(cs_list)) {
    snp_idx <- cs_list[[i]]
    purity_min <- purity_mat[i, 'min.abs.corr']
    top_snp <- snp_idx[which.max(fit$pip[snp_idx])]
    top_pip <- fit$pip[top_snp]

    cat(sprintf('Credible set %d: size=%d, purity=%.3f, top SNP=%s (PIP=%.3f)\n',
                i, length(snp_idx), purity_min, gwas_df$SNP[top_snp], top_pip))

    if (purity_min >= 0.5) {
        valid_sets[[length(valid_sets) + 1]] <- list(
            set_id = i, snps = snp_idx, purity = purity_min,
            top_snp = top_snp, top_pip = top_pip
        )
        gwas_df$in_credible_set[snp_idx] <- i
    } else {
        cat(sprintf('  -> dropped (purity < 0.5; LD-confounded)\n'))
    }
}

cat(sprintf('\nValid credible sets after purity filter: %d\n', length(valid_sets)))

# ----- Report top PIP variants ---------------------------------------------
top_pip_df <- gwas_df[order(-gwas_df$PIP), ][1:15, c('SNP', 'CHR', 'POS', 'Z', 'P', 'PIP', 'in_credible_set')]
cat('\nTop 15 variants by PIP:\n')
print(top_pip_df, row.names = FALSE)

cat(sprintf('\nVariants with PIP > 0.95: %d\n', sum(fit$pip > 0.95)))
cat(sprintf('Variants with PIP > 0.50: %d\n', sum(fit$pip > 0.50)))

# ----- Save outputs ---------------------------------------------------------
write.table(gwas_df, 'finemap_pips.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

credset_out <- do.call(rbind, lapply(valid_sets, function(v) {
    data.frame(
        credible_set = v$set_id,
        size = length(v$snps),
        purity = v$purity,
        top_snp = gwas_df$SNP[v$top_snp],
        top_pip = v$top_pip,
        snps = paste(gwas_df$SNP[v$snps], collapse = ',')
    )
}))
write.table(credset_out, 'finemap_credible_sets.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

cat('\nWrote finemap_pips.tsv and finemap_credible_sets.tsv\n')
