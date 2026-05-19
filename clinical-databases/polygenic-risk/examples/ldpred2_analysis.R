# Reference: bigsnpr 1.12+, R 4.3+ | Verify snp_ldpred2_auto signature if version differs
# LDpred2-auto PRS construction with sample-overlap detection, LD-mismatch QC,
# and ancestry-conditional Z normalization.
# Pin allow_jump_sign=FALSE and shrink_corr=0.95 per Prive 2022 misspecification paper.

library(bigsnpr)
library(data.table)

run_ldpred2_auto <- function(target_rds, sumstats_path, ld_ref_path, n_cores = 8,
                              hapmap3_only = TRUE) {
    # 1. Load target genotypes (bigSNP object from snp_readBed)
    obj <- snp_attach(target_rds)
    G <- obj$genotypes
    map <- obj$map
    setnames(map, c("chr", "rsid", "genetic.dist", "pos", "a1", "a0"))

    # 2. Load GWAS summary statistics
    sumstats <- fread(sumstats_path)
    # Required columns: chr, pos, a0 (ref/non-effect), a1 (effect), beta, beta_se, n_eff, p

    # 3. Restrict to HapMap3 if using PRS-CS/LDpred2 default reference (recommended for stability)
    if (hapmap3_only) {
        hm3 <- fread("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3.csv")
        sumstats <- sumstats[rsid %in% hm3$rsid]
    }

    # 4. Match variants strand-aware (handles A/T C/G ambiguity by allele frequency)
    df_beta <- snp_match(sumstats, map, strand_flip = TRUE, join_by_pos = FALSE)
    message("Matched ", nrow(df_beta), " variants from ", nrow(sumstats), " sumstats")

    # 5. Load LD correlation matrix
    # Prefer UKB-LD (~40k samples) over 1KG-EUR (~489); larger ref reduces s misspecification
    corr_full <- readRDS(ld_ref_path)  # output of snp_cor or precomputed UKB LD
    corr <- corr_full[df_beta$`_NUM_ID_`, df_beta$`_NUM_ID_`]

    # 6. LDSC heritability + LD-mismatch diagnostic
    ldsc_res <- snp_ldsc2(corr, df_beta)
    h2_est <- ldsc_res[["h2"]]
    s_est <- ldsc_res[["int"]]  # LD-mismatch / intercept; large => sample overlap
    message(sprintf("LDSC h2 = %.3f; intercept (s) = %.3f", h2_est, s_est))
    if (s_est > 0.05) {
        warning("LDSC intercept > 0.05 -- possible sample overlap or LD mismatch. ",
                "Run EraSOR or bivariate LDSC for confirmation.")
    }

    # 7. LDpred2-auto with multiple chains
    # CRITICAL: allow_jump_sign = FALSE and shrink_corr = 0.95 per Prive 2022
    multi_auto <- snp_ldpred2_auto(
        corr, df_beta,
        h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.2, 30),
        burn_in = 500, num_iter = 200,
        allow_jump_sign = FALSE,   # required
        shrink_corr = 0.95,        # required
        ncores = n_cores
    )

    # 8. Filter divergent chains (Prive 2022 quality filter)
    range_auto <- sapply(multi_auto, function(x) diff(range(x$corr_est)))
    keep <- range_auto > (0.95 * quantile(range_auto, 0.95, na.rm = TRUE))
    message("Keeping ", sum(keep), " / ", length(keep), " chains after convergence filter")

    beta_auto <- sapply(multi_auto[keep], function(x) x$beta_est)
    beta_final <- rowMeans(beta_auto)

    # 9. Score
    prs <- big_prodMat(G, beta_final, ind.col = df_beta$`_NUM_ID_`)

    list(prs = prs, beta = beta_final, h2 = h2_est, s = s_est, n_chains = sum(keep))
}


ancestry_conditional_z <- function(prs, pcs, n_pcs = 10) {
    # Recalibrate PRS by removing ancestry effects (Sun 2021; Ding 2023 continuous-ancestry)
    # IMPORTANT: PCs must be computed in the TEST cohort, NOT discovery cohort
    pc_design <- cbind(intercept = 1, pcs[, 1:n_pcs])
    mean_fit <- lm.fit(pc_design, prs)
    residuals <- prs - mean_fit$fitted.values

    # Variance also varies by PC (Sun 2021)
    log_var_fit <- lm.fit(pc_design, log(residuals^2 + 1e-12))
    sd_est <- sqrt(exp(log_var_fit$fitted.values))

    z <- residuals / sd_est
    list(z = z, percentile = pnorm(z) * 100)
}


sample_overlap_check_command <- function() {
    # EraSOR / bivariate LDSC intercept for sample-overlap detection
    cat("# Bivariate LDSC intercept (required before any PRS evaluation):\n")
    cat("ldsc.py \\\n")
    cat("    --rg target_sumstats.sumstats.gz,discovery_sumstats.sumstats.gz \\\n")
    cat("    --ref-ld-chr eur_w_ld_chr/ \\\n")
    cat("    --w-ld-chr eur_w_ld_chr/ \\\n")
    cat("    --out overlap_check\n\n")
    cat("# In the *.log, inspect 'gcov_int' (genetic-covariance intercept).\n")
    cat("# |gcov_int| > 0.05 with target n >= 1000 indicates substantial overlap.\n")
    cat("# Reject PRS evaluation if confirmed; use disjoint test set.\n")
}


if (sys.nframe() == 0) {
    message("LDpred2-auto workflow (production version requires real genotype + sumstats):")
    message("  result <- run_ldpred2_auto(")
    message("      target_rds = 'target.rds',")
    message("      sumstats_path = 'gwas.sumstats',")
    message("      ld_ref_path = 'ukb_ld_eur.rds',")
    message("      n_cores = 8)")
    message("  # Apply ancestry-conditional Z normalization on result$prs using test-cohort PCs")
    message("  # Report HR per SD + absolute-risk integration; cite Hingorani 2023 BMJ Med caveats")
    cat("\n--- Sample overlap detection (run BEFORE evaluation) ---\n")
    sample_overlap_check_command()
}
