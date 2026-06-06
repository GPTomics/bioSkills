# Reference: R 4.3+ (base stats); production pmp 1.14+, statTarget 1.30+, imputeLCMD 2.1+ | Verify API if version differs
# End-to-end untargeted-metabolomics preprocessing on synthetic data, demonstrating the
# load-bearing claims: D-ratio filtering (technical vs biological variance), QC-anchored drift
# correction validated on HELD-OUT QCs, PQN dilution normalization, and mechanism-aware imputation
# (left-censored MNAR vs sporadic MAR). Self-contained in base R so it runs without Bioconductor;
# the SKILL.md shows the production pmp/statTarget/imputeLCMD equivalents.

set.seed(1)

n_features <- 200
n_bio <- 60
n_qc <- 12

# Injection schedule: QCs interspersed every ~5 biological samples, bracketing both ends.
order_bio <- sort(sample(seq_len(n_bio + n_qc), n_bio))
order_qc <- setdiff(seq_len(n_bio + n_qc), order_bio)
n_total <- n_bio + n_qc

group <- rep(c('control', 'case'), length.out = n_bio)
true_abundance <- 2^rnorm(n_features, mean = 10, sd = 2)

# Biological effect: 20 of 200 features modestly up in cases (the signal to preserve).
# Kept to 10% of features and a 1.5x fold so the PQN median-quotient recovers dilution, not biology.
effect <- rep(1, n_features)
effect[1:20] <- 1.5

# Per-feature monotonic drift vs injection order, strong enough that correction genuinely helps.
drift_slope <- rnorm(n_features, mean = 0, sd = 0.02)

# Per-sample dilution (random, unrelated to group) -- the quantity PQN should recover.
# QC injections share one homogeneous pool, so their dilution is fixed at 1.
dilution <- setNames(2^rnorm(n_total, sd = 0.4), seq_len(n_total))
dilution[as.character(order_qc)] <- 1

build_sample <- function(inj_order, fold) {
    drift <- 1 + drift_slope * inj_order
    noise <- 2^rnorm(n_features, sd = 0.1)
    true_abundance * fold * drift * noise * dilution[as.character(inj_order)]
}

mat <- matrix(NA_real_, nrow = n_total, ncol = n_features)
is_qc <- logical(n_total)
sample_group <- character(n_total)
inj <- integer(n_total)

for (i in seq_len(n_bio)) {
    o <- order_bio[i]
    fold <- ifelse(group[i] == 'case', effect, 1)
    mat[o, ] <- build_sample(o, fold)
    sample_group[o] <- group[i]
    inj[o] <- o
}
for (o in order_qc) {
    mat[o, ] <- build_sample(o, 1)
    is_qc[o] <- TRUE
    sample_group[o] <- 'QC'
    inj[o] <- o
}

# Inject left-censored missingness below an LOD threshold. The 15 lowest-abundance features are
# censored hard (>50% missing) so the detection filter drops them; the next 15 stay recoverable.
lod_features <- order(true_abundance)[1:30]
for (i in seq_along(lod_features)) {
    f <- lod_features[i]
    censor_q <- if (i <= 15) 0.65 else 0.35
    lod <- quantile(mat[, f], censor_q, na.rm = TRUE)
    mat[mat[, f] < lod, f] <- NA
}

colnames(mat) <- paste0('M', seq_len(n_features))
qc_rows <- which(is_qc)
bio_rows <- which(!is_qc)

cat(sprintf('Start: %d samples (%d QC), %d features\n', n_total, n_qc, n_features))

# --- 1. Detection-rate filter BEFORE imputation (never impute a mostly-missing feature) ---
detection_min <- 0.5
detect_rate <- colMeans(!is.na(mat))
mat <- mat[, detect_rate >= detection_min]
cat(sprintf('After detection-rate>=%.0f%% filter: %d features\n', detection_min * 100, ncol(mat)))

# --- 2. Within-batch drift correction: per-feature LOESS of QC intensity vs injection order ---
# Validated on HELD-OUT QCs (split QCs into fit/test) -- the only non-circular QC metric.
qc_fit <- qc_rows[seq(1, length(qc_rows), by = 2)]
qc_test <- setdiff(qc_rows, qc_fit)

drift_correct <- function(data, inj, qc_fit_rows) {
    out <- data
    for (f in colnames(data)) {
        y <- data[qc_fit_rows, f]
        x <- inj[qc_fit_rows]
        ok <- !is.na(y)
        if (sum(ok) < 4)
            next
        fit <- loess(y[ok] ~ x[ok], span = 0.75, degree = 2)
        # Samples beyond the QC range are EXTRAPOLATED (unstable); fall back to the QC median there.
        trend <- suppressWarnings(predict(fit, inj))
        trend[is.na(trend) | trend <= 0] <- median(y[ok], na.rm = TRUE)
        out[, f] <- data[, f] / trend * median(y[ok], na.rm = TRUE)
    }
    out
}

mat_corr <- drift_correct(mat, inj, qc_fit)

rsd <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
heldout_rsd_before <- median(apply(mat[qc_test, ], 2, rsd), na.rm = TRUE)
heldout_rsd_after <- median(apply(mat_corr[qc_test, ], 2, rsd), na.rm = TRUE)
cat(sprintf('Held-out QC median RSD: %.3f -> %.3f (correction is real only if this drops)\n',
            heldout_rsd_before, heldout_rsd_after))

# --- 3. Robust D-ratio + RSD filter (technical SD / biological SD), MAD-based for skewed data ---
dratio_max <- 0.5
rsd_max <- 0.3
sd_qc <- apply(mat_corr[qc_rows, ], 2, mad, na.rm = TRUE)
sd_bio <- apply(mat_corr[bio_rows, ], 2, mad, na.rm = TRUE)
dratio <- sd_qc / sd_bio
qc_rsd <- apply(mat_corr[qc_rows, ], 2, rsd)
keep <- dratio <= dratio_max & qc_rsd <= rsd_max
keep[is.na(keep)] <- FALSE
mat_corr <- mat_corr[, keep]
cat(sprintf('After D-ratio<=%.1f & RSD<=%.0f%% filter: %d features\n',
            dratio_max, rsd_max * 100, ncol(mat_corr)))

# --- 4. Mechanism-aware imputation: left-censored holes drawn from a low truncated distribution ---
# MNAR (below-LOD) -> draw plausible low values, NOT a single half-min constant (which zeroes
# variance and inflates false significance). A MAR method (kNN/RF) here would pull censored
# group means up and erase on/off biology.
impute_left_censored <- function(data) {
    out <- data
    for (f in colnames(data)) {
        miss <- is.na(data[, f])
        if (!any(miss))
            next
        obs <- data[!miss, f]
        floor_val <- quantile(obs, 0.01, na.rm = TRUE)
        draws <- rnorm(sum(miss), mean = floor_val, sd = sd(obs, na.rm = TRUE) * 0.1)
        out[miss, f] <- pmax(draws, min(obs, na.rm = TRUE) * 0.5)
    }
    out
}

mat_imp <- impute_left_censored(mat_corr)
cat(sprintf('Missing values remaining after imputation: %d\n', sum(is.na(mat_imp))))

# --- 5. PQN sample normalization: median feature quotient vs a QC reference spectrum ---
pqn_normalize <- function(data, reference_rows) {
    reference <- apply(data[reference_rows, , drop = FALSE], 2, median, na.rm = TRUE)
    quotients <- sweep(data, 2, reference, '/')
    factors <- apply(quotients, 1, median, na.rm = TRUE)
    list(data = data / factors, factors = factors)
}

pqn <- pqn_normalize(mat_imp, qc_rows)
mat_norm <- pqn$data

# Guardrail: the PQN factor should track true dilution, NOT group. A high factor-vs-group
# correlation means the normalization is eating the biological effect.
factor_dilution_cor <- suppressWarnings(cor(pqn$factors[bio_rows], dilution[bio_rows]))
factor_group_cor <- suppressWarnings(cor(pqn$factors[bio_rows],
                                          as.integer(sample_group[bio_rows] == 'case')))
cat(sprintf('PQN factor vs true dilution: %.3f (want high -- PQN recovered the dilution)\n',
            factor_dilution_cor))
verdict <- if (abs(factor_group_cor) > 0.3) 'TRIPPED: a strong group effect leaked into the PQN factor; normalize to a measured external quantity instead' else 'ok: factor is independent of group'
cat(sprintf('PQN factor vs group: %.3f (%s)\n', factor_group_cor, verdict))

out_dir <- tempdir()
out_file <- file.path(out_dir, 'normalized_feature_table.csv')
write.csv(mat_norm, out_file)
cat(sprintf('Wrote %s (%d samples x %d features)\n', out_file, nrow(mat_norm), ncol(mat_norm)))

file.remove(out_file)
cat('Cleaned up temporary output\n')
