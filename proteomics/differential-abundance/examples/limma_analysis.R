# Reference: limma 3.58+, ashr 2.2+ | Verify API if version differs
# Differential protein abundance with limma empirical-Bayes moderation, treat() min-FC, and ashr FC shrinkage.
# At n=3-5 the per-protein variance is unusable raw; trend=TRUE + robust=TRUE moderation is the load-bearing step.
library(limma)

set.seed(0)
n_proteins <- 400
n_true <- 40
n_per_group <- 4  # small n typical of proteomics: only ~6 residual df before moderation
LFC_THRESHOLD <- log2(1.2)  # 1.2-fold floor; treat tests against this null without double-filter FDR inflation

# Simulate a log2 intensity matrix: first n_true proteins up in Treatment (0.4 within-group spread, clean demo)
log2_matrix <- matrix(rnorm(n_proteins * 2 * n_per_group, 20, 0.4), nrow = n_proteins)
log2_matrix[1:n_true, (n_per_group + 1):(2 * n_per_group)] <-
    log2_matrix[1:n_true, (n_per_group + 1):(2 * n_per_group)] + 1.5  # true up-regulation in case
rownames(log2_matrix) <- sprintf('P%04d', seq_len(n_proteins))

sample_info <- data.frame(condition = factor(rep(c('Control', 'Treatment'), each = n_per_group),
                                             levels = c('Control', 'Treatment')))

# Normalize log2 intensities (median centering); summarization/normalization mechanics live in quantification
log2_norm <- normalizeBetweenArrays(log2_matrix, method = 'scale')

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(log2_norm, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)  # trend mandatory for intensity data; robust Winsorizes outliers

results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')  # columns include adj.P.Val, no $FDR
results$significant <- results$adj.P.Val < 0.05

# Minimum-fold-change testing: treat()+topTreat(), never topTable(lfc=...) nor a post-hoc volcano double filter
fit_treat <- treat(fit2, lfc = LFC_THRESHOLD)
treat_results <- topTreat(fit_treat, coef = 1, number = Inf)  # topTreat omits the B column

# Fold-change shrinkage for effect-size recovery (report alongside raw logFC; use raw FC for GSEA)
if (requireNamespace('ashr', quietly = TRUE)) {
    se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]
    shrunk <- ashr::ash(fit2$coefficients[, 1], se, mixcompdist = 'normal')
    results$logFC_shrunk <- shrunk$result$PosteriorMean[match(rownames(results), rownames(fit2$coefficients))]
}

cat('Tested:', nrow(results), '\n')
cat('Significant (adj.P.Val < 0.05):', sum(results$significant), '\n')
cat('Pass 1.2-fold treat() (adj.P.Val < 0.05):', sum(treat_results$adj.P.Val < 0.05), '\n')
print(head(results[order(results$adj.P.Val), c('logFC', 'adj.P.Val')]))
