# Reference: limma 3.58+, ashr 2.2+, mixOmics 6.28+ | Verify API if version differs
library(limma)
library(mixOmics)

data <- as.matrix(read.csv('feature_intensities.csv', row.names = 1))
sample_info <- read.csv('sample_info.csv')

cat('Samples:', ncol(data), '\n')
cat('Features:', nrow(data), '\n')
cat('Groups:', levels(factor(sample_info$group)), '\n')

# Preprocess: log2 transform
log2_data <- log2(data)
log2_data[!is.finite(log2_data)] <- NA

# PQN normalization
reference <- apply(log2_data, 1, median, na.rm = TRUE)
quotients <- sweep(log2_data, 1, reference, '/')
norm_factors <- apply(quotients, 2, median, na.rm = TRUE)
normalized <- sweep(log2_data, 2, norm_factors, '/')

# Filter features present in >50% of samples
valid <- rowSums(!is.na(normalized)) > ncol(normalized) * 0.5
normalized <- normalized[valid, ]
cat('Features after filtering:', nrow(normalized), '\n')

# 1. limma differential analysis
sample_info$group <- factor(sample_info$group)
design <- model.matrix(~0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

imputed <- normalized
imputed[is.na(imputed)] <- min(imputed, na.rm = TRUE) - 1

fit <- lmFit(imputed, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')

# FC shrinkage with ashr for more accurate effect size estimates
if (requireNamespace('ashr', quietly = TRUE)) {
    library(ashr)
    se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]
    shrunk <- ash(fit2$coefficients[, 1], se, mixcompdist = 'normal')
    results$logFC_raw <- results$logFC
    results$logFC_shrunk <- shrunk$result$PosteriorMean[match(rownames(results), names(fit2$coefficients[, 1]))]
    results$lfsr <- shrunk$result$lfsr[match(rownames(results), names(fit2$coefficients[, 1]))]
}

results$significant <- results$adj.P.Val < 0.05
cat('\nSignificant (adj.P.Val < 0.05):', sum(results$significant), '\n')

# 2. PCA
pca <- prcomp(t(normalized), scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:2] * 100
cat('\nPCA variance explained:')
cat('\n  PC1:', round(var_explained[1], 1), '%')
cat('\n  PC2:', round(var_explained[2], 1), '%\n')

# 3. PLS-DA
plsda <- plsda(t(normalized), sample_info$group, ncomp = 2)
vip_scores <- vip(plsda)
top_vip <- sort(vip_scores[, 2], decreasing = TRUE)[1:10]
cat('\nTop 10 VIP features:\n')
print(round(top_vip, 2))

# 4. Combine limma + VIP results
results$feature <- rownames(results)
results$vip <- vip_scores[results$feature, 2]
results$significant_both <- results$adj.P.Val < 0.05 & results$vip > 1

cat('\nFeatures with adj.P.Val<0.05 AND VIP>1:', sum(results$significant_both, na.rm = TRUE), '\n')

write.csv(results, 'statistical_results.csv', row.names = FALSE)
cat('\nResults saved to statistical_results.csv\n')
