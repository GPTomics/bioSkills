# Reference: limma 3.58+, ashr 2.2+ | Verify API if version differs
# Differential protein abundance with limma and ashr FC shrinkage
library(limma)

protein_matrix <- as.matrix(read.csv('protein_intensities.csv', row.names = 1))
sample_info <- read.csv('sample_info.csv')

# Preprocess: log2 transform raw intensities
log2_matrix <- log2(protein_matrix)
log2_matrix[!is.finite(log2_matrix)] <- NA

# Normalize: median centering via limma scale method
log2_norm <- normalizeBetweenArrays(log2_matrix, method = 'scale')

# Design and contrasts
sample_info$condition <- factor(sample_info$condition, levels = c('Control', 'Treatment'))
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(log2_norm, design)
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
    results$logFC <- shrunk$result$PosteriorMean[match(rownames(results), names(fit2$coefficients[, 1]))]
    results$lfsr <- shrunk$result$lfsr[match(rownames(results), names(fit2$coefficients[, 1]))]
}

results$protein <- rownames(results)
results$significant <- results$adj.P.Val < 0.05

cat('Total proteins tested:', nrow(results), '\n')
cat('Significant (adj.p < 0.05):', sum(results$significant), '\n')
cat('Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('Down-regulated:', sum(results$significant & results$logFC < 0), '\n')

write.csv(results, 'differential_results.csv', row.names = FALSE)
