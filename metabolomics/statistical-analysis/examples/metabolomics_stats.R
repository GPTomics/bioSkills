# Reference: R stats (base), mixOmics 6.28+, ggplot2 3.5+ | Verify API if version differs
library(tidyverse)
library(mixOmics)

data <- as.matrix(read.csv('normalized_data.csv', row.names = 1))
groups <- factor(read.csv('sample_info.csv')$group)

cat('Samples:', nrow(data), '\n')
cat('Features:', ncol(data), '\n')
cat('Groups:', levels(groups), '\n')

# Log2 transform before statistical testing
data_log2 <- log2(data)

# 1. Univariate analysis (Welch's t-test on log2-transformed data)
ttest <- function(x, groups) {
    treat_vals <- x[groups == levels(groups)[2]]
    ctrl_vals <- x[groups == levels(groups)[1]]
    test <- t.test(treat_vals, ctrl_vals)
    c(pvalue = test$p.value, log2fc = mean(treat_vals) - mean(ctrl_vals))
}

univariate <- as.data.frame(t(apply(data_log2, 2, ttest, groups = groups)))
univariate$padj <- p.adjust(univariate$pvalue, method = 'BH')
univariate$feature <- rownames(univariate)

sig_univariate <- univariate[univariate$padj < 0.05, ]
cat('\nSignificant features (Welch t-test, padj<0.05):', nrow(sig_univariate), '\n')

# 2. PCA
pca <- prcomp(data_log2, scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:2] * 100

cat('\nPCA variance explained:')
cat('\n  PC1:', round(var_explained[1], 1), '%')
cat('\n  PC2:', round(var_explained[2], 1), '%\n')

# 3. PLS-DA
plsda <- plsda(data_log2, groups, ncomp = 2)
vip_scores <- vip(plsda)
top_vip <- sort(vip_scores[, 2], decreasing = TRUE)[1:10]

cat('\nTop 10 VIP features:\n')
print(round(top_vip, 2))

# 4. Combine results
results <- univariate
results$vip <- vip_scores[results$feature, 2]
results$significant <- results$padj < 0.05 & results$vip > 1

cat('\nFeatures with padj<0.05 AND VIP>1:', sum(results$significant), '\n')

write.csv(results, 'statistical_results.csv', row.names = FALSE)
cat('\nResults saved to statistical_results.csv\n')
