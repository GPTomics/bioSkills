library(tidyverse)
library(mixOmics)

# Load data
data <- as.matrix(read.csv('normalized_data.csv', row.names = 1))
groups <- factor(read.csv('sample_info.csv')$group)

cat('Samples:', nrow(data), '\n')
cat('Features:', ncol(data), '\n')
cat('Groups:', levels(groups), '\n')

# 1. Univariate analysis (t-test)
ttest <- function(x, groups) {
    test <- t.test(x ~ groups)
    c(pvalue = test$p.value, diff = diff(test$estimate))
}

univariate <- as.data.frame(t(apply(data, 2, ttest, groups = groups)))
univariate$fdr <- p.adjust(univariate$pvalue, method = 'BH')
univariate$feature <- rownames(univariate)

sig_univariate <- univariate[univariate$fdr < 0.05, ]
cat('\nSignificant features (t-test, FDR<0.05):', nrow(sig_univariate), '\n')

# 2. PCA
pca <- prcomp(data, scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:2] * 100

cat('\nPCA variance explained:')
cat('\n  PC1:', round(var_explained[1], 1), '%')
cat('\n  PC2:', round(var_explained[2], 1), '%\n')

# 3. PLS-DA
plsda <- plsda(data, groups, ncomp = 2)
vip_scores <- vip(plsda)
top_vip <- sort(vip_scores[, 2], decreasing = TRUE)[1:10]

cat('\nTop 10 VIP features:\n')
print(round(top_vip, 2))

# 4. Combine results
results <- univariate
results$vip <- vip_scores[results$feature, 2]
results$significant <- results$fdr < 0.05 & results$vip > 1

cat('\nFeatures with FDR<0.05 AND VIP>1:', sum(results$significant), '\n')

# Save
write.csv(results, 'statistical_results.csv', row.names = FALSE)
cat('\nResults saved to statistical_results.csv\n')
