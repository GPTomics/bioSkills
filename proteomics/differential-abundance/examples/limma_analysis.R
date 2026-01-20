# Differential protein abundance with limma
library(limma)

protein_matrix <- read.csv('normalized_imputed.csv', row.names = 1)
sample_info <- read.csv('sample_info.csv')

sample_info$condition <- factor(sample_info$condition, levels = c('Control', 'Treatment'))
design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(as.matrix(protein_matrix), design)

contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf, adjust.method = 'BH')
results$protein <- rownames(results)

results$significant <- abs(results$logFC) > 1 & results$adj.P.Val < 0.05

cat('Total proteins tested:', nrow(results), '\n')
cat('Significant (|log2FC| > 1, adj.p < 0.05):', sum(results$significant), '\n')
cat('Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('Down-regulated:', sum(results$significant & results$logFC < 0), '\n')

write.csv(results, 'differential_results.csv', row.names = FALSE)
