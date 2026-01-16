# Filter and summarize DESeq2 results

library(DESeq2)
library(dplyr)

# Simulate DESeq2 results (in practice, use your dds object)
set.seed(42)
n_genes <- 5000
res_df <- data.frame(
    baseMean = 10^runif(n_genes, 1, 4),
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 2),
    pvalue = 10^(-runif(n_genes, 0, 5)),
    row.names = paste0('gene', 1:n_genes)
)
res_df$padj <- p.adjust(res_df$pvalue, method = 'BH')
res_df$gene <- rownames(res_df)

# Summary statistics
cat('=== Summary ===\n')
cat(sprintf('Total genes tested: %d\n', sum(!is.na(res_df$padj))))
cat(sprintf('Significant (padj < 0.05): %d\n', sum(res_df$padj < 0.05, na.rm = TRUE)))
cat(sprintf('Up-regulated: %d\n', sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)))
cat(sprintf('Down-regulated: %d\n', sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)))

# Filter significant genes
sig_genes <- res_df %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)

cat(sprintf('\nSignificant with |log2FC| > 1: %d genes\n', nrow(sig_genes)))

# Top 10 genes
cat('\n=== Top 10 Genes ===\n')
print(head(sig_genes[, c('gene', 'baseMean', 'log2FoldChange', 'padj')], 10))

# Export
write.csv(sig_genes, 'significant_genes.csv', row.names = FALSE)
cat('\nResults saved to significant_genes.csv\n')
