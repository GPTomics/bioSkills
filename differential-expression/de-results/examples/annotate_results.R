# Add gene annotations to DE results

library(dplyr)

# Example: annotate with custom annotation file
# In practice, use org.Hs.eg.db or biomaRt

# Simulate results
set.seed(42)
n_genes <- 100
res_df <- data.frame(
    gene_id = paste0('ENSG', sprintf('%011d', 1:n_genes)),
    baseMean = 10^runif(n_genes, 1, 4),
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    padj = 10^(-runif(n_genes, 0, 3))
)

# Simulate annotation file
gene_info <- data.frame(
    gene_id = paste0('ENSG', sprintf('%011d', 1:n_genes)),
    symbol = paste0('GENE', 1:n_genes),
    description = paste('Description for gene', 1:n_genes)
)

# Merge annotations
res_annotated <- res_df %>%
    left_join(gene_info, by = 'gene_id') %>%
    arrange(padj)

# View annotated results
cat('=== Annotated Results (top 10) ===\n')
print(head(res_annotated[, c('gene_id', 'symbol', 'log2FoldChange', 'padj')], 10))

# Filter significant and export
sig_annotated <- res_annotated %>%
    filter(padj < 0.05) %>%
    select(gene_id, symbol, description, log2FoldChange, padj)

write.csv(sig_annotated, 'annotated_significant.csv', row.names = FALSE)
cat('\nAnnotated results saved to annotated_significant.csv\n')
