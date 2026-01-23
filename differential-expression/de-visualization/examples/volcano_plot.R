# Create volcano plot from DESeq2 results

library(DESeq2)
library(ggplot2)
library(ggrepel)

# Simulate DESeq2 results (in practice, use your own dds object)
set.seed(42)
n_genes <- 1000
res_df <- data.frame(
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    pvalue = 10^(-runif(n_genes, 0, 5)),
    padj = p.adjust(10^(-runif(n_genes, 0, 5)), method = 'BH'),
    baseMean = 10^runif(n_genes, 1, 4),
    row.names = paste0('gene', 1:n_genes)
)
res_df$gene <- rownames(res_df)

# Define significance
# padj < 0.05: Standard FDR cutoff (5% false discovery rate)
# |log2FC| > 1: Requires 2-fold change - conservative, reduces noise
# For exploratory analysis, use |log2FC| > 0.5 (1.4-fold) or padj < 0.1
res_df$significant <- ifelse(
    res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
    ifelse(res_df$log2FoldChange > 1, 'Up', 'Down'),
    'NS'
)

# Top genes to label
top_genes <- head(res_df[order(res_df$padj), ], 10)

# Create volcano plot
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c('Down' = 'blue', 'NS' = 'grey60', 'Up' = 'red')) +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'grey40') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey40') +
    geom_text_repel(data = top_genes, aes(label = gene), color = 'black',
                    size = 3, max.overlaps = 15) +
    labs(x = 'log2 Fold Change', y = '-log10(p-value)',
         title = 'Volcano Plot', color = 'Significance') +
    theme_bw() +
    theme(legend.position = 'bottom')

print(p)
ggsave('volcano_plot.png', p, width = 8, height = 6, dpi = 300)
