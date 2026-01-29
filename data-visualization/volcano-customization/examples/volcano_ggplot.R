library(ggplot2)
library(ggrepel)
library(dplyr)

# --- ALTERNATIVE: Use real DE results ---
# For realistic volcanos, use DESeq2/edgeR output:
#
# library(DESeq2)
# library(airway)
# data('airway')
# dds <- DESeqDataSet(airway, design = ~ dex)
# dds <- DESeq(dds)
# res <- results(dds)
# df <- as.data.frame(res) %>%
#     tibble::rownames_to_column('gene') %>%
#     filter(!is.na(padj))

set.seed(42)

# Simulate DE results with realistic distribution
# Most genes show no change (centered at 0), few are significant
n_genes <- 5000
df <- data.frame(
    gene = paste0('Gene', 1:n_genes),
    log2FoldChange = c(
        rnorm(200, mean = 2, sd = 0.5),    # Upregulated
        rnorm(150, mean = -1.8, sd = 0.6), # Downregulated
        rnorm(4650, mean = 0, sd = 0.5)    # No change
    ),
    stringsAsFactors = FALSE
)

# P-values correlate with fold change magnitude (realistic)
# Larger effects tend to have smaller p-values
df$pvalue <- 10^(-abs(df$log2FoldChange) * runif(n_genes, 0.5, 3) - runif(n_genes, 0, 2))
df$pvalue <- pmin(df$pvalue, 1)

# Adjust for multiple testing
# padj will be larger than pvalue, some lose significance
df$padj <- p.adjust(df$pvalue, method = 'BH')

# Thresholds for significance
# FC > 1: Standard 2-fold change cutoff
# padj < 0.05: Standard FDR threshold
fc_threshold <- 1
padj_threshold <- 0.05

df <- df %>%
    mutate(significance = case_when(
        padj < padj_threshold & log2FoldChange > fc_threshold ~ 'Up',
        padj < padj_threshold & log2FoldChange < -fc_threshold ~ 'Down',
        TRUE ~ 'NS'
    ))

# Count significant genes
cat('Significant genes:\n')
cat(sprintf('  Upregulated: %d\n', sum(df$significance == 'Up')))
cat(sprintf('  Downregulated: %d\n', sum(df$significance == 'Down')))
cat(sprintf('  Not significant: %d\n', sum(df$significance == 'NS')))

# Basic volcano with threshold lines
p1 <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(Up = '#E64B35', Down = '#4DBBD5', NS = 'gray70'),
                       name = 'Regulation') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'gray40') +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold),
               linetype = 'dashed', color = 'gray40') +
    labs(x = expression(log[2]~Fold~Change),
         y = expression(-log[10]~p-value),
         title = 'Differential Expression') +
    theme_classic(base_size = 12) +
    theme(legend.position = 'right')

ggsave('volcano_basic.pdf', p1, width = 8, height = 6)

# Volcano with top gene labels
# Label top 15 by p-value among significant genes
top_genes <- df %>%
    filter(significance != 'NS') %>%
    arrange(pvalue) %>%
    head(15)

p2 <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significance), alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c(Up = '#E64B35', Down = '#4DBBD5', NS = 'gray70')) +
    geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = 'gray50',
        segment.size = 0.3
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'gray40') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'gray40') +
    labs(x = expression(log[2]~Fold~Change),
         y = expression(-log[10]~p-value),
         title = 'DE Genes with Labels') +
    theme_classic() +
    theme(legend.position = 'none')

ggsave('volcano_labeled.pdf', p2, width = 8, height = 6)

# Volcano highlighting specific genes
genes_of_interest <- c('Gene1', 'Gene201', 'Gene350', 'Gene25', 'Gene180')
highlight_df <- df %>% filter(gene %in% genes_of_interest)

p3 <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significance), alpha = 0.3, size = 1.5) +
    scale_color_manual(values = c(Up = '#E64B35', Down = '#4DBBD5', NS = 'gray70')) +
    geom_point(data = highlight_df, color = 'black', size = 3, shape = 21,
               fill = 'yellow', stroke = 1) +
    geom_text_repel(data = highlight_df, aes(label = gene),
                    fontface = 'bold', size = 4,
                    nudge_y = 2, box.padding = 0.5) +
    labs(x = expression(log[2]~Fold~Change),
         y = expression(-log[10]~p-value),
         title = 'Genes of Interest Highlighted') +
    theme_classic() +
    theme(legend.position = 'none')

ggsave('volcano_highlighted.pdf', p3, width = 8, height = 6)

message('Volcano plots saved: volcano_basic.pdf, volcano_labeled.pdf, volcano_highlighted.pdf')
