library(ggplot2)
library(patchwork)
library(dplyr)

theme_publication <- function(base_size = 10) {
    theme_bw(base_size = base_size) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 0.5),
        axis.text = element_text(color = 'black'),
        strip.background = element_blank(),
        plot.title = element_text(size = 11, face = 'bold')
    )
}

set.seed(42)
n <- 100
df <- data.frame(
    gene = paste0('Gene', 1:n),
    log2FC = rnorm(n, 0, 2),
    pvalue = 10^(-runif(n, 0, 10)),
    expression = rnorm(n, 10, 3),
    group = sample(c('Control', 'Treatment'), n, replace = TRUE)
)
df$padj <- p.adjust(df$pvalue, method = 'BH')
df$significant <- df$padj < 0.05 & abs(df$log2FC) > 1

p_volcano <- ggplot(df, aes(log2FC, -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c('grey60', 'red3')) +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~P-value),
         title = 'Differential Expression') +
    theme_publication() +
    theme(legend.position = 'none')

pca_data <- data.frame(
    PC1 = rnorm(20),
    PC2 = rnorm(20),
    group = rep(c('Control', 'Treatment'), each = 10)
)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95) +
    scale_color_manual(values = c('Control' = '#4DBBD5', 'Treatment' = '#E64B35')) +
    labs(title = 'PCA', color = 'Group') +
    theme_publication()

expr_df <- data.frame(
    gene = rep(paste0('Gene', 1:5), each = 10),
    expression = c(rnorm(10, 8), rnorm(10, 12), rnorm(10, 6),
                   rnorm(10, 10), rnorm(10, 9)),
    group = rep(c('Control', 'Treatment'), 25)
)

p_boxplot <- ggplot(expr_df, aes(gene, expression, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, size = 1) +
    scale_fill_manual(values = c('Control' = '#4DBBD5', 'Treatment' = '#E64B35')) +
    labs(x = NULL, y = 'Expression', title = 'Top DE Genes', fill = 'Group') +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

bar_df <- data.frame(
    category = c('Upregulated', 'Downregulated', 'Not Significant'),
    count = c(sum(df$significant & df$log2FC > 0),
              sum(df$significant & df$log2FC < 0),
              sum(!df$significant))
)

p_bar <- ggplot(bar_df, aes(category, count, fill = category)) +
    geom_col() +
    scale_fill_manual(values = c('Upregulated' = '#E64B35', 'Downregulated' = '#4DBBD5',
                                  'Not Significant' = 'grey60')) +
    labs(x = NULL, y = 'Number of Genes', title = 'DE Summary') +
    theme_publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))

combined <- (p_volcano | p_pca) / (p_boxplot | p_bar) +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = 'collect') &
    theme(
        plot.tag = element_text(face = 'bold', size = 12),
        legend.position = 'bottom'
    )

ggsave('Figure1.pdf', combined, width = 10, height = 8, units = 'in')
ggsave('Figure1.png', combined, width = 10, height = 8, units = 'in', dpi = 300)

message('Figure saved: Figure1.pdf and Figure1.png')
