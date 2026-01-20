library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)

set.seed(42)
n <- 500

de_results <- data.frame(
    gene = paste0('Gene', 1:n),
    log2FoldChange = rnorm(n, 0, 1.5),
    pvalue = 10^(-runif(n, 0, 8)),
    baseMean = 10^runif(n, 1, 4)
)
de_results$padj <- p.adjust(de_results$pvalue, method = 'BH')

volcano_plot <- function(res, fdr = 0.05, lfc = 1, top_n = 10) {
    res <- res %>%
        mutate(
            significance = case_when(
                padj < fdr & log2FoldChange > lfc ~ 'Up',
                padj < fdr & log2FoldChange < -lfc ~ 'Down',
                TRUE ~ 'NS'
            ),
            label = ifelse(rank(padj) <= top_n & significance != 'NS', gene, '')
        )

    ggplot(res, aes(log2FoldChange, -log10(pvalue), color = significance)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(aes(label = label), color = 'black', size = 3, max.overlaps = 20) +
        scale_color_manual(values = c('Up' = '#E64B35', 'Down' = '#4DBBD5', 'NS' = 'grey60')) +
        geom_vline(xintercept = c(-lfc, lfc), linetype = 'dashed', color = 'grey40') +
        geom_hline(yintercept = -log10(fdr), linetype = 'dashed', color = 'grey40') +
        labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~P-value),
             title = 'Volcano Plot') +
        theme_bw() + theme(panel.grid = element_blank())
}

ma_plot <- function(res, fdr = 0.05) {
    res <- res %>%
        mutate(significant = padj < fdr & !is.na(padj))

    ggplot(res, aes(log10(baseMean), log2FoldChange, color = significant)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c('FALSE' = 'grey60', 'TRUE' = '#E64B35')) +
        geom_hline(yintercept = 0, color = 'black', linewidth = 0.5) +
        labs(x = expression(Log[10]~Mean~Expression), y = expression(Log[2]~Fold~Change),
             title = 'MA Plot') +
        theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
}

pca_data <- data.frame(
    PC1 = c(rnorm(10, -2), rnorm(10, 2)),
    PC2 = c(rnorm(10, 0), rnorm(10, 0)),
    condition = rep(c('Control', 'Treatment'), each = 10),
    sample = paste0('Sample', 1:20)
)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 'dashed') +
    scale_color_manual(values = c('Control' = '#4DBBD5', 'Treatment' = '#E64B35')) +
    labs(x = 'PC1 (45.2%)', y = 'PC2 (18.7%)', title = 'PCA Plot') +
    theme_bw() + theme(panel.grid = element_blank())

p_volcano <- volcano_plot(de_results)
p_ma <- ma_plot(de_results)

combined <- (p_volcano | p_ma) / pca_plot +
    plot_annotation(tag_levels = 'A')

ggsave('omics_plots.pdf', combined, width = 12, height = 10)
message('Saved: omics_plots.pdf')
