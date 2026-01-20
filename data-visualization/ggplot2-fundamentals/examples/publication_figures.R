library(ggplot2)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

theme_publication <- function(base_size = 12) {
    theme_bw(base_size = base_size) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 0.5),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5)
    )
}

publication_colors <- c(
    'Control' = '#4DBBD5',
    'Treatment' = '#E64B35',
    'Up' = '#E64B35',
    'Down' = '#4DBBD5',
    'NS' = 'grey60'
)

create_volcano <- function(res, fdr_threshold = 0.05, lfc_threshold = 1, top_n = 10) {
    res <- res %>%
        mutate(
            significance = case_when(
                padj < fdr_threshold & log2FoldChange > lfc_threshold ~ 'Up',
                padj < fdr_threshold & log2FoldChange < -lfc_threshold ~ 'Down',
                TRUE ~ 'NS'
            ),
            label = ifelse(rank(padj) <= top_n & padj < fdr_threshold, gene, '')
        )

    ggplot(res, aes(log2FoldChange, -log10(pvalue), color = significance)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(aes(label = label), color = 'black', size = 3, max.overlaps = 20) +
        scale_color_manual(values = publication_colors) +
        geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                   linetype = 'dashed', color = 'grey40') +
        geom_hline(yintercept = -log10(fdr_threshold),
                   linetype = 'dashed', color = 'grey40') +
        labs(x = expression(Log[2]~Fold~Change),
             y = expression(-Log[10]~P-value),
             color = 'Significance') +
        theme_publication()
}

create_boxplot <- function(df, x_var, y_var, fill_var = NULL) {
    p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]]))

    if (!is.null(fill_var)) {
        p <- p + aes(fill = .data[[fill_var]])
    }

    p +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
        scale_fill_brewer(palette = 'Set2') +
        labs(x = NULL) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

create_pca_plot <- function(pca_df, color_var, shape_var = NULL) {
    p <- ggplot(pca_df, aes(PC1, PC2, color = .data[[color_var]]))

    if (!is.null(shape_var)) {
        p <- p + aes(shape = .data[[shape_var]])
    }

    p +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_brewer(palette = 'Set1') +
        labs(x = paste0('PC1 (', round(pca_df$var_explained[1], 1), '%)'),
             y = paste0('PC2 (', round(pca_df$var_explained[2], 1), '%)')) +
        theme_publication()
}

save_publication_figure <- function(plot, filename, width = 7, height = 5) {
    ggsave(paste0(filename, '.pdf'), plot, width = width, height = height, units = 'in')
    ggsave(paste0(filename, '.png'), plot, width = width, height = height, units = 'in', dpi = 300)
}

create_multi_panel <- function(p1, p2, p3, p4 = NULL) {
    if (is.null(p4)) {
        combined <- (p1 | p2) / p3 +
            plot_annotation(tag_levels = 'A') +
            plot_layout(heights = c(1, 1))
    } else {
        combined <- (p1 | p2) / (p3 | p4) +
            plot_annotation(tag_levels = 'A')
    }
    combined & theme_publication()
}
