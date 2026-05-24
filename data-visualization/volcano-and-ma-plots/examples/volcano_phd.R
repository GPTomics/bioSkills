# Reference: DESeq2 1.42+, EnhancedVolcano 1.20+, ggplot2 3.5+, ggrepel 0.9.5+ | Verify API if version differs

# PhD-level volcano + MA plot from DESeq2 results.
# Encodes the four correctness traps: (1) shrunken LFC, (2) padj on y-axis,
# (3) label by combined rank, (4) max.overlaps = Inf for ggrepel.

library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

# 1. SHRINKAGE -- apeglm by default (Zhu et al 2019)
# coef = name from resultsNames(dds); ashr is the alternative when contrast= is needed
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
# alt: res <- lfcShrink(dds, contrast = c('condition', 'treated', 'control'), type = 'ashr')

# 2. CATEGORICAL SIGNIFICANCE -- threshold on padj NOT pvalue
fdr <- 0.05
lfc_threshold <- 1
res_df <- as.data.frame(res) %>%
    rownames_to_column('gene') %>%
    mutate(
        significance = case_when(
            is.na(padj) ~ 'NS',
            padj < fdr & log2FoldChange >  lfc_threshold ~ 'Up',
            padj < fdr & log2FoldChange < -lfc_threshold ~ 'Down',
            TRUE ~ 'NS'),
        neg_log10_p = -log10(pvalue))

# 3. LABEL SELECTION -- combined rank, not pure top-N-by-p
genes_of_interest <- c('TP53', 'MYC', 'BRCA1')
top_by_rank <- res_df %>%
    filter(significance != 'NS') %>%
    mutate(rank_score = -log10(pvalue) * abs(log2FoldChange)) %>%
    arrange(desc(rank_score)) %>% head(10) %>% pull(gene)
labels <- union(genes_of_interest, top_by_rank)
res_df$label <- ifelse(res_df$gene %in% labels, res_df$gene, '')

# 4. PLOT -- Okabe-Ito categorical palette (Wong 2011)
okabe_ito <- c(Up = '#D55E00', Down = '#0072B2', NS = '#999999')

# Optional: cap y-axis when extreme p dominate
y_cap <- 50  # if max(-log10 p) > 50 the cap improves readability without losing data

p_volcano <- ggplot(res_df, aes(log2FoldChange, neg_log10_p, color = significance)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = okabe_ito, name = NULL) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
    # threshold line at the padj cutoff projected to the raw-p axis is approximate;
    # for an exact line on padj you must change the y-axis to padj directly
    geom_hline(yintercept = -log10(fdr), linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
    geom_text_repel(aes(label = label), color = 'black', size = 3,
                    max.overlaps = Inf,                    # critical -- default 10 drops labels
                    box.padding = 0.4, segment.size = 0.2, min.segment.length = 0) +
    coord_cartesian(ylim = c(0, y_cap)) +                  # keep all points; clip display only
    labs(x = expression(log[2]~'fold change (shrunken)'),
         y = expression(-log[10]~italic(p))) +
    theme_classic(base_size = 10) +
    theme(panel.grid = element_blank())

# 5. MA PLOT as the shrinkage diagnostic
# fan shape with extreme |LFC| at low baseMean indicates inadequate shrinkage
p_ma <- ggplot(res_df, aes(log10(baseMean), log2FoldChange,
                            color = significance != 'NS')) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c(`FALSE` = '#999999', `TRUE` = '#D55E00'),
                       guide = 'none') +
    geom_hline(yintercept = 0, color = 'black', linewidth = 0.4) +
    labs(x = expression(log[10]~'mean normalized count'),
         y = expression(log[2]~'fold change (shrunken)')) +
    theme_classic(base_size = 10) +
    theme(panel.grid = element_blank())

# 6. SAVE -- cairo_pdf embeds TrueType; >5000 features -> rasterize the point layer
# ggsave with raster requires ggplot2 3.5+ and ragg/Cairo backend
ggsave('volcano.pdf', p_volcano, width = 89, height = 90, units = 'mm',
       device = cairo_pdf)
ggsave('ma_plot.pdf', p_ma, width = 89, height = 70, units = 'mm',
       device = cairo_pdf)

# 7. ENHANCEDVOLCANO equivalent -- aware of selectLab + pCutoff trap
# genes in selectLab that fail pCutoff/FCcutoff are silently unlabeled
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange', y = 'padj',                       # y on padj NOT pvalue
    pCutoff = fdr, FCcutoff = lfc_threshold,
    selectLab = labels,
    drawConnectors = TRUE, widthConnectors = 0.3,
    maxoverlapsConnectors = Inf,                            # match ggrepel guidance
    col = c('grey60', '#0072B2', '#56B4E9', '#D55E00'),
    pointSize = 1.5, labSize = 3, colAlpha = 0.6,
    legendPosition = 'right')
