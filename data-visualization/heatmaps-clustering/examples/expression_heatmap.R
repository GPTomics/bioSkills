library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

set.seed(42)
n_genes <- 100
n_samples <- 20

mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
rownames(mat) <- paste0('Gene', 1:n_genes)
colnames(mat) <- paste0('Sample', 1:n_samples)

metadata <- data.frame(
    sample = colnames(mat),
    condition = rep(c('Control', 'Treatment'), each = 10),
    batch = rep(c('A', 'B'), 10),
    row.names = colnames(mat)
)

gene_info <- data.frame(
    gene = rownames(mat),
    pathway = sample(c('Metabolism', 'Signaling', 'Immune'), n_genes, replace = TRUE),
    log2FC = rnorm(n_genes, 0, 1.5),
    row.names = rownames(mat)
)

col_fun <- colorRamp2(c(-2, 0, 2), c('#4DBBD5', 'white', '#E64B35'))

ha_col <- HeatmapAnnotation(
    Condition = metadata$condition,
    Batch = metadata$batch,
    col = list(
        Condition = c(Control = '#4DBBD5', Treatment = '#E64B35'),
        Batch = c(A = '#00A087', B = '#3C5488')
    ),
    annotation_name_side = 'left'
)

ha_row <- rowAnnotation(
    Pathway = gene_info$pathway,
    LogFC = anno_barplot(gene_info$log2FC, baseline = 0,
                          gp = gpar(fill = ifelse(gene_info$log2FC > 0, '#E64B35', '#4DBBD5')),
                          width = unit(2, 'cm')),
    col = list(Pathway = c(Metabolism = '#8491B4', Signaling = '#91D1C2', Immune = '#F39B7F'))
)

ht <- Heatmap(mat,
              name = 'Z-score',
              col = col_fun,
              top_annotation = ha_col,
              left_annotation = ha_row,
              row_split = gene_info$pathway,
              column_split = metadata$condition,
              cluster_row_slices = FALSE,
              cluster_column_slices = FALSE,
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_rot = 45,
              row_title_rot = 0,
              column_title = 'Expression Heatmap',
              heatmap_legend_param = list(title = 'Z-score', direction = 'horizontal'))

pdf('expression_heatmap.pdf', width = 12, height = 10)
draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
dev.off()

message('Heatmap saved: expression_heatmap.pdf')
