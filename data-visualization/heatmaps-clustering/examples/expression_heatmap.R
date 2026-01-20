library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# --- ALTERNATIVE: Use real Bioconductor datasets ---
# For realistic expression heatmaps, try these publicly available datasets:
#
# library(airway)
# data('airway')
# vst_data <- assay(vst(airway))
# top_var <- head(order(rowVars(vst_data), decreasing = TRUE), 100)
# mat <- t(scale(t(vst_data[top_var, ])))  # Z-score normalize
# metadata <- as.data.frame(colData(airway))
#
# Or use the ALL leukemia dataset:
# library(ALL)
# data('ALL')
# mat <- exprs(ALL)[1:100, ]
# metadata <- pData(ALL)

# Simulated data with realistic expression patterns
set.seed(42)
n_genes <- 100
n_samples <- 20

# Create base expression with biological structure
mat <- matrix(rnorm(n_genes * n_samples, 0, 0.5), nrow = n_genes)

# Add treatment effect for subset of genes (simulating DE genes)
treatment_responsive <- 1:30
mat[treatment_responsive, 11:20] <- mat[treatment_responsive, 11:20] + rnorm(30 * 10, 1.5, 0.3)

# Add pathway-correlated expression patterns
immune_genes <- 31:50
signaling_genes <- 51:70
mat[immune_genes, ] <- mat[immune_genes, ] + matrix(rep(rnorm(n_samples, 0, 0.8), 20), nrow = 20, byrow = TRUE)
mat[signaling_genes, ] <- mat[signaling_genes, ] + matrix(rep(rnorm(n_samples, 0, 0.6), 20), nrow = 20, byrow = TRUE)

rownames(mat) <- paste0('Gene', 1:n_genes)
colnames(mat) <- paste0('Sample', 1:n_samples)

metadata <- data.frame(
    sample = colnames(mat),
    condition = rep(c('Control', 'Treatment'), each = 10),
    batch = rep(c('A', 'B'), 10),
    row.names = colnames(mat)
)

pathway_assignment <- c(rep('Immune', 20), rep('Signaling', 20), rep('Metabolism', 60))
gene_info <- data.frame(
    gene = rownames(mat),
    pathway = pathway_assignment[1:n_genes],
    log2FC = c(rnorm(30, 1.2, 0.5), rnorm(70, 0, 0.8)),  # DE genes have positive log2FC
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
