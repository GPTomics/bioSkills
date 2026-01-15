# Create heatmap of differentially expressed genes

library(pheatmap)
library(RColorBrewer)

# Simulate expression data (in practice, use vst(dds) output)
set.seed(42)
n_genes <- 50
n_samples <- 6

# Simulated normalized, scaled expression matrix
mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes,
              dimnames = list(paste0('gene', 1:n_genes),
                             paste0('sample', 1:n_samples)))

# Add some pattern (first half of genes up in treated)
mat[1:25, 4:6] <- mat[1:25, 4:6] + 2
mat[26:50, 4:6] <- mat[26:50, 4:6] - 2

# Sample annotation
annotation_col <- data.frame(
    condition = factor(rep(c('control', 'treated'), each = 3)),
    row.names = colnames(mat)
)

# Color palette
colors <- colorRampPalette(rev(brewer.pal(n = 7, name = 'RdBu')))(100)

# Create heatmap
pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         color = colors,
         border_color = NA,
         main = 'Top Differentially Expressed Genes',
         filename = 'heatmap.png',
         width = 8,
         height = 10)

cat('Heatmap saved to heatmap.png\n')
