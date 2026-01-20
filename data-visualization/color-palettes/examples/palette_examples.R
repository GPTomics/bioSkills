library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)

# --- ALTERNATIVE: Use real expression data ---
# For realistic gene expression visualization, try these Bioconductor datasets:
#
# library(airway)
# data('airway')
# vst_data <- vst(airway)
# pca <- prcomp(t(assay(vst_data)))
# df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
#                  condition = colData(airway)$dex,
#                  cell_line = colData(airway)$cell)

# Simulated data mimicking PCA of expression data
set.seed(42)
n_samples <- 100

# Simulate clusters like you'd see in a real PCA (treatment vs control separation)
group <- sample(c('Control', 'DrugA', 'DrugB', 'DrugC', 'DrugD'), n_samples, replace = TRUE)
group_centers <- list(Control = c(0, 0), DrugA = c(3, 1), DrugB = c(-2, 2), DrugC = c(1, -3), DrugD = c(-3, -2))
df <- data.frame(
    x = sapply(group, function(g) rnorm(1, group_centers[[g]][1], 0.8)),
    y = sapply(group, function(g) rnorm(1, group_centers[[g]][2], 0.8)),
    value = rnorm(n_samples, 0, 1.5),  # Continuous value like log2FC
    group = factor(group)
)
names(df)[1:2] <- c('PC1', 'PC2')

p1 <- ggplot(df, aes(PC1, PC2, color = value)) +
    geom_point(size = 3) +
    scale_color_viridis_c(option = 'viridis') +
    labs(title = 'Viridis (Sequential)', x = 'PC1', y = 'PC2', color = 'log2FC') +
    theme_minimal()

p2 <- ggplot(df, aes(PC1, PC2, color = value)) +
    geom_point(size = 3) +
    scale_color_gradient2(low = '#4DBBD5', mid = 'white', high = '#E64B35', midpoint = 0) +
    labs(title = 'Custom Diverging', x = 'PC1', y = 'PC2', color = 'log2FC') +
    theme_minimal()

p3 <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = 'Set1') +
    labs(title = 'Set1 (Qualitative)', x = 'PC1', y = 'PC2', color = 'Treatment') +
    theme_minimal()

npg_colors <- c('#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F')
p4 <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = npg_colors) +
    labs(title = 'NPG-style Colors', x = 'PC1', y = 'PC2', color = 'Treatment') +
    theme_minimal()

combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = 'Color Palette Examples',
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = 'bold')))

ggsave('palette_examples.pdf', combined, width = 12, height = 10)
message('Saved: palette_examples.pdf')

cat('\nRecommended palettes:\n')
cat('Sequential: viridis, magma, plasma, Blues, YlOrRd\n')
cat('Diverging: RdBu, coolwarm, PuOr\n')
cat('Qualitative: Set1, Dark2, npg, tab10\n')
