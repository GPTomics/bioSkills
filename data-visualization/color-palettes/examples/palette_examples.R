library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)

set.seed(42)
df <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    value = rnorm(100),
    group = sample(LETTERS[1:5], 100, replace = TRUE)
)

p1 <- ggplot(df, aes(x, y, color = value)) +
    geom_point(size = 3) +
    scale_color_viridis_c(option = 'viridis') +
    labs(title = 'Viridis (Sequential)') +
    theme_minimal()

p2 <- ggplot(df, aes(x, y, color = value)) +
    geom_point(size = 3) +
    scale_color_gradient2(low = '#4DBBD5', mid = 'white', high = '#E64B35', midpoint = 0) +
    labs(title = 'Custom Diverging') +
    theme_minimal()

p3 <- ggplot(df, aes(x, y, color = group)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = 'Set1') +
    labs(title = 'Set1 (Qualitative)') +
    theme_minimal()

npg_colors <- c('#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F')
p4 <- ggplot(df, aes(x, y, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = npg_colors) +
    labs(title = 'NPG-style Colors') +
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
