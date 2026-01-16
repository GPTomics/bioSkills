library(phyloseq)
library(ggplot2)
library(vegan)

abundance <- read.table('merged_abundance.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

species <- abundance[grepl('s__', rownames(abundance)) & !grepl('t__', rownames(abundance)), ]
rownames(species) <- sapply(strsplit(rownames(species), '\\|'), tail, 1)
rownames(species) <- gsub('s__', '', rownames(species))

otu <- otu_table(as.matrix(species), taxa_are_rows = TRUE)

sample_metadata <- data.frame(
    Sample = colnames(species),
    Group = rep(c('Control', 'Treatment'), length.out = ncol(species)),
    row.names = colnames(species)
)
samp <- sample_data(sample_metadata)

ps <- phyloseq(otu, samp)

top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:10]
ps_top <- prune_taxa(top_taxa, ps)

p1 <- plot_bar(ps_top, fill = 'taxa_names') +
    geom_bar(stat = 'identity', position = 'stack') +
    theme_minimal() +
    labs(x = 'Sample', y = 'Relative Abundance (%)', title = 'Species Composition') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
ggsave('stacked_bar_phyloseq.png', p1, width = 10, height = 6, dpi = 300)

ord <- ordinate(ps, method = 'PCoA', distance = 'bray')
p2 <- plot_ordination(ps, ord, color = 'Group') +
    geom_point(size = 4) +
    stat_ellipse() +
    theme_minimal() +
    labs(title = 'PCoA of Sample Composition (Bray-Curtis)')
ggsave('pcoa_phyloseq.png', p2, width = 8, height = 6, dpi = 300)

alpha_div <- estimate_richness(ps, measures = c('Shannon', 'Simpson', 'Observed'))
alpha_div$Group <- sample_data(ps)$Group

p3 <- ggplot(alpha_div, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, size = 2) +
    theme_minimal() +
    labs(title = 'Alpha Diversity by Group', y = 'Shannon Index')
ggsave('alpha_diversity.png', p3, width = 6, height = 5, dpi = 300)

cat('Visualizations saved: stacked_bar_phyloseq.png, pcoa_phyloseq.png, alpha_diversity.png\n')
