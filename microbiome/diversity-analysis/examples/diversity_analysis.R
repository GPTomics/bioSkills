# Alpha and beta diversity analysis with phyloseq
library(phyloseq)
library(vegan)
library(ggplot2)

# Public microbiome datasets:
# - Qiita: https://qiita.ucsd.edu (16S/metagenomics studies)
# - Earth Microbiome Project: https://earthmicrobiome.org/data-and-code
# - HMP: https://hmpdacc.org (Human Microbiome Project)
# - Bioconductor: microbiomeDataSets package with curated phyloseq objects
# - MicrobiomeDB: https://microbiomedb.org
# Example: data('GlobalPatterns', package = 'phyloseq') for built-in demo

seqtab <- readRDS('seqtab_nochim.rds')
taxa <- readRDS('taxa.rds')
metadata <- read.csv('sample_metadata.csv', row.names = 1)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), tax_table(taxa), sample_data(metadata))
taxa_names(ps) <- paste0('ASV', seq(ntaxa(ps)))
cat('Samples:', nsamples(ps), '| ASVs:', ntaxa(ps), '\n')

# Rarefaction
min_depth <- min(sample_sums(ps))
cat('Minimum sample depth:', min_depth, '\n')
ps_rare <- rarefy_even_depth(ps, sample.size = min_depth, rngseed = 42)

# Alpha diversity
alpha <- estimate_richness(ps_rare, measures = c('Observed', 'Chao1', 'Shannon', 'Simpson'))
alpha <- cbind(alpha, sample_data(ps_rare))

cat('\nAlpha diversity by group:\n')
print(aggregate(Shannon ~ Group, data = alpha, FUN = function(x) c(mean = mean(x), sd = sd(x))))

# Statistical test
kw <- kruskal.test(Shannon ~ Group, data = alpha)
cat('\nKruskal-Wallis test (Shannon):\n')
cat(sprintf('  Chi-squared = %.2f, p = %.4f\n', kw$statistic, kw$p.value))

# Alpha diversity plot
p_alpha <- ggplot(alpha, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    theme_minimal() +
    labs(title = 'Shannon Diversity', y = 'Shannon Index') +
    theme(legend.position = 'none')
ggsave('alpha_diversity.pdf', p_alpha, width = 6, height = 5)

# Beta diversity (Bray-Curtis)
bray <- phyloseq::distance(ps_rare, method = 'bray')

# PERMANOVA
meta_df <- data.frame(sample_data(ps_rare))
# permutations=999: Standard for PERMANOVA. Use 9999 for publication; 99 for quick tests.
perm <- adonis2(bray ~ Group, data = meta_df, permutations = 999)
cat('\nPERMANOVA results:\n')
print(perm)

# PCoA ordination
pcoa <- ordinate(ps_rare, method = 'PCoA', distance = bray)
p_beta <- plot_ordination(ps_rare, pcoa, color = 'Group') +
    stat_ellipse(level = 0.95) +
    theme_minimal() +
    labs(title = sprintf('PCoA (Bray-Curtis)\nPERMANOVA R2=%.2f, p=%.3f',
                        perm$R2[1], perm$`Pr(>F)`[1]))
ggsave('beta_diversity_pcoa.pdf', p_beta, width = 7, height = 6)

cat('\nPlots saved: alpha_diversity.pdf, beta_diversity_pcoa.pdf\n')
