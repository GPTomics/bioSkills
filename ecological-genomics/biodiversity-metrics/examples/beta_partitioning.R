library(betapart)
library(vegan)
library(ggplot2)

# --- Example community matrix: sites (rows) x species (columns) ---
set.seed(42)
n_sites <- 12
n_species <- 30
community_matrix <- matrix(rpois(n_sites * n_species, lambda = 3),
                           nrow = n_sites, ncol = n_species)
rownames(community_matrix) <- paste0('site_', LETTERS[1:n_sites])
colnames(community_matrix) <- paste0('sp_', 1:n_species)

# Remove empty species columns
community_matrix <- community_matrix[, colSums(community_matrix) > 0]

# --- Presence/absence conversion ---
pa_matrix <- ifelse(community_matrix > 0, 1, 0)

# --- Pairwise Sorensen decomposition ---
# beta.sim: turnover (species replacement between sites)
# beta.sne: nestedness (richness difference / subset pattern)
# beta.sor: total beta diversity (turnover + nestedness)
pair_sor <- beta.pair(pa_matrix, index.family = 'sorensen')

cat('Mean pairwise turnover (Sorensen):', mean(as.vector(pair_sor$beta.sim)), '\n')
cat('Mean pairwise nestedness (Sorensen):', mean(as.vector(pair_sor$beta.sne)), '\n')
cat('Mean pairwise total beta (Sorensen):', mean(as.vector(pair_sor$beta.sor)), '\n')

# --- Pairwise Jaccard decomposition ---
pair_jac <- beta.pair(pa_matrix, index.family = 'jaccard')

# --- Multi-site decomposition ---
# Summarizes beta diversity across all sites simultaneously
multi <- beta.multi(pa_matrix, index.family = 'sorensen')
cat('\nMulti-site turnover:', multi$beta.SIM, '\n')
cat('Multi-site nestedness:', multi$beta.SNE, '\n')
cat('Multi-site total:', multi$beta.SOR, '\n')
cat('Turnover proportion:', round(multi$beta.SIM / multi$beta.SOR, 3), '\n')

# --- Abundance-based beta (Bray-Curtis decomposition) ---
# beta.bray.bal: balanced variation in abundance (analogous to turnover)
# beta.bray.gra: abundance gradient (analogous to nestedness)
pair_abund <- beta.pair.abund(community_matrix, index.family = 'bray')

cat('\nMean balanced variation (Bray):', mean(as.vector(pair_abund$beta.bray.bal)), '\n')
cat('Mean abundance gradient (Bray):', mean(as.vector(pair_abund$beta.bray.gra)), '\n')

# --- Visualization: turnover vs nestedness proportions ---
beta_df <- data.frame(
    turnover = as.vector(pair_sor$beta.sim),
    nestedness = as.vector(pair_sor$beta.sne),
    total = as.vector(pair_sor$beta.sor)
)
beta_df$turn_prop <- ifelse(beta_df$total > 0,
                            beta_df$turnover / beta_df$total, 0.5)

p1 <- ggplot(beta_df, aes(x = total, y = turn_prop)) +
    geom_point(alpha = 0.4, size = 2) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey40') +
    labs(x = 'Total beta diversity (Sorensen)',
         y = 'Turnover proportion',
         title = 'Beta Diversity Decomposition') +
    annotate('text', x = max(beta_df$total) * 0.8, y = 0.85, label = 'Turnover-dominated') +
    annotate('text', x = max(beta_df$total) * 0.8, y = 0.15, label = 'Nestedness-dominated') +
    ylim(0, 1) +
    theme_bw()
ggsave('beta_decomposition.pdf', p1, width = 7, height = 6)

# --- Cluster dendrogram from turnover distances ---
hc_turn <- hclust(pair_sor$beta.sim, method = 'average')
pdf('turnover_dendrogram.pdf', width = 8, height = 5)
plot(hc_turn, main = 'Site Clustering by Turnover (Sorensen)',
     xlab = '', sub = '', ylab = 'Turnover dissimilarity')
dev.off()

# --- NMDS ordination of beta diversity components ---
nmds_total <- metaMDS(pair_sor$beta.sor, k = 2, trymax = 100)
# stress < 0.1: good representation; 0.1-0.2: acceptable; > 0.2: poor
cat('\nNMDS stress:', nmds_total$stress, '\n')

nmds_df <- data.frame(NMDS1 = nmds_total$points[, 1],
                       NMDS2 = nmds_total$points[, 2],
                       site = rownames(nmds_total$points))

p2 <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, label = site)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    labs(title = sprintf('NMDS of Total Beta Diversity (stress = %.3f)',
                         nmds_total$stress)) +
    theme_bw()
ggsave('beta_nmds.pdf', p2, width = 7, height = 6)
