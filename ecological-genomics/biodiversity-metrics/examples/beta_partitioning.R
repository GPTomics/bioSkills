# Reference: betapart 1.6+, vegan 2.6+, ggplot2 3.5+ | Verify API if version differs
# Beta diversity partitioning into turnover and nestedness components.
# Reports the Baselga partition; flags that the partition is not unique
# (Podani/Carvalho alternative gives different ecological interpretation).
# Applies Hellinger transformation before any ordination.

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
community_matrix <- community_matrix[, colSums(community_matrix) > 0]

# --- Presence/absence conversion (required for Baselga partition) ---
pa_matrix <- ifelse(community_matrix > 0, 1, 0)

# --- BASELGA partition: Sorensen = turnover + nestedness ---
# Baselga 2010 Glob Ecol Biogeogr 19:134-143
# beta.sim: Simpson turnover (pure species replacement)
# beta.sne: nestedness component (richness-difference)
# beta.sor: total Sorensen beta diversity
# NOTE: this is ONE of multiple valid partitions; report alternative below
pair_sor <- beta.pair(pa_matrix, index.family = 'sorensen')

cat('--- Baselga Sorensen partition (pairwise) ---\n')
cat('Mean turnover (beta.sim):', round(mean(as.vector(pair_sor$beta.sim)), 3), '\n')
cat('Mean nestedness (beta.sne):', round(mean(as.vector(pair_sor$beta.sne)), 3), '\n')
cat('Mean total (beta.sor):', round(mean(as.vector(pair_sor$beta.sor)), 3), '\n')

# --- Jaccard family (alternative metric, same Baselga framework) ---
pair_jac <- beta.pair(pa_matrix, index.family = 'jaccard')

# --- Multi-site decomposition ---
multi <- beta.multi(pa_matrix, index.family = 'sorensen')
cat('\n--- Multi-site Baselga partition ---\n')
cat('Multi-site turnover (beta.SIM):', round(multi$beta.SIM, 3), '\n')
cat('Multi-site nestedness (beta.SNE):', round(multi$beta.SNE, 3), '\n')
cat('Multi-site total (beta.SOR):', round(multi$beta.SOR, 3), '\n')
cat('Turnover proportion:', round(multi$beta.SIM / multi$beta.SOR, 3), '\n')

# --- Abundance-based decomposition (Bray-Curtis) ---
# beta.bray.bal: balanced variation (analogous to turnover)
# beta.bray.gra: abundance gradient (analogous to nestedness)
# NOTE: Bray-Curtis is NOT a true metric; use Sorensen if downstream method
# requires triangle inequality (e.g., Ward clustering)
pair_abund <- beta.pair.abund(community_matrix, index.family = 'bray')
cat('\n--- Abundance-based (Bray) ---\n')
cat('Mean balanced variation:', round(mean(as.vector(pair_abund$beta.bray.bal)), 3), '\n')
cat('Mean abundance gradient:', round(mean(as.vector(pair_abund$beta.bray.gra)), 3), '\n')

# --- DOCUMENT the partition-non-uniqueness caveat ---
# The Baselga partition is not the only valid decomposition.
# Podani & Schmera 2011 Oikos 120:1625-1638 proposed an alternative
# (richness-difference + replacement) that can give different conclusions.
# For postdoc-grade reporting: compute BOTH and report which dominates in each.
# Carvalho et al. 2012 betapart::beta.pair.abund() with index.family='ruzicka' gives
# yet another framework. The choice is ecological, not just statistical.
cat('\n--- Partition non-uniqueness flag ---\n')
cat('Baselga partition reported above; for sensitivity, also compute',
    'Podani/Carvalho framework (e.g., adespatial::beta.div.comp(...)).\n',
    'The two frameworks can disagree on which component dominates.\n', sep = '')

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
         y = 'Turnover proportion (Baselga partition)',
         title = 'Beta Diversity Decomposition (Baselga)') +
    annotate('text', x = max(beta_df$total) * 0.8, y = 0.85,
             label = 'Turnover-dominated') +
    annotate('text', x = max(beta_df$total) * 0.8, y = 0.15,
             label = 'Nestedness-dominated') +
    ylim(0, 1) +
    theme_bw()
ggsave('beta_decomposition.pdf', p1, width = 7, height = 6)

# --- Hellinger-transformed PCA (THE postdoc-grade community ordination) ---
# Legendre & Gallagher 2001 Oecologia 129:271-280
# Hellinger transform solves the double-zero problem:
# raw PCA on counts is dominated by sample-size effects, not biology
species_hell <- decostand(community_matrix, method = 'hellinger')
pca_hell <- rda(species_hell)
hell_scores <- scores(pca_hell, display = 'sites', scaling = 2)
hell_df <- data.frame(
    PC1 = hell_scores[, 1],
    PC2 = hell_scores[, 2],
    site = rownames(hell_scores)
)
p_pca <- ggplot(hell_df, aes(x = PC1, y = PC2, label = site)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    labs(title = 'Hellinger-Transformed PCA',
         subtitle = 'Hellinger transform mandatory; raw-count PCA is malpractice') +
    theme_bw()
ggsave('hellinger_pca.pdf', p_pca, width = 7, height = 6)

# --- NMDS ordination of total beta diversity ---
# Use Sorensen (a true metric) rather than Bray-Curtis (not a metric)
# for NMDS when downstream interpretation requires triangle inequality
nmds_total <- metaMDS(pair_sor$beta.sor, k = 2, trymax = 100)
# stress < 0.1: good representation; 0.1-0.2: acceptable; > 0.2: poor
cat('\nNMDS stress (Sorensen distances):', round(nmds_total$stress, 3), '\n')

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
