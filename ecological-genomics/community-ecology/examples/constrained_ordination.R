library(vegan)
library(ggplot2)

# --- Simulated ecological data ---
set.seed(42)
n_sites <- 30
n_species <- 25

env_data <- data.frame(
    temperature = rnorm(n_sites, mean = 15, sd = 5),
    precipitation = rnorm(n_sites, mean = 800, sd = 200),
    soil_pH = rnorm(n_sites, mean = 6.5, sd = 1),
    elevation = rnorm(n_sites, mean = 500, sd = 200),
    habitat = sample(c('forest', 'grassland', 'wetland'), n_sites, replace = TRUE)
)
rownames(env_data) <- paste0('site_', 1:n_sites)

species_matrix <- matrix(rpois(n_sites * n_species, lambda = 5),
                         nrow = n_sites, ncol = n_species)
rownames(species_matrix) <- paste0('site_', 1:n_sites)
colnames(species_matrix) <- paste0('sp_', 1:n_species)
species_matrix <- species_matrix[, colSums(species_matrix) > 0]

# --- Step 1: Check gradient length with DCA ---
dca <- decorana(species_matrix)
# Axis 1 gradient length > 3 SD: unimodal (CCA); < 3 SD: linear (RDA)
gradient_length <- diff(range(scores(dca, display = 'sites', choices = 1)))
cat('DCA axis 1 length:', gradient_length, 'SD\n')

# --- Step 2: CCA (for unimodal responses) ---
cca_result <- cca(species_matrix ~ temperature + precipitation + soil_pH + elevation,
                  data = env_data)

cat('\n--- CCA Results ---\n')
cat('Total inertia:', cca_result$tot.chi, '\n')
cat('Constrained inertia:', cca_result$CCA$tot.chi, '\n')
cat('Proportion explained:', round(cca_result$CCA$tot.chi / cca_result$tot.chi, 3), '\n')

# Permutation tests
cat('\nGlobal test:\n')
anova(cca_result, permutations = 999)

cat('\nMarginal tests (per variable):\n')
anova(cca_result, by = 'margin', permutations = 999)

# --- Step 3: RDA with Hellinger transformation ---
species_hell <- decostand(species_matrix, method = 'hellinger')

rda_full <- rda(species_hell ~ temperature + precipitation + soil_pH + elevation,
                data = env_data)
cat('\n--- RDA Results ---\n')
cat('Adjusted R-squared:', round(RsquareAdj(rda_full)$adj.r.squared, 3), '\n')

# --- Step 4: Forward selection ---
rda_null <- rda(species_hell ~ 1, data = env_data)
rda_sel <- ordiR2step(rda_null, scope = formula(rda_full),
                      direction = 'forward', permutations = 999)

cat('\nSelected model formula:\n')
print(rda_sel$call)
cat('Selected model adj R-squared:', round(RsquareAdj(rda_sel)$adj.r.squared, 3), '\n')

# VIF check: > 10 indicates problematic collinearity
cat('\nVIF for selected variables:\n')
print(vif.cca(rda_sel))

# --- Step 5: Variance partitioning ---
env_numeric <- env_data[, c('temperature', 'precipitation', 'soil_pH', 'elevation')]
spatial_coords <- data.frame(x = rnorm(n_sites), y = rnorm(n_sites))

vp <- varpart(species_hell, env_numeric, spatial_coords)
cat('\n--- Variance Partitioning ---\n')
cat('[a] Unique environment:', round(vp$part$fract$Adj.R.squared[1], 3), '\n')
cat('[b] Shared:', round(vp$part$fract$Adj.R.squared[2], 3), '\n')
cat('[c] Unique spatial:', round(vp$part$fract$Adj.R.squared[3], 3), '\n')
cat('[d] Residual:', round(vp$part$fract$Adj.R.squared[4], 3), '\n')

pdf('variance_partition.pdf', width = 6, height = 6)
plot(vp, digits = 2, bg = c('skyblue', 'tomato'))
dev.off()

# --- Step 6: Publication-quality CCA triplot ---
site_scores <- as.data.frame(scores(cca_result, display = 'sites', scaling = 2))
site_scores$habitat <- env_data$habitat

sp_scores <- as.data.frame(scores(cca_result, display = 'species', scaling = 2))
env_scores <- as.data.frame(scores(cca_result, display = 'bp', scaling = 2))
env_scores$variable <- rownames(env_scores)

# Arrow scaling factor for visibility
arrow_scale <- 3

p <- ggplot() +
    geom_point(data = site_scores, aes(x = CCA1, y = CCA2, color = habitat),
               size = 3, alpha = 0.8) +
    geom_text(data = sp_scores, aes(x = CCA1, y = CCA2, label = rownames(sp_scores)),
              size = 2, color = 'grey50', alpha = 0.6) +
    geom_segment(data = env_scores,
                 aes(x = 0, y = 0, xend = CCA1 * arrow_scale, yend = CCA2 * arrow_scale),
                 arrow = arrow(length = unit(0.2, 'cm')), color = 'red', linewidth = 0.8) +
    geom_text(data = env_scores,
              aes(x = CCA1 * arrow_scale * 1.15, y = CCA2 * arrow_scale * 1.15,
                  label = variable),
              color = 'red', size = 3.5, fontface = 'bold') +
    labs(x = paste0('CCA1 (', round(cca_result$CCA$eig[1] / cca_result$tot.chi * 100, 1), '%)'),
         y = paste0('CCA2 (', round(cca_result$CCA$eig[2] / cca_result$tot.chi * 100, 1), '%)'),
         title = 'CCA Triplot: Species-Environment Relationships') +
    theme_bw() +
    theme(legend.position = 'bottom')

ggsave('cca_triplot.pdf', p, width = 9, height = 8)
