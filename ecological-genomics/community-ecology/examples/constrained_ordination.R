# Reference: vegan 2.6+, ggplot2 3.5+ | Verify API if version differs
# Constrained ordination workflow: DCA gradient-length check -> CCA or RDA decision,
# forward selection with adjusted R-squared, variance partitioning, triplot.
# PERMANOVA is paired with PERMDISP (Anderson & Walsh 2013).
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

# --- Step 1: DCA gradient-length decision ---
# > 3 SD: unimodal (CCA); < 3 SD: linear (RDA on Hellinger-transformed); 2-3 SD: gray zone
# decorana prints the canonical "Axis lengths" line in SD units; read off the printout
dca <- decorana(species_matrix)
print(dca)  # Inspect "Axis lengths:" row for the canonical SD-unit gradient lengths
# Approximate programmatic access via the range of detrended axis-1 scores
gradient_length_approx <- diff(range(scores(dca, display = 'sites', choices = 1)))
cat('DCA axis 1 gradient length (approx, SD):', round(gradient_length_approx, 2), '\n')

# --- Step 2: CCA for unimodal responses ---
cca_result <- cca(species_matrix ~ temperature + precipitation + soil_pH + elevation,
                  data = env_data)
cat('\n--- CCA Results ---\n')
cat('Constrained inertia / total:',
    round(cca_result$CCA$tot.chi / cca_result$tot.chi, 3), '\n')
print(anova(cca_result, by = 'margin', permutations = 999))

# --- Step 3: RDA with mandatory Hellinger transformation ---
# Hellinger transform (Legendre & Gallagher 2001) is non-negotiable for community data
species_hell <- decostand(species_matrix, method = 'hellinger')
rda_full <- rda(species_hell ~ temperature + precipitation + soil_pH + elevation,
                data = env_data)
cat('\n--- RDA Results (Hellinger-transformed) ---\n')
# Adjusted R^2 (Peres-Neto 2006 Ecology 87:2614); raw R^2 is biased upward
cat('Adjusted R^2:', round(RsquareAdj(rda_full)$adj.r.squared, 3), '\n')

# --- Step 4: Forward selection with adjusted-R^2 criterion ---
rda_null <- rda(species_hell ~ 1, data = env_data)
rda_sel <- ordiR2step(rda_null, scope = formula(rda_full),
                      direction = 'forward', permutations = 999)
cat('\nSelected model:\n')
print(rda_sel$call)
# VIF > 10 indicates problematic multicollinearity
cat('\nVIF for selected variables:\n')
print(vif.cca(rda_sel))

# --- Step 5: PERMANOVA + PERMDISP (Anderson & Walsh 2013) ---
# Distance matrix (use Bray-Curtis on RAW abundances OR Hellinger-distance on transformed)
bray_dist <- vegdist(species_matrix, method = 'bray')

# adonis2 is the modern PERMANOVA; adonis() is deprecated
# by='margin' is correct for unbalanced designs
permanova <- adonis2(bray_dist ~ habitat + soil_pH, data = env_data,
                     by = 'margin', permutations = 999)
cat('\n--- PERMANOVA (adonis2, by=margin) ---\n')
print(permanova)

# MANDATORY companion: PERMDISP via betadisper
# Tests dispersion homogeneity; without this, PERMANOVA significance is confounded
disp <- betadisper(bray_dist, env_data$habitat)
disp_test <- permutest(disp, permutations = 999)
cat('\n--- PERMDISP (betadisper) ---\n')
print(disp_test)

# Interpretation rule (Anderson & Walsh 2013):
# PERMANOVA p<0.05 AND betadisper p>0.05 -> location difference is real
# PERMANOVA p<0.05 AND betadisper p<0.05 -> CANNOT distinguish location from dispersion
permanova_sig <- permanova[['Pr(>F)']][1] < 0.05
disp_sig <- disp_test$tab[['Pr(>F)']][1] < 0.05
if (permanova_sig && disp_sig) {
    cat('\nWARNING: PERMANOVA and PERMDISP both significant.\n')
    cat('Location-vs-dispersion confounded; cannot conclude pure centroid difference.\n')
}

# --- Step 6: Variance partitioning with adjusted R^2 ---
env_numeric <- env_data[, c('temperature', 'precipitation', 'soil_pH', 'elevation')]
spatial_coords <- data.frame(x = rnorm(n_sites), y = rnorm(n_sites))
vp <- varpart(species_hell, env_numeric, spatial_coords)
cat('\n--- Variance Partitioning (Adjusted R^2) ---\n')
print(vp$part$fract)
pdf('variance_partition.pdf', width = 6, height = 6)
plot(vp, digits = 2, bg = c('skyblue', 'tomato'))
dev.off()

# --- Step 7: Publication-quality CCA triplot ---
site_scores <- as.data.frame(scores(cca_result, display = 'sites', scaling = 2))
site_scores$habitat <- env_data$habitat
sp_scores <- as.data.frame(scores(cca_result, display = 'species', scaling = 2))
env_scores <- as.data.frame(scores(cca_result, display = 'bp', scaling = 2))
env_scores$variable <- rownames(env_scores)
arrow_scale <- 3  # for visibility in the plot

p <- ggplot() +
    geom_point(data = site_scores, aes(x = CCA1, y = CCA2, color = habitat),
               size = 3, alpha = 0.8) +
    geom_text(data = sp_scores, aes(x = CCA1, y = CCA2,
                                     label = rownames(sp_scores)),
              size = 2, color = 'grey50', alpha = 0.6) +
    geom_segment(data = env_scores,
                 aes(x = 0, y = 0,
                     xend = CCA1 * arrow_scale, yend = CCA2 * arrow_scale),
                 arrow = arrow(length = unit(0.2, 'cm')),
                 color = 'red', linewidth = 0.8) +
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
