library(indicspecies)
library(vegan)

# --- Simulated community data ---
set.seed(42)
n_sites <- 40
n_species <- 30

species_matrix <- matrix(rpois(n_sites * n_species, lambda = 3),
                         nrow = n_sites, ncol = n_species)
rownames(species_matrix) <- paste0('site_', 1:n_sites)
colnames(species_matrix) <- paste0('sp_', 1:n_species)

# Habitat groups (10 sites each)
site_groups <- factor(rep(c('forest', 'grassland', 'wetland', 'alpine'), each = 10))

# Inject some habitat-specific species to create realistic indicator patterns
species_matrix[1:10, 1:3] <- species_matrix[1:10, 1:3] + rpois(30, lambda = 20)
species_matrix[11:20, 4:6] <- species_matrix[11:20, 4:6] + rpois(30, lambda = 15)
species_matrix[21:30, 7:9] <- species_matrix[21:30, 7:9] + rpois(30, lambda = 18)
species_matrix[31:40, 10:12] <- species_matrix[31:40, 10:12] + rpois(30, lambda = 12)

# --- IndVal indicator species analysis ---
# func='IndVal.g': group-equalized IndVal; handles unbalanced group sizes
# duleg=TRUE: tests individual groups only (not combinations)
# 999 permutations for significance testing
mp <- multipatt(species_matrix, site_groups, func = 'IndVal.g',
                duleg = TRUE, control = how(nperm = 999))

cat('--- Indicator Species Analysis (IndVal.g) ---\n')
summary(mp)

# --- Extract significant indicators ---
# p.value < 0.05: standard significance threshold for indicator associations
sig_indicators <- mp$sign[mp$sign$p.value < 0.05, ]
sig_indicators <- sig_indicators[order(sig_indicators$p.value), ]

cat('\n--- Significant Indicators ---\n')
cat('Number of significant indicator species:', nrow(sig_indicators), '\n\n')
print(sig_indicators)

# --- Indicator strength (stat = sqrt(IndVal)) ---
# stat > 0.7: strong indicator; 0.5-0.7: moderate; < 0.5: weak
strong <- mp$sign[mp$sign$p.value < 0.05 & mp$sign$stat > 0.7, ]
cat('\nStrong indicators (stat > 0.7):\n')
print(strong)

# --- Point-biserial correlation alternative ---
# r.g provides correlation coefficients (easier to interpret effect size)
pb <- multipatt(species_matrix, site_groups, func = 'r.g',
                duleg = TRUE, control = how(nperm = 999))

cat('\n--- Point-Biserial Correlation ---\n')
summary(pb)

# --- Indicator species for group combinations ---
# duleg=FALSE: tests all possible group combinations
# Useful for finding species shared between habitat types
mp_combo <- multipatt(species_matrix, site_groups, func = 'IndVal.g',
                      duleg = FALSE, control = how(nperm = 999))

cat('\n--- Indicators for Group Combinations ---\n')
summary(mp_combo, indvalcomp = TRUE)

# --- Specificity and fidelity components ---
# IndVal = Specificity * Fidelity
# Specificity (A): proportion of individuals in target group
# Fidelity (B): proportion of target-group sites where species occurs
cat('\nSpecificity (A) and Fidelity (B) components:\n')
summary(mp, indvalcomp = TRUE)
