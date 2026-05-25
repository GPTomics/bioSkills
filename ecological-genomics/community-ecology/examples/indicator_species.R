# Reference: indicspecies 1.7+, vegan 2.6+, ggplot2 3.5+ | Verify API if version differs
# Indicator species analysis using IndVal.g (group-size-equalized; required for
# unbalanced designs per De Caceres & Legendre 2009). Basic 'IndVal' is biased
# toward larger groups.
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

# Inject habitat-specific species to create realistic indicator patterns
species_matrix[1:10, 1:3] <- species_matrix[1:10, 1:3] + rpois(30, lambda = 20)
species_matrix[11:20, 4:6] <- species_matrix[11:20, 4:6] + rpois(30, lambda = 15)
species_matrix[21:30, 7:9] <- species_matrix[21:30, 7:9] + rpois(30, lambda = 18)
species_matrix[31:40, 10:12] <- species_matrix[31:40, 10:12] + rpois(30, lambda = 12)

# --- IndVal.g indicator species analysis (group-size-equalized) ---
# func='IndVal.g': the MODERN default per De Caceres & Legendre 2009 Ecology 90:3566
#   Basic 'IndVal' is biased toward larger groups; do NOT use for unbalanced designs
# duleg=TRUE: tests individual groups only (faster, simpler interpretation)
# duleg=FALSE: tests species against all group combinations
# 999 permutations: standard significance testing
mp <- multipatt(species_matrix, site_groups,
                func = 'IndVal.g',
                duleg = TRUE,
                control = how(nperm = 999))

cat('--- Indicator Species Analysis (IndVal.g, duleg=TRUE) ---\n')
summary(mp, alpha = 0.05)

# --- Extract significant indicators sorted by p-value ---
# For multi-comparison correction across many species, prefer FDR over Bonferroni
sig <- mp$sign[!is.na(mp$sign$p.value) & mp$sign$p.value < 0.05, ]
sig <- sig[order(sig$p.value), ]
cat('\nSignificant indicators (p < 0.05), sorted by p-value:\n')
print(sig)

# --- Strong vs moderate vs weak indicators ---
# stat = sqrt(IndVal.g): combined measure of specificity * fidelity
# stat > 0.7: strong; 0.5-0.7: moderate; < 0.5: weak
strong <- mp$sign[!is.na(mp$sign$p.value) & mp$sign$p.value < 0.05 &
                  mp$sign$stat > 0.7, ]
cat('\nStrong indicators (stat > 0.7):\n')
print(strong)

# --- Specificity (A) and Fidelity (B) decomposition ---
# IndVal = A * B
# A (specificity): proportion of individuals of species i found in target group
# B (fidelity): proportion of target-group sites where species i occurs
# High A + low B: dominant where present but rare overall; "concentration indicator"
# Low A + high B: widespread in target group but not exclusive; "frequency indicator"
# Strong indicators have BOTH high A and high B
cat('\n--- Specificity (A) and Fidelity (B) components ---\n')
summary(mp, indvalcomp = TRUE)

# --- Point-biserial correlation alternative ---
# func='r.g': group-equalized correlation coefficient
# Often easier to interpret as a standardized effect size
pb <- multipatt(species_matrix, site_groups,
                func = 'r.g',
                duleg = TRUE,
                control = how(nperm = 999))
cat('\n--- Point-Biserial Correlation (r.g, group-equalized) ---\n')
summary(pb, alpha = 0.05)

# --- Combination indicators ---
# duleg=FALSE: tests species against all possible group combinations
# Identifies species characteristic of MULTIPLE related habitats
# (e.g., "wetland or alpine" but not forest or grassland)
mp_combo <- multipatt(species_matrix, site_groups,
                      func = 'IndVal.g',
                      duleg = FALSE,
                      control = how(nperm = 999))
cat('\n--- Combination Indicators (duleg=FALSE) ---\n')
summary(mp_combo, indvalcomp = TRUE)
