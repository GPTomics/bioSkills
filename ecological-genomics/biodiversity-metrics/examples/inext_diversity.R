# Reference: iNEXT 3.0+, ggplot2 3.5+, vegan 2.6+ | Verify API if version differs
# Hill-number diversity with coverage-based rarefaction/extrapolation.
# Reports effective species counts (numbers-equivalents), NOT raw entropies.
# Bounds extrapolation at 2x reference sample size (Chao et al. 2014 doubling rule).

library(iNEXT)
library(ggplot2)
library(vegan)

# --- Example abundance data: species counts per site ---
abundance_data <- list(
    forest = c(150, 80, 45, 30, 20, 15, 10, 8, 5, 3, 2, 1, 1, 1),
    grassland = c(100, 90, 70, 50, 40, 30, 25, 20, 15, 10, 8, 5, 3, 2, 1),
    wetland = c(200, 50, 20, 10, 5, 3, 2, 1, 1)
)

# --- Singleton/doubleton diagnostic BEFORE computing Chao1 ---
# Chao1 assumes singletons reflect undersampling of rare species
# In amplicon/eDNA data, singletons are dominated by PCR error -> Chao1 is unreliable
# In rigorous specimen-based data, singletons are usually real biology -> Chao1 is informative
f1_obs <- sapply(abundance_data, function(x) sum(x == 1))
f2_obs <- sapply(abundance_data, function(x) sum(x == 2))
cat('Singletons (f1):', f1_obs, '\n')
cat('Doubletons (f2):', f2_obs, '\n')
cat('Chao1 bound term f1^2/(2*f2):',
    f1_obs^2 / (2 * pmax(f2_obs, 0.5)), '\n')

# --- Hill numbers q=0,1,2 with coverage-based rarefaction-extrapolation ---
# q=0 = species richness (rare-species sensitive)
# q=1 = Shannon-equivalent exp(H) = effective species count (geometric mean weighting)
# q=2 = Simpson-equivalent 1/D = effective dominant-species count
# nboot=200: publication-quality bootstrap CIs (default 50 is too few)
result <- iNEXT(abundance_data, q = c(0, 1, 2), datatype = 'abundance',
                nboot = 200)

# --- Rarefaction/extrapolation curves (diversity vs sample size) ---
# Default endpoint = 2*max(sample size); enforces the doubling rule
p1 <- ggiNEXT(result, type = 1) +
    theme_bw() +
    labs(title = 'Rarefaction/Extrapolation Curves (Hill q=0,1,2)')
ggsave('rarefaction_curves.pdf', p1, width = 10, height = 6)

# --- Sample completeness profile (Good's coverage vs sample size) ---
# Shows whether the sample is approaching saturation
p2 <- ggiNEXT(result, type = 2) +
    theme_bw() +
    labs(title = 'Sample Completeness (Good\'s Coverage)')
ggsave('sample_completeness.pdf', p2, width = 10, height = 6)

# --- Coverage-based rarefaction (diversity vs coverage) ---
# THE postdoc-grade cross-site comparison plot
# Sample-size rarefaction is biased when assemblages differ in true diversity
p3 <- ggiNEXT(result, type = 3) +
    theme_bw() +
    labs(title = 'Coverage-Based Hill-Number Diversity')
ggsave('coverage_based_diversity.pdf', p3, width = 10, height = 6)

# --- Point estimates at standardized coverage ---
# coverage=0.95: postdoc-grade default for cross-site comparison
# 0.99: high-precision; 0.85: minimum to accept a site into comparison
est_95 <- estimateD(abundance_data, datatype = 'abundance',
                    base = 'coverage', level = 0.95)
cat('\nDiversity at C=0.95 (effective species counts):\n')
est_95

# --- Asymptotic estimates ---
# Chao1 (q=0), Chao-Shannon (q=1), Chao-Simpson (q=2) WITH 95% CIs
# CRITICAL: report as LOWER BOUND ('at least X species'), not point estimate
cat('\nAsymptotic diversity (LOWER BOUNDS):\n')
result$AsyEst

# --- Diagnostic: is coverage adequate to trust the bound? ---
# Coverage < 0.85 means Chao1 CI will be very wide and the bound uninformative
asy_with_coverage <- merge(result$AsyEst,
                            est_95[est_95$Order.q == 0, c('Assemblage', 'SC')],
                            by = 'Assemblage')
cat('\nAsymptotic estimates with sample coverage:\n')
asy_with_coverage

# --- Multiple coverage levels for sensitivity ---
coverages <- c(0.85, 0.95, 0.99)
for (cov in coverages) {
    est <- estimateD(abundance_data, datatype = 'abundance',
                     base = 'coverage', level = cov)
    cat(sprintf('\n--- Coverage = %.0f%% ---\n', cov * 100))
    print(est[est$Order.q == 0, c('Assemblage', 'qD', 'qD.LCL', 'qD.UCL')])
}

# --- Incidence-based analysis for amplicon/replicated data ---
# Format: first element = number of replicates, remaining = incidence frequencies
# Chao2 is the incidence equivalent of Chao1 and is more robust to PCR-error singletons
incidence_data <- list(
    site_A = c(20, 18, 15, 12, 10, 8, 6, 4, 3, 2, 1, 1),  # 20 replicates
    site_B = c(15, 14, 12, 10, 8, 6, 5, 3, 2, 1)           # 15 replicates
)

result_inc <- iNEXT(incidence_data, q = c(0, 1, 2),
                    datatype = 'incidence_freq', nboot = 200)
ggiNEXT(result_inc, type = 3) + theme_bw() +
    labs(title = 'Incidence-Based Coverage Comparison')

# --- vegan equivalents for completeness ---
# Convert abundance_data list to a community matrix
species_names <- paste0('sp', 1:max(sapply(abundance_data, length)))
community_matrix <- t(sapply(abundance_data, function(x) {
    out <- rep(0, length(species_names))
    out[seq_along(x)] <- x
    out
}))
colnames(community_matrix) <- species_names

# Numbers-equivalents (Hill numbers from vegan)
richness <- specnumber(community_matrix)             # q=0
shannon_eff <- exp(diversity(community_matrix, 'shannon'))  # q=1 effective
invsimpson <- diversity(community_matrix, 'invsimpson')     # q=2 effective
cat('\nvegan-equivalent Hill numbers (effective species counts):\n')
data.frame(q0 = richness, q1 = round(shannon_eff, 2),
           q2 = round(invsimpson, 2))
