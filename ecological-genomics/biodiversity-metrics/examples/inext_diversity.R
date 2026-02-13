# Reference: ggplot2 3.5+, vegan 2.6+ | Verify API if version differs
library(iNEXT)
library(ggplot2)

# --- Example abundance data: species counts per site ---
abundance_data <- list(
    forest = c(150, 80, 45, 30, 20, 15, 10, 8, 5, 3, 2, 1, 1, 1),
    grassland = c(100, 90, 70, 50, 40, 30, 25, 20, 15, 10, 8, 5, 3, 2, 1),
    wetland = c(200, 50, 20, 10, 5, 3, 2, 1, 1)
)

# --- iNEXT coverage-based rarefaction/extrapolation ---
# q=c(0,1,2): richness, Shannon exponential, Simpson inverse
# nboot=200: bootstrap replicates for 95% confidence intervals
# endpoint: extrapolate to 2x observed sample size for each assemblage
result <- iNEXT(abundance_data, q = c(0, 1, 2), datatype = 'abundance',
                nboot = 200)

# --- Rarefaction/extrapolation curves (diversity vs sample size) ---
p1 <- ggiNEXT(result, type = 1) +
    theme_bw() +
    labs(title = 'Rarefaction/Extrapolation Curves')
ggsave('rarefaction_curves.pdf', p1, width = 10, height = 6)

# --- Sample completeness profiles (coverage vs sample size) ---
p2 <- ggiNEXT(result, type = 2) +
    theme_bw() +
    labs(title = 'Sample Completeness')
ggsave('sample_completeness.pdf', p2, width = 10, height = 6)

# --- Coverage-based rarefaction (diversity vs coverage) ---
# Preferred comparison method: standardizes by sampling completeness
p3 <- ggiNEXT(result, type = 3) +
    theme_bw() +
    labs(title = 'Coverage-Based Diversity Comparison')
ggsave('coverage_based_diversity.pdf', p3, width = 10, height = 6)

# --- Point estimates at standardized coverage ---
# coverage=0.95: 95% of individuals belong to detected species
# This level represents adequate-to-good sampling for most communities
est_95 <- estimateD(abundance_data, datatype = 'abundance',
                    base = 'coverage', level = 0.95)
cat('Diversity at 95% coverage:\n')
est_95

# --- Asymptotic diversity estimates ---
# Chao1 (q=0), Chao-Shannon (q=1), Chao-Simpson (q=2)
# Represents estimated true diversity if sampling were exhaustive
cat('\nAsymptotic diversity estimates:\n')
result$AsyEst

# --- Custom comparison at multiple coverage levels ---
coverages <- c(0.90, 0.95, 0.99)
for (cov in coverages) {
    est <- estimateD(abundance_data, datatype = 'abundance',
                     base = 'coverage', level = cov)
    cat(sprintf('\n--- Coverage = %.0f%% ---\n', cov * 100))
    print(est[est$Order.q == 0, c('Assemblage', 'qD', 'qD.LCL', 'qD.UCL')])
}

# --- Incidence-based analysis (detection/non-detection across replicates) ---
# Format: first element = number of replicates, remaining = incidence frequencies
incidence_data <- list(
    site_A = c(20, 18, 15, 12, 10, 8, 6, 4, 3, 2, 1, 1),  # 20 replicates
    site_B = c(15, 14, 12, 10, 8, 6, 5, 3, 2, 1)           # 15 replicates
)

result_inc <- iNEXT(incidence_data, q = c(0, 1, 2), datatype = 'incidence_freq',
                    nboot = 200)
ggiNEXT(result_inc, type = 3) + theme_bw()
