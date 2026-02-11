library(GONE2)
library(ggplot2)

# --- GONE2: Recent Ne Trajectory from Linkage Disequilibrium ---
# Requires PLINK bed/bim/fam files
# Minimum: 10,000 SNPs and 50 diploid individuals for reliable estimates

# hc=0.05: maximum recombination distance in Morgans (~5 cM)
# Captures Ne changes over ~200 recent generations
# Smaller values (0.01) focus on the most recent generations
# Larger values (0.1) extend further back but with less resolution
gone_result <- gone('genotypes', hc = 0.05, num_threads = 4)

cat('--- GONE2 Ne Trajectory ---\n')
cat('Generations estimated:', length(gone_result$generation), '\n')
cat('Most recent Ne:', gone_result$Ne[1], '\n')
cat('Oldest Ne:', gone_result$Ne[length(gone_result$Ne)], '\n')

# Ne trajectory data
ne_df <- data.frame(generation = gone_result$generation, Ne = gone_result$Ne)

# --- Conservation threshold assessment ---
# Ne < 50: critical risk of inbreeding depression (Franklin 1980, 50/500 rule)
# Ne < 500: insufficient for long-term adaptive potential
# Ne > 500: genetically viable population
current_ne <- gone_result$Ne[1]
if (current_ne < 50) {
    cat('\nWARNING: Ne < 50 - critical inbreeding risk\n')
} else if (current_ne < 500) {
    cat('\nCAUTION: Ne < 500 - limited adaptive potential\n')
} else {
    cat('\nNe > 500 - adequate for long-term viability\n')
}

# --- Plot Ne trajectory ---
p1 <- ggplot(ne_df, aes(x = generation, y = Ne)) +
    geom_line(linewidth = 1.2, color = 'blue') +
    geom_hline(yintercept = 50, linetype = 'dashed', color = 'red') +
    geom_hline(yintercept = 500, linetype = 'dashed', color = 'orange') +
    annotate('text', x = max(ne_df$generation) * 0.8, y = 50, label = 'Ne = 50 (critical)',
             color = 'red', vjust = -0.5, size = 3) +
    annotate('text', x = max(ne_df$generation) * 0.8, y = 500, label = 'Ne = 500 (viable)',
             color = 'orange', vjust = -0.5, size = 3) +
    scale_y_log10() +
    labs(x = 'Generations ago', y = 'Effective population size (Ne)',
         title = 'Recent Ne Trajectory (GONE2)') +
    theme_bw()
ggsave('gone2_ne_trajectory.pdf', p1, width = 9, height = 6)

# --- Detect population decline ---
if (nrow(ne_df) >= 10) {
    recent_ne <- mean(ne_df$Ne[1:5])
    historical_ne <- mean(ne_df$Ne[(nrow(ne_df) - 4):nrow(ne_df)])
    decline_ratio <- recent_ne / historical_ne
    cat(sprintf('\nRecent/historical Ne ratio: %.2f\n', decline_ratio))
    if (decline_ratio < 0.5) {
        cat('Strong decline detected (>50% reduction)\n')
    } else if (decline_ratio < 0.8) {
        cat('Moderate decline detected (20-50% reduction)\n')
    } else {
        cat('No major decline detected\n')
    }
}

# --- NeEstimator LD Method (R interface) ---
# For quick contemporary Ne estimate without the full GONE2 trajectory
# This section shows how to prepare input and parse output

# Generate genepop format for NeEstimator
# NeEstimator CLI command:
# NeEstimator --input genotypes.gen --method LD --pcrit 0.02 --output ne_ld.txt
#
# pcrit=0.02: exclude alleles with frequency < 2%
# Lower pcrit (0.01) includes more rare alleles but increases downward bias
# Higher pcrit (0.05) reduces bias but loses information
#
# Parse NeEstimator output:
# Ne point estimate, 95% CI (jackknife), parametric CI

# --- Stairway Plot 2 blueprint preparation ---
# Generates the blueprint file from a site frequency spectrum (SFS)
# SFS can be computed from VCF using tools like easySFS or dadi

# Example blueprint content:
blueprint <- c(
    'popid: my_species',
    'nseq: 100',           # 2 * number of diploid individuals
    'L: 50000000',         # total callable sites (monomorphic + polymorphic)
    'whether_folded: true', # true if ancestral allele is unknown
    'SFS: 5000 3000 2000 1500 1000 800 600 500 400 350',  # folded SFS entries
    'smallest_size_of_SFS_bin_used_for_estimation: 1',
    'largest_size_of_SFS_bin_used_for_estimation: 49',
    'pct_training: 0.67',
    'nrand: 10 20 30 40',  # random break points to try
    'project_dir: stairway_output',
    'stairway_plot_dir: /path/to/stairway_plot_v2',
    'ninput: 200',         # bootstrap replicates
    'random_seed: 12345',
    # mutation_rate: per-generation per-site
    # 1.4e-8: typical vertebrate rate; adjust for target taxon
    'mu: 1.4e-8',
    # generation_time: years per generation
    'year_per_generation: 5'
)

writeLines(blueprint, 'stairway_blueprint.txt')
cat('\nStairway Plot 2 blueprint written to stairway_blueprint.txt\n')
cat('Run: java -cp stairway_plot_v2.jar Stairbuilder stairway_blueprint.txt\n')

# --- Export GONE2 results ---
write.csv(ne_df, 'gone2_ne_estimates.csv', row.names = FALSE)
cat('GONE2 results written to gone2_ne_estimates.csv\n')
