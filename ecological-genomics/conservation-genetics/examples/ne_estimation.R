library(ggplot2)

# --- GONE2: Recent Ne Trajectory from Linkage Disequilibrium ---
# GONE2 is a standalone CLI tool (esrud/GONE2), not an R package
# Requires PLINK bed/bim/fam or VCF files
# Minimum: 10,000 SNPs and 50 diploid individuals for reliable estimates

# Run GONE2 from command line first:
# ./gone2 -t 4 -u 0.05 genotypes.vcf
# -t 4: threads; -u 0.05: upper recombination rate bound (default; pairs with r > 0.05 excluded)
# Captures Ne changes over ~200 recent generations
# Smaller -u (0.01) focuses on the most recent generations
# Larger -u (0.1) extends further back but with less resolution

# Parse GONE2 output file (tab-separated: generation, Ne columns)
ne_df <- read.table('OUTPUT_GONE2', header = TRUE, sep = '\t')

cat('--- GONE2 Ne Trajectory ---\n')
cat('Generations estimated:', nrow(ne_df), '\n')
cat('Most recent Ne:', ne_df$Ne[1], '\n')
cat('Oldest Ne:', ne_df$Ne[nrow(ne_df)], '\n')

# --- Conservation threshold assessment ---
# Ne < 50: critical risk of inbreeding depression (Franklin 1980, 50/500 rule)
# Ne < 500: insufficient for long-term adaptive potential
# Ne > 500: genetically viable population
current_ne <- ne_df$Ne[1]
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

# --- NeEstimator LD Method ---
# For quick contemporary Ne estimate without the full GONE2 trajectory
# NeEstimator v2 uses an option file (.ne2), not command-line flags
# Build from GitHub: bunop/NeEstimator2.X (requires JDK 1.8+ and Apache Ant)
# Run: java -jar NeEstimator.jar option_file.ne2
#
# Key options in .ne2 file:
#   Input file: genepop or FSTAT format
#   Method: LD (linkage disequilibrium)
#   Pcrit: 0.02 (exclude alleles < 2% frequency)
#     Lower pcrit (0.01) includes more rare alleles but increases downward bias
#     Higher pcrit (0.05) reduces bias but loses information
#
# Output: tab-separated with Ne point estimate, 95% CI (jackknife + parametric)

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
