# Reference: GONE2, NeEstimator V2.1+, SNeP 1.1+, ggplot2 3.5+ | Verify API if version differs
# Ne estimation across time horizons:
#   - GONE2 for recent trajectory (cM column MUST be populated; not PLINK default 0)
#   - NeEstimator V2 for contemporary Ne (option-file API; line order matters)
#   - SNeP for physical-linkage-corrected LDNe on RAD-seq / WGS data
# Modern 100/1000 Ne rule (Frankham 2014); 50/500 was 1980-era convention.

library(ggplot2)

# --- GONE2: Recent Ne Trajectory ---
# CRITICAL pre-check: BIM/MAP cM column must be POPULATED (not all zeros)
# PLINK default uses physical position; GONE2 silently fails with cM=0
# Verify with: head genotypes.bim   -> column 3 should NOT be all zeros

# Run GONE2 from command line:
#   ./gone2 -t 4 -u 0.05 genotypes.vcf
# -t 4: threads
# -u 0.05: upper recombination rate bound (default; pairs with r > 0.05 excluded)
# Smaller -u (0.01) focuses on most recent generations; larger (0.1) extends back
# Minimum: 10,000 SNPs and 50 diploid individuals for reliable trajectory

# Parse GONE2 output
ne_df <- read.table('OUTPUT_GONE2', header = TRUE, sep = '\t')
cat('--- GONE2 Recent Ne Trajectory ---\n')
cat('Generations estimated:', nrow(ne_df), '\n')
cat('Most recent Ne:', ne_df$Ne[1], '\n')
cat('Oldest Ne in window:', ne_df$Ne[nrow(ne_df)], '\n')

# --- Modern 100/1000 thresholds (Frankham 2014; revises 1980s 50/500 convention) ---
current_ne <- ne_df$Ne[1]
if (current_ne < 100) {
    cat('\nWARNING: Ne <', 100, '- insufficient for short-term inbreeding protection\n')
    cat('         (Frankham 2014 revised 50/500 -> 100/1000)\n')
} else if (current_ne < 1000) {
    cat('\nCAUTION: Ne < 1000 - insufficient for long-term adaptive maintenance\n')
} else {
    cat('\nNe >= 1000 - adequate for long-term adaptive potential\n')
}

# --- Plot Ne trajectory ---
p1 <- ggplot(ne_df, aes(x = generation, y = Ne)) +
    geom_line(linewidth = 1.2, color = 'blue') +
    geom_hline(yintercept = 100, linetype = 'dashed', color = 'red') +
    geom_hline(yintercept = 1000, linetype = 'dashed', color = 'orange') +
    annotate('text', x = max(ne_df$generation) * 0.8, y = 100,
             label = 'Ne = 100 (Frankham 2014)',
             color = 'red', vjust = -0.5, size = 3) +
    annotate('text', x = max(ne_df$generation) * 0.8, y = 1000,
             label = 'Ne = 1000 (long-term)',
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
}

# --- NeEstimator V2 Option-File API ---
# NeEstimator V2 is option-file driven, NOT CLI-flag driven
# Line order in the .ne2 file matters; defaults can cause silent failures
# The INFO line in output confirms which methods actually ran

# Build .ne2 option file
option_file <- c(
    '1 0',                                 # 1=GenePop input format
    'input_genotypes.gen',                 # input file path
    '1 0 0 0',                             # methods: LD only; HetExcess/Coancestry/Temporal off
    '3',                                   # number of Pcrit cutoffs
    '0.05 0.02 0.01',                      # Pcrit values: 0.02 standard; 0.05 conservative
    '0',                                   # mating system: 0 = random mating; 1 = monogamy
    'output_results.txt'                   # output file
)
writeLines(option_file, 'ne_option.ne2')

# Run NeEstimator V2:
#   java -jar NeEstimator.jar ne_option.ne2
# Build from GitHub: bunop/NeEstimator2.X (JDK 1.8+ and Apache Ant)
# Verify INFO line in output_results.txt; confirms LD method ran

# --- SNeP: Physical-Linkage-Corrected LDNe for RAD-seq / WGS ---
# Naive LDNe overestimates drift (downward Ne bias) when SNPs are physically linked
# SNeP applies Waples & Do's chromosomal-linkage correction
# Multi-threaded; supports several corrections (sample size, mutation, phasing, recombination)

# Run from command line:
#   ./SNeP1.1 -ped genotypes.ped -map genotypes.map -threads 4 \
#             -mutationrate 1.4e-8 -out snep_results.txt
# Output: Ne estimates per recombination-rate bin

# --- Stairway Plot 2 blueprint preparation (for SFS-based deep Ne) ---
# Stairway Plot 2 works from the SFS (folded or unfolded) without parametric model
# Generate SFS first with easySFS or vcf2sfs from VCF
blueprint <- c(
    'popid: my_species',
    'nseq: 100',           # 2 * number of diploid individuals
    'L: 50000000',         # total callable sites (monomorphic + polymorphic)
    'whether_folded: true',
    'SFS: 5000 3000 2000 1500 1000 800 600 500 400 350',
    'smallest_size_of_SFS_bin_used_for_estimation: 1',
    'largest_size_of_SFS_bin_used_for_estimation: 49',
    'pct_training: 0.67',
    'nrand: 10 20 30 40',  # random break points to try
    'project_dir: stairway_output',
    'stairway_plot_dir: /path/to/stairway_plot_v2',
    'ninput: 200',         # bootstrap replicates
    'random_seed: 12345',
    # mu: per-generation per-site mutation rate (1.4e-8 typical vertebrate)
    'mu: 1.4e-8',
    'year_per_generation: 5'
)
writeLines(blueprint, 'stairway_blueprint.txt')
cat('\nStairway Plot 2 blueprint written. Run:\n')
cat('  java -cp stairway_plot_v2.jar Stairbuilder stairway_blueprint.txt\n')
cat('  bash stairway_blueprint.sh\n')

# --- Ne/Nc conversion warning ---
# Frankham 1995 reported Ne/Nc = 0.1 median in terrestrial vertebrates
# Hauser & Carvalho 2008 documented Ne/Nc spanning 2-6 orders of magnitude (10^-2 to 10^-6) in marine fish
# Do NOT default to 0.1 for high-fecundity / sweepstakes-reproduction species
cat('\n--- Ne/Nc conversion caveat ---\n')
cat('For marine fish with sweepstakes recruitment, Hauser & Carvalho 2008\n',
    'documented Ne/Nc spanning 2-6 orders of magnitude smaller than census\n',
    '(i.e., Ne/Nc ratio from 10^-2 to 10^-6 across species).\n',
    'For long-lived mammals/birds, 0.1-0.3 is typical (Frankham 1995 median).\n',
    'Cite Hauser & Carvalho 2008 Fish Fish 9:333-362 for taxonomic variation.\n',
    sep = '')

# --- Export results ---
write.csv(ne_df, 'gone2_ne_estimates.csv', row.names = FALSE)
