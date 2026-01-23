# BCR repertoire analysis with Immcantation

library(alakazam)
library(shazam)
library(dplyr)
library(ggplot2)

# Load AIRR-formatted data
# Required columns: sequence_id, sequence, v_call, j_call, junction, junction_aa
db <- readChangeoDb('bcr_data_airr.tsv')
cat('Loaded', nrow(db), 'sequences\n')

# Clonal clustering
# Groups sequences that likely derive from same ancestor B cell
# Threshold 0.15 = 15% nucleotide distance in junction region
# Lower threshold (0.10) for IgM, higher (0.20) for highly mutated IgG
library(scoper)
CLONE_THRESHOLD <- 0.15
db <- hierarchicalClones(
    db,
    threshold = CLONE_THRESHOLD,
    method = 'nt',
    linkage = 'single'
)

n_clones <- length(unique(db$clone_id))
cat('Identified', n_clones, 'clonal lineages\n')

# Clone size distribution
clone_sizes <- db %>%
    group_by(clone_id) %>%
    summarize(size = n()) %>%
    arrange(desc(size))

cat('Largest clone has', max(clone_sizes$size), 'sequences\n')

# Calculate somatic hypermutation frequencies
# Uses S5F mutation model (Smith et al. 2009)
# Compares sequence to germline to count mutations
db <- observedMutations(
    db,
    sequenceColumn = 'sequence_alignment',
    germlineColumn = 'germline_alignment_d_mask',
    regionDefinition = IMGT_V,  # IMGT V region definition
    mutationDefinition = MUTATION_SCHEMES$S5F
)

# Mutation frequency summary
# mu_freq_seq_r: Replacement mutation frequency
# mu_freq_seq_s: Silent mutation frequency
# Higher R mutations may indicate selection
mutation_summary <- db %>%
    summarize(
        mean_mu_r = mean(mu_freq_seq_r, na.rm = TRUE),
        mean_mu_s = mean(mu_freq_seq_s, na.rm = TRUE),
        median_mu_r = median(mu_freq_seq_r, na.rm = TRUE)
    )

cat('\nMutation summary:\n')
print(mutation_summary)

# Plot mutation frequency distribution
ggplot(db, aes(x = mu_freq_seq_r)) +
    geom_histogram(bins = 50, fill = 'steelblue', alpha = 0.7) +
    labs(
        x = 'Replacement Mutation Frequency',
        y = 'Count',
        title = 'Somatic Hypermutation Distribution'
    ) +
    theme_minimal()
ggsave('mutation_distribution.pdf', width = 8, height = 6)

# V gene usage
v_usage <- countGenes(db, gene = 'v_call', groups = NULL, mode = 'gene')
cat('\nTop V genes:\n')
print(head(v_usage, 10))

# Test for selection (if sufficient data)
# BASELINe compares observed R/S ratio to expected under neutrality
# Requires clones with multiple mutations
if (sum(db$mu_count_seq_r >= 5, na.rm = TRUE) >= 50) {
    cat('\nRunning selection analysis...\n')

    baseline <- estimateBaseline(
        db %>% filter(mu_count_seq_r >= 5),  # Need sufficient mutations
        sequenceColumn = 'sequence_alignment',
        germlineColumn = 'germline_alignment_d_mask',
        testStatistic = 'focused',
        regionDefinition = IMGT_V,
        nproc = 4
    )

    selection <- summarizeBaseline(baseline, returnType = 'df')
    cat('Selection analysis complete\n')
    # Positive sigma in CDR = positive selection (antigen-driven)
    # Negative sigma in FWR = negative selection (maintain structure)
}

cat('\nAnalysis complete!\n')
