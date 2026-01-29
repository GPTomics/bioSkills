library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)

# --- ALTERNATIVE: Use real data ---
# Download VCFs from TCGA/GDC:
# https://portal.gdc.cancer.gov/
#
# Or use PCAWG mutational signatures paper data

# For demonstration, simulate mutation data
set.seed(42)

# Simulate 96-context mutation matrix for 5 samples
# In practice, use read_vcfs_as_granges() to load VCF files
n_samples <- 5
n_contexts <- 96

message('Simulating mutation matrix for demonstration...')
message('For real data, load VCF files with read_vcfs_as_granges()')

# Get context names from MutationalPatterns
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
contexts <- c(
    'A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C',
    # ... full 96 contexts would be here
    # Using placeholder for demonstration
    paste0('context', 1:96)
)

# Simulate counts - mixture of signatures
# Sample 1-2: APOBEC-like (SBS2/13)
# Sample 3-4: UV-like (SBS7)
# Sample 5: Age-related (SBS1/5)

sim_mat <- matrix(
    c(
        rpois(96, lambda = 20),  # Sample 1
        rpois(96, lambda = 25),  # Sample 2
        rpois(96, lambda = 30),  # Sample 3
        rpois(96, lambda = 15),  # Sample 4
        rpois(96, lambda = 10)   # Sample 5
    ),
    nrow = 96,
    ncol = 5
)
colnames(sim_mat) <- paste0('Sample', 1:5)
rownames(sim_mat) <- contexts

# --- Real workflow starts here ---

# Load COSMIC signatures (version 3.2)
# Contains 79 single base substitution signatures
cosmic_sigs <- get_known_signatures(muttype = 'snv')
message(sprintf('Loaded %d COSMIC signatures', ncol(cosmic_sigs)))

# Fit samples to COSMIC signatures
# Uses non-negative least squares (NNLS)
message('Fitting samples to COSMIC signatures...')
fit_result <- fit_to_signatures(sim_mat, cosmic_sigs)

# Extract contributions
contributions <- fit_result$contribution
reconstructed <- fit_result$reconstructed

# Calculate cosine similarity between original and reconstructed
cos_sim_samples <- diag(cos_sim_matrix(sim_mat, reconstructed))
message(sprintf('Mean reconstruction accuracy (cosine): %.3f', mean(cos_sim_samples)))

# Get top signatures per sample
message('\nTop signatures per sample:')
for (sample in colnames(contributions)) {
    sample_contrib <- contributions[, sample]
    total <- sum(sample_contrib)

    # Filter to signatures with >5% contribution
    significant <- sample_contrib[sample_contrib / total > 0.05]
    significant <- sort(significant, decreasing = TRUE)

    message(sprintf('\n%s (total mutations: %d):', sample, total))
    for (sig in names(significant)[1:min(3, length(significant))]) {
        pct <- significant[sig] / total * 100
        message(sprintf('  %s: %.1f%%', sig, pct))
    }
}

# Signature interpretation reference
signature_etiology <- c(
    'SBS1' = 'Spontaneous deamination (age-related)',
    'SBS2' = 'APOBEC activity',
    'SBS3' = 'HR deficiency (BRCA1/2)',
    'SBS4' = 'Tobacco smoking',
    'SBS5' = 'Unknown (age-related)',
    'SBS6' = 'MMR deficiency',
    'SBS7a' = 'UV exposure',
    'SBS7b' = 'UV exposure',
    'SBS13' = 'APOBEC activity',
    'SBS10a' = 'POLE mutation',
    'SBS10b' = 'POLE mutation'
)

message('\n=== Signature Etiology Reference ===')
for (sig in names(signature_etiology)) {
    message(sprintf('%s: %s', sig, signature_etiology[sig]))
}

# --- Visualization ---
# Uncomment to generate plots

# Plot 96-profile (mutation spectrum)
# pdf('mutation_spectrum.pdf', width = 12, height = 8)
# plot_96_profile(sim_mat)
# dev.off()

# Plot signature contribution
# pdf('signature_contributions.pdf', width = 10, height = 6)
# plot_contribution(contributions, cosmic_sigs, mode = 'relative')
# dev.off()

# Cosine similarity heatmap
# pdf('cosine_similarity.pdf', width = 10, height = 8)
# cos_sim <- cos_sim_matrix(sim_mat, cosmic_sigs)
# plot_cosine_heatmap(cos_sim)
# dev.off()

message('\nAnalysis complete.')
message('For clinical interpretation, check for:')
message('  - SBS3: Consider BRCA testing and PARP inhibitors')
message('  - SBS6/15/26/44: Consider MSI testing and immunotherapy')
message('  - SBS4: Evidence of tobacco exposure')
