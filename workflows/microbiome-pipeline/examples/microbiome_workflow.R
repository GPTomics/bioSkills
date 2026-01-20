# Complete 16S microbiome workflow: FASTQ to differential abundance
library(dada2)
library(phyloseq)
library(ALDEx2)
library(vegan)
library(ggplot2)

# === CONFIGURATION ===
path <- 'raw_reads'
silva_train <- 'silva_nr99_v138.1_train_set.fa.gz'
silva_species <- 'silva_species_assignment_v138.1.fa.gz'
metadata_file <- 'sample_metadata.csv'
output_dir <- 'microbiome_results'
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, 'plots'), showWarnings = FALSE)

cat('=== Microbiome Pipeline ===\n')

# === 1. READ FILES ===
cat('\n1. Reading input files...\n')
fnFs <- sort(list.files(path, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = '_R2_001.fastq.gz', full.names = TRUE))
sample_names <- sapply(strsplit(basename(fnFs), '_'), `[`, 1)
cat('  Samples:', length(fnFs), '\n')

# Quality profiles
pdf(file.path(output_dir, 'plots/quality_profiles.pdf'))
plotQualityProfile(fnFs[1:min(2, length(fnFs))])
plotQualityProfile(fnRs[1:min(2, length(fnRs))])
dev.off()

# === 2. FILTER & TRIM ===
cat('\n2. Filtering and trimming...\n')
dir.create('filtered', showWarnings = FALSE)
filtFs <- file.path('filtered', paste0(sample_names, '_F_filt.fastq.gz'))
filtRs <- file.path('filtered', paste0(sample_names, '_R_filt.fastq.gz'))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
cat('  Reads passing filter:', round(100 * sum(out[, 2]) / sum(out[, 1]), 1), '%\n')

# === 3. LEARN ERRORS & DENOISE ===
cat('\n3. Learning errors and denoising...\n')
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# === 4. MERGE & CHIMERAS ===
cat('\n4. Merging pairs and removing chimeras...\n')
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = FALSE)
seqtab <- makeSequenceTable(mergers)
cat('  ASVs before chimera removal:', ncol(seqtab), '\n')

seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE, verbose = FALSE)
cat('  ASVs after chimera removal:', ncol(seqtab_nochim), '\n')
cat('  Reads retained:', round(100 * sum(seqtab_nochim) / sum(seqtab), 1), '%\n')

# Track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochim))
colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim')
write.csv(track, file.path(output_dir, 'read_tracking.csv'))

# === 5. ASSIGN TAXONOMY ===
cat('\n5. Assigning taxonomy...\n')
taxa <- assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE)
if (file.exists(silva_species)) taxa <- addSpecies(taxa, silva_species)
cat('  Genus assignment rate:', round(100 * sum(!is.na(taxa[, 'Genus'])) / nrow(taxa), 1), '%\n')

# === 6. CREATE PHYLOSEQ ===
cat('\n6. Creating phyloseq object...\n')
metadata <- read.csv(metadata_file, row.names = 1)
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), tax_table(taxa), sample_data(metadata))
taxa_names(ps) <- paste0('ASV', seq(ntaxa(ps)))
cat('  Final ASVs:', ntaxa(ps), '\n')
cat('  Final samples:', nsamples(ps), '\n')
saveRDS(ps, file.path(output_dir, 'phyloseq_object.rds'))

# === 7. DIVERSITY ===
cat('\n7. Calculating diversity...\n')
ps_rare <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 42, verbose = FALSE)

alpha_div <- estimate_richness(ps_rare, measures = c('Observed', 'Shannon', 'Simpson'))
alpha_div <- cbind(alpha_div, sample_data(ps_rare))
write.csv(alpha_div, file.path(output_dir, 'alpha_diversity.csv'))

p_alpha <- ggplot(alpha_div, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    theme_minimal() + labs(title = 'Shannon Diversity') + theme(legend.position = 'none')
ggsave(file.path(output_dir, 'plots/alpha_diversity.pdf'), p_alpha, width = 6, height = 5)

bray_dist <- phyloseq::distance(ps_rare, method = 'bray')
meta_df <- data.frame(sample_data(ps_rare))
permanova <- adonis2(bray_dist ~ Group, data = meta_df, permutations = 999)
cat('  PERMANOVA R2:', round(permanova$R2[1], 3), 'p =', permanova$`Pr(>F)`[1], '\n')

pcoa <- ordinate(ps_rare, method = 'PCoA', distance = bray_dist)
p_beta <- plot_ordination(ps_rare, pcoa, color = 'Group') + stat_ellipse(level = 0.95) +
    theme_minimal() + labs(title = sprintf('PCoA (PERMANOVA R2=%.2f, p=%.3f)', permanova$R2[1], permanova$`Pr(>F)`[1]))
ggsave(file.path(output_dir, 'plots/beta_diversity_pcoa.pdf'), p_beta, width = 7, height = 6)

# === 8. DIFFERENTIAL ABUNDANCE ===
cat('\n8. Differential abundance testing...\n')
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) > 0.1 * nsamples(ps), TRUE)
otu <- as.data.frame(t(otu_table(ps_filt)))
groups <- as.character(sample_data(ps_filt)$Group)

aldex_out <- aldex(otu, groups, mc.samples = 128, test = 'welch', effect = TRUE, include.sample.summary = FALSE)
aldex_out$ASV <- rownames(aldex_out)
aldex_out$significant <- aldex_out$we.eBH < 0.05 & abs(aldex_out$effect) > 1
tax_df <- as.data.frame(tax_table(ps_filt))
aldex_out <- merge(aldex_out, tax_df, by.x = 'ASV', by.y = 'row.names', all.x = TRUE)
write.csv(aldex_out, file.path(output_dir, 'aldex2_results.csv'), row.names = FALSE)

cat('  Differential taxa:', sum(aldex_out$significant), '\n')

p_effect <- ggplot(aldex_out, aes(x = effect, y = -log10(we.eBH))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    scale_color_manual(values = c('grey60', 'firebrick')) +
    theme_minimal() + labs(title = 'ALDEx2 Effect Plot', x = 'Effect Size', y = '-log10(Adjusted P)')
ggsave(file.path(output_dir, 'plots/aldex2_effect_plot.pdf'), p_effect, width = 7, height = 6)

# === SUMMARY ===
cat('\n=== Pipeline Complete ===\n')
cat('Results saved to:', output_dir, '\n')
