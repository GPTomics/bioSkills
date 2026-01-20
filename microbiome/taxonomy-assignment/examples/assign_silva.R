# Assign taxonomy to ASVs using SILVA database
library(dada2)

seqtab_nochim <- readRDS('seqtab_nochim.rds')
cat('ASVs to classify:', ncol(seqtab_nochim), '\n')

# Download SILVA from: https://zenodo.org/record/4587955
silva_train <- 'silva_nr99_v138.1_train_set.fa.gz'
silva_species <- 'silva_species_assignment_v138.1.fa.gz'

if (!file.exists(silva_train)) {
    stop('Download SILVA training set from https://zenodo.org/record/4587955')
}

# Assign taxonomy (genus level)
cat('Assigning taxonomy...\n')
taxa <- assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE, minBoot = 80)

# Add species-level assignments where possible
cat('Adding species assignments...\n')
if (file.exists(silva_species)) {
    taxa <- addSpecies(taxa, silva_species)
}

# Summarize results
cat('\nTaxonomy assignment summary:\n')
for (rank in colnames(taxa)) {
    assigned <- sum(!is.na(taxa[, rank]))
    cat(sprintf('  %s: %d/%d (%.1f%%)\n', rank, assigned, nrow(taxa), 100 * assigned / nrow(taxa)))
}

# Create readable output
taxa_df <- as.data.frame(taxa)
taxa_df$ASV <- rownames(taxa_df)
taxa_df$Sequence <- rownames(seqtab_nochim)
taxa_df <- taxa_df[, c('ASV', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Sequence')]
rownames(taxa_df) <- paste0('ASV', seq_len(nrow(taxa_df)))

write.csv(taxa_df, 'taxonomy_silva.csv', row.names = FALSE)
cat('\nSaved taxonomy_silva.csv\n')

# Save for phyloseq
saveRDS(taxa, 'taxa.rds')
