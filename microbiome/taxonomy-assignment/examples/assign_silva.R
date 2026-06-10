# Reference: DADA2 1.30+, DECIPHER 2.30+ | Verify API if version differs
# Assign taxonomy to 16S ASVs with DADA2 naive Bayes (genus) + addSpecies (exact-match species).
# Honest message: a short 16S read licenses genus at best; species comes only from an exact match.
library(dada2)

seqtab_nochim <- readRDS('seqtab_nochim.rds')
cat('ASVs to classify:', ncol(seqtab_nochim), '\n')

# The reference must match the marker AND ideally the primer region. A FULL-LENGTH SILVA training
# set applied to V4 (~250 bp) reads mismatches k-mer composition and degrades calls (Werner 2012;
# Bokulich 2018). For V4 data prefer a region-extracted reference (see assign_qiime2_region.sh).
silva_train <- 'silva_nr99_v138.1_train_set.fa.gz'        # DADA2-formatted, rank-labelled headers
silva_species <- 'silva_species_assignment_v138.1.fa.gz'  # species-level reference for exact match

minBoot <- 50  # DADA2 default and the RDP recommendation for reads <=250 nt; tutorials often use
               # 80 (a stricter CHOICE). Raising it truncates to shallower-but-reliable ranks;
               # ranks below the floor are returned as NA, not guessed.

taxa <- assignTaxonomy(seqtab_nochim, silva_train, minBoot = minBoot, tryRC = TRUE, multithread = TRUE)

# addSpecies assigns species ONLY by exact (100%) match against the species reference - it does
# not infer species from a noisy read. Most ASVs stay NA at species; that is the honest 16S result.
if (file.exists(silva_species)) {
    taxa <- addSpecies(taxa, silva_species)
}

cat('\nAssigned fraction per rank (NA = unassigned at threshold, kept honestly):\n')
for (rank in colnames(taxa)) {
    assigned <- sum(!is.na(taxa[, rank]))
    cat(sprintf('  %-8s %d/%d (%.1f%%)\n', rank, assigned, nrow(taxa), 100 * assigned / nrow(taxa)))
}

# State the conditioning choices alongside the labels (classifier + database release + region).
attr(taxa, 'classifier') <- 'DADA2 assignTaxonomy (RDP naive Bayes)'
attr(taxa, 'reference') <- 'SILVA 138.1'
attr(taxa, 'region') <- 'verify the reference region matches the amplicon primers'

saveRDS(taxa, 'taxa.rds')
cat('\nSaved taxa.rds (genus via naive Bayes, species via exact match only)\n')
