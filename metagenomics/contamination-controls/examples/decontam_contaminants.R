# Reference: decontam 1.22+, phyloseq 1.46+ | Verify API if version differs
# Flag reagent/kitome contaminants from a shotgun classifier table using BOTH decontam signals:
# frequency (contaminants scale inversely with input DNA) and prevalence (enriched in blanks).
# decontam runs on the FEATURE TABLE (Bracken/MetaPhlAn output), not raw reads.
library(decontam)

# seqtab: samples x taxa matrix of abundances; dna_conc: per-sample DNA concentration (Qubit/PicoGreen);
# is_blank: logical, TRUE for extraction/no-template blanks; batch_id: lot/run (the kitome differs by lot).
# (Load these from your abundance table and sample metadata.)

# Combined signal, standard threshold.
contam <- isContaminant(seqtab, conc = dna_conc, neg = is_blank,
                        method = 'combined', threshold = 0.1, batch = batch_id)

# Low-biomass studies: prevalence method at the more aggressive 0.5 - calls anything more prevalent in
# blanks than samples a contaminant. ALWAYS inspect the calls before removing.
contam_lowbio <- isContaminant(seqtab, neg = is_blank,
                               method = 'prevalence', threshold = 0.5, batch = batch_id)

kitome_genera <- c('Bradyrhizobium', 'Ralstonia', 'Burkholderia', 'Pseudomonas',
                   'Acinetobacter', 'Sphingomonas', 'Methylobacterium', 'Stenotrophomonas')
flagged <- rownames(contam)[contam$contaminant]
cat('Flagged contaminants:', length(flagged), '\n')
cat('Of which canonical kitome genera:',
    sum(sapply(kitome_genera, function(g) any(grepl(g, flagged)))), '\n')

# Remove only after inspection - over-aggressive removal deletes real low-abundance taxa.
seqtab_clean <- seqtab[, !contam$contaminant]
