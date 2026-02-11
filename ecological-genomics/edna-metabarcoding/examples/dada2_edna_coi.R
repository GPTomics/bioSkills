library(dada2)
library(decontam)

# --- Input setup ---
input_dir <- 'cutadapt_trimmed'
fwd_files <- sort(list.files(input_dir, pattern = '_R1.fastq.gz', full.names = TRUE))
rev_files <- sort(list.files(input_dir, pattern = '_R2.fastq.gz', full.names = TRUE))
sample_names <- gsub('_R1.fastq.gz', '', basename(fwd_files))

filt_dir <- file.path(input_dir, 'filtered')
filt_fwd <- file.path(filt_dir, paste0(sample_names, '_F_filt.fastq.gz'))
filt_rev <- file.path(filt_dir, paste0(sample_names, '_R_filt.fastq.gz'))

# --- Quality inspection ---
plotQualityProfile(fwd_files[1:4])
plotQualityProfile(rev_files[1:4])

# --- Filter and trim ---
# maxEE=c(2,2): max expected errors per read; balances sensitivity and error removal
# truncLen: set based on quality profiles; COI mlCOIintF/jgHCO2198 produces ~313bp
# minLen=200: removes truncated fragments below useful COI barcode length
out <- filterAndTrim(fwd_files, filt_fwd, rev_files, filt_rev,
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 200, rm.phix = TRUE, multithread = TRUE)

exists_filter <- file.exists(filt_fwd) & file.exists(filt_rev)
filt_fwd <- filt_fwd[exists_filter]
filt_rev <- filt_rev[exists_filter]
sample_names <- sample_names[exists_filter]

# --- Learn error rates ---
err_fwd <- learnErrors(filt_fwd, multithread = TRUE)
err_rev <- learnErrors(filt_rev, multithread = TRUE)
plotErrors(err_fwd, nominalQ = TRUE)

# --- Denoise ---
dada_fwd <- dada(filt_fwd, err = err_fwd, multithread = TRUE)
dada_rev <- dada(filt_rev, err = err_rev, multithread = TRUE)

# --- Merge paired ends ---
# minOverlap=20: COI amplicon overlap region; increase if short overlap expected
merged <- mergePairs(dada_fwd, filt_fwd, dada_rev, filt_rev, minOverlap = 20)

# --- Build sequence table ---
seqtab <- makeSequenceTable(merged)
rownames(seqtab) <- sample_names
cat('ASV length distribution:\n')
table(nchar(getSequences(seqtab)))

# 280-340 bp: expected COI mlCOIintF/jgHCO2198 amplicon range
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% 280:340]

# --- Remove chimeras ---
# consensus method: per-sample chimera detection pooled across samples
# typical chimera rate 5-15%; >20% suggests library prep problems
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE)
chimera_pct <- round((1 - sum(seqtab_nochim) / sum(seqtab)) * 100, 1)
cat('Chimera rate:', chimera_pct, '%\n')

# --- Track reads through pipeline ---
track <- cbind(out[exists_filter, ], sapply(dada_fwd, function(x) sum(getUniques(x))),
               sapply(merged, function(x) sum(getUniques(x))),
               rowSums(seqtab_nochim))
colnames(track) <- c('input', 'filtered', 'denoised', 'merged', 'nonchim')
rownames(track) <- sample_names
track

# --- Contamination filtering with decontam ---
# Identify negative controls in the sample metadata
is_negative <- grepl('blank|negative|NTC', sample_names, ignore.case = TRUE)

if (sum(is_negative) > 0) {
    # prevalence method: ASVs more prevalent in negatives are contaminants
    # threshold=0.5: balanced false positive/negative rate
    contam <- isContaminant(seqtab_nochim, neg = is_negative,
                            method = 'prevalence', threshold = 0.5)
    cat('Contaminant ASVs identified:', sum(contam$contaminant), '\n')
    seqtab_clean <- seqtab_nochim[!is_negative, !contam$contaminant]
} else {
    cat('No negative controls found; skipping decontam.\n')
    seqtab_clean <- seqtab_nochim
}

# --- Taxonomy assignment ---
# MIDORI2 database formatted for DADA2 (download from midori2.info)
# minBoot=80: genus-level bootstrap confidence; use 50 for family-level
refdb <- 'MIDORI2_DADA2_COI.fasta.gz'
taxa <- assignTaxonomy(seqtab_clean, refdb, minBoot = 80, multithread = TRUE)

# --- Species-level assignment (exact matching) ---
# Only attempts species assignment for ASVs with genus already assigned
taxa_species <- addSpecies(taxa, 'MIDORI2_DADA2_COI_species.fasta.gz')

# --- Export results ---
asv_ids <- paste0('ASV_', seq_len(ncol(seqtab_clean)))
asv_seqs <- colnames(seqtab_clean)
colnames(seqtab_clean) <- asv_ids
rownames(taxa_species) <- asv_ids

write.csv(seqtab_clean, 'asv_abundance_table.csv')
write.csv(taxa_species, 'taxonomy_table.csv')
writeLines(paste0('>', asv_ids, '\n', asv_seqs), 'asv_sequences.fasta')

cat('Output files: asv_abundance_table.csv, taxonomy_table.csv, asv_sequences.fasta\n')
