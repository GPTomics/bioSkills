# DADA2 amplicon processing workflow for paired-end 16S data
library(dada2)

path <- 'raw_reads'
fnFs <- sort(list.files(path, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = '_R2_001.fastq.gz', full.names = TRUE))
sample_names <- sapply(strsplit(basename(fnFs), '_'), `[`, 1)

cat('Processing', length(fnFs), 'samples\n')

# Quality profiles (check first 2 samples)
pdf('quality_profiles.pdf')
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
dev.off()

# Setup filtered files
dir.create('filtered', showWarnings = FALSE)
filtFs <- file.path('filtered', paste0(sample_names, '_F_filt.fastq.gz'))
filtRs <- file.path('filtered', paste0(sample_names, '_R_filt.fastq.gz'))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Filter and trim (adjust truncLen based on quality profiles)
# truncLen: Set where Q-score drops below ~25-30; forward usually longer than reverse
# maxEE: Max expected errors; 2 is DADA2 default, use 1 for stricter filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Denoise
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
cat('ASVs before chimera removal:', ncol(seqtab), '\n')

# Remove chimeras; 'consensus' method is more conservative than 'pooled', reducing false positives
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE, verbose = TRUE)
cat('ASVs after chimera removal:', ncol(seqtab_nochim), '\n')
cat('Reads retained:', round(100 * sum(seqtab_nochim) / sum(seqtab), 1), '%\n')

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab_nochim))
colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim')
rownames(track) <- sample_names
write.csv(track, 'read_tracking.csv')

# Save sequence table
saveRDS(seqtab_nochim, 'seqtab_nochim.rds')
cat('Saved seqtab_nochim.rds\n')
