# Reference: DADA2 1.30+, cutadapt 4.6+ | Verify API if version differs
# Multi-run DADA2 amplicon workflow: per-run error model -> mergeSequenceTables -> ONE chimera removal.
# Primers MUST be removed (cutadapt --discard-untrimmed) BEFORE this script runs; see remove_primers.sh.
# Leaving primers on corrupts the error model and fakes chimeras, so the order is non-negotiable.
library(dada2)

# truncLen is a DETECTION BUDGET, not just a quality cut: truncLen_F + truncLen_R must exceed
# amplicon_length + minOverlap (DADA2 default minOverlap = 12) or denoised pairs cannot merge.
# These values suit 16S V4 (~253 bp) on 2x250; for V3-V4 (~460 bp) keep more length and loosen maxEE_R.
truncLen <- c(240, 200)
maxEE <- c(2, 2)         # expected-errors filter (on the TRUNCATED read); beats a hard Q cutoff. Default 2.
truncQ <- 2             # truncate each read at the first base with Q <= 2

run_dirs <- c('run1', 'run2')   # one directory of primer-trimmed FASTQs per sequencing run

process_run <- function(run_path) {
    fnFs <- sort(list.files(run_path, pattern='_R1_001.fastq.gz', full.names=TRUE))
    fnRs <- sort(list.files(run_path, pattern='_R2_001.fastq.gz', full.names=TRUE))
    sample_names <- sapply(strsplit(basename(fnFs), '_'), `[`, 1)
    filtFs <- file.path(run_path, 'filtered', paste0(sample_names, '_F.fastq.gz'))
    filtRs <- file.path(run_path, 'filtered', paste0(sample_names, '_R.fastq.gz'))
    names(filtFs) <- sample_names
    names(filtRs) <- sample_names

    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen, maxEE=maxEE,
                         truncQ=truncQ, maxN=0, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

    errF <- learnErrors(filtFs, multithread=TRUE)   # fit THIS run only - error rates are run-specific
    errR <- learnErrors(filtRs, multithread=TRUE)

    # Inspect the fit: observed points must track the fitted line and fall with Q. On NovaSeq/NextSeq
    # binned quality the fit can go non-monotonic - enforce monotonicity before trusting the denoising.
    ggsave(file.path(run_path, 'error_fit_F.png'), plotErrors(errF, nominalQ=TRUE), width=8, height=6)

    dadaFs <- dada(filtFs, err=errF, multithread=TRUE)   # pool='pseudo' for rare-ASV sensitivity
    dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
    seqtab <- makeSequenceTable(mergers)

    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
    colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged')
    rownames(track) <- sample_names
    list(seqtab=seqtab, track=track)
}

per_run <- lapply(run_dirs, process_run)

# Combine AFTER per-run inference: mergeSequenceTables joins by the exact ASV sequence string,
# which is only possible because ASVs are exact sequences (the operational payoff of "ASVs replace OTUs").
seqtab_all <- mergeSequenceTables(tables=lapply(per_run, `[[`, 'seqtab'))

# ONE chimera removal on the combined table. Chimeras are many ASVs but few READS; a large read loss
# here is a leftover-primer smell (degenerate bases look chimeric), not a real chimera storm.
seqtab_nochim <- removeBimeraDenovo(seqtab_all, method='consensus', multithread=TRUE, verbose=TRUE)
cat('reads retained after chimera removal:', round(100 * sum(seqtab_nochim) / sum(seqtab_all), 1), '%\n')

# ASV length distribution - off-target lengths flag mis-merges or non-target amplification.
table(nchar(getSequences(seqtab_nochim)))

track_all <- do.call(rbind, lapply(per_run, `[[`, 'track'))
write.csv(track_all, 'read_tracking.csv')
saveRDS(seqtab_nochim, 'seqtab_nochim.rds')   # carry 'run' as a batch covariate into downstream stats
