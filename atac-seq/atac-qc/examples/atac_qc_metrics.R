#!/usr/bin/env Rscript
# ATAC-seq QC metrics

library(ATACseqQC)
library(GenomicAlignments)

calculate_atac_qc <- function(bam_file, peaks_file = NULL, output_prefix = 'atac_qc') {
    cat('=== ATAC-seq QC ===\n')

    # Fragment size distribution
    cat('Fragment size distribution...\n')
    frag_sizes <- fragSizeDist(bam_file, output_prefix)

    # TSS enrichment (requires TxDb)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    gal <- readGAlignmentPairs(bam_file, param = ScanBamParam(mapqFilter = 30))
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
    tsse <- TSSEscore(gal, txs)
    cat(sprintf('TSS enrichment score: %.2f\n', tsse$TSSEscore))

    # FRiP (Fraction of Reads in Peaks)
    if (!is.null(peaks_file) && file.exists(peaks_file)) {
        library(GenomicRanges)
        peaks <- import(peaks_file)
        reads_in_peaks <- countOverlaps(gal, peaks)
        frip <- sum(reads_in_peaks > 0) / length(gal)
        cat(sprintf('FRiP: %.3f\n', frip))
    }

    # NFR ratio
    nfr <- sum(width(gal) < 100)
    mono <- sum(width(gal) >= 180 & width(gal) <= 247)
    nfr_ratio <- nfr / (nfr + mono)
    cat(sprintf('NFR ratio: %.3f\n', nfr_ratio))

    # Save metrics
    metrics <- data.frame(
        sample = output_prefix,
        tss_enrichment = tsse$TSSEscore,
        nfr_ratio = nfr_ratio,
        total_fragments = length(gal)
    )
    write.csv(metrics, paste0(output_prefix, '_metrics.csv'), row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    calculate_atac_qc(args[1], args[2], args[3])
}
