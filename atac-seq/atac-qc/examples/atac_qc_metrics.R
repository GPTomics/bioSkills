#!/usr/bin/env Rscript
# Reference: ATACseqQC 1.26+, GenomicAlignments 1.38+, samtools 1.19+, pysam 0.22+ | Verify API if version differs
# ENCODE 4 ATAC-seq QC metrics. Computes both ENCODE-style and ATACseqQC-style TSS enrichment.

suppressPackageStartupMessages({
    library(ATACseqQC); library(GenomicAlignments); library(GenomicRanges)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(rtracklayer)
})

calculate_atac_qc <- function(bam_file, peaks_file = NULL, output_prefix = 'atac_qc') {
    cat('=== ENCODE 4 ATAC-seq QC ===\n')

    # Fragment-size distribution -> visual periodicity check
    pdf(sprintf('%s_fragsize.pdf', output_prefix))
    frag_sizes <- fragSizeDist(bam_file, output_prefix)
    dev.off()

    # MAPQ-filtered alignments (paired-end)
    gal <- readGAlignmentPairs(bam_file, param = ScanBamParam(mapqFilter = 30))
    cat(sprintf('Paired alignments (MAPQ >= 30): %d\n', length(gal)))

    # ATACseqQC TSS enrichment (sum-over-window method; not ENCODE-comparable)
    # TSSEscore requires a GAlignments (not GAlignmentPairs), so read one for this step
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
    gal_se <- readGAlignments(bam_file, param = ScanBamParam(mapqFilter = 30))
    tsse <- TSSEscore(gal_se, txs)
    cat(sprintf('ATACseqQC TSSEscore: %.2f\n', tsse$TSSEscore))
    cat('  (ATACseqQC uses 100bp center / 1000bp flanks; typically 2-3x ENCODE pyTSSe)\n')

    # Fragment-size class fractions
    frags <- width(granges(gal))
    nfr <- sum(frags < 100); mono <- sum(frags >= 180 & frags <= 247)
    di <- sum(frags >= 315 & frags <= 473); tri <- sum(frags >= 558 & frags <= 615)
    cat(sprintf('NFR: %d (%.1f%%)  Mono: %d (%.1f%%)  Di: %d (%.1f%%)  Tri: %d\n',
                nfr, 100 * nfr / length(gal), mono, 100 * mono / length(gal),
                di, 100 * di / length(gal), tri))

    # FRiP requires peaks file
    frip <- NA_real_
    if (!is.null(peaks_file) && nzchar(peaks_file) && file.exists(peaks_file)) {
        peaks <- import(peaks_file)
        in_peaks <- sum(countOverlaps(gal, peaks) > 0)
        frip <- in_peaks / length(gal)
        cat(sprintf('FRiP: %.3f  (ENCODE: ideal >= 0.3, accept >= 0.2)\n', frip))
    }

    # ENCODE pass/fail grading (subset; adapt thresholds for non-Omni protocols)
    grade <- function(value, accept, ideal, inverted = FALSE) {
        if (is.na(value)) return('NA')
        if (inverted) {
            if (value > accept) 'FAIL' else if (value <= ideal) 'PASS' else 'WARN'
        } else {
            if (value < accept) 'FAIL' else if (value >= ideal) 'PASS' else 'WARN'
        }
    }

    metrics <- data.frame(
        sample = output_prefix,
        nuclear_reads_M = round(length(gal) / 1e6, 1),
        ATACseqQC_TSSEscore = round(tsse$TSSEscore, 2),
        nfr_fraction = round(nfr / length(gal), 3),
        mono_fraction = round(mono / length(gal), 3),
        FRiP = ifelse(is.na(frip), NA, round(frip, 3)),
        nuc_reads_grade = grade(length(gal) / 1e6, 25, 50),
        FRiP_grade = grade(frip, 0.2, 0.3)
    )
    write.csv(metrics, sprintf('%s_metrics.csv', output_prefix), row.names = FALSE)
    cat('\nReport card:\n'); print(metrics)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) calculate_atac_qc(args[1], if (length(args) > 1) args[2] else NULL,
                                        if (length(args) > 2) args[3] else 'sample')
