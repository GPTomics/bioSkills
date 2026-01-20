#!/usr/bin/env Rscript
# Nucleosome positioning analysis from ATAC-seq

library(ATACseqQC)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

analyze_nucleosomes <- function(bam_file, output_prefix = 'nucleosome') {
    cat('Loading BAM file...\n')

    # Fragment size distribution
    cat('Calculating fragment size distribution...\n')
    frag_sizes <- fragSizeDist(bam_file, output_prefix)

    # Plot fragment sizes
    pdf(paste0(output_prefix, '_fragsize.pdf'), width = 8, height = 6)
    plot(frag_sizes, type = 'l', xlab = 'Fragment Size (bp)', ylab = 'Proportion',
         main = 'ATAC-seq Fragment Size Distribution')
    abline(v = c(100, 180, 247, 315, 473), lty = 2, col = 'grey')
    text(x = c(50, 150, 210, 280, 400), y = max(frag_sizes$proportion) * 0.9,
         labels = c('NFR', 'Mono', '', 'Di', ''), cex = 0.8)
    dev.off()

    # Read BAM
    gal <- readGAlignmentPairs(bam_file, param = ScanBamParam(mapqFilter = 30))

    # Split by nucleosome occupancy
    cat('Splitting reads by fragment size...\n')
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

    # NFR: < 100bp
    # Mono: 180-247bp
    # Di: 315-473bp
    nfr <- gal[width(gal) < 100]
    mono <- gal[width(gal) >= 180 & width(gal) <= 247]
    di <- gal[width(gal) >= 315 & width(gal) <= 473]

    cat(sprintf('NFR reads: %d\n', length(nfr)))
    cat(sprintf('Mono-nucleosome reads: %d\n', length(mono)))
    cat(sprintf('Di-nucleosome reads: %d\n', length(di)))

    # TSS enrichment score
    cat('Calculating TSS enrichment...\n')
    tsse <- TSSEscore(gal, txs)
    cat(sprintf('TSS enrichment score: %.2f\n', tsse$TSSEscore))

    # Save results
    results <- list(
        fragment_sizes = frag_sizes,
        nfr_count = length(nfr),
        mono_count = length(mono),
        di_count = length(di),
        tss_enrichment = tsse$TSSEscore
    )

    saveRDS(results, paste0(output_prefix, '_results.rds'))
    cat('Saved:', paste0(output_prefix, '_results.rds\n'))

    return(results)
}

# Run if executed directly
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    bam_file <- args[1]
    output_prefix <- if (length(args) > 1) args[2] else 'nucleosome'
    analyze_nucleosomes(bam_file, output_prefix)
} else {
    cat('Usage: Rscript nucleosome_analysis.R input.bam [output_prefix]\n')
}
