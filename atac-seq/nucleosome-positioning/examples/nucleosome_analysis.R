#!/usr/bin/env Rscript
# Reference: ATACseqQC 1.26+, GenomicAlignments 1.38+, TxDb.Hsapiens.UCSC.hg38.knownGene 3.18+, BSgenome.Hsapiens.UCSC.hg38 1.4+ | Verify API if version differs
# Nucleosome positioning from ATAC-seq: V-plot diagnostic, fragment classes, TSS-anchored signal, +1 nucleosome aggregate.

suppressPackageStartupMessages({
    library(ATACseqQC); library(GenomicAlignments); library(GenomicRanges)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(BSgenome.Hsapiens.UCSC.hg38)
    library(ggplot2)
})

analyze_nucleosomes <- function(bam_file, output_prefix='nucleosome', upstream=1000, downstream=1000) {
    cat('=== Nucleosome positioning analysis ===\n')

    # Fragment-size distribution + diagnostic plot with Buenrostro 2013 windows
    pdf(sprintf('%s_fragsize.pdf', output_prefix), 8, 6)
    fragSizeDist(bam_file, output_prefix)
    dev.off()

    # MAPQ-filtered alignments
    gal <- readGAlignmentPairs(bam_file, param=ScanBamParam(mapqFilter=30))
    cat(sprintf('Paired alignments (MAPQ>=30): %d\n', length(gal)))

    # Fragment-class counts (Buenrostro 2013 / ENCODE convention)
    frags <- width(gal)
    nfr  <- gal[frags < 100]
    mono <- gal[frags >= 180 & frags <= 247]
    di   <- gal[frags >= 315 & frags <= 473]
    cat(sprintf('NFR (<100): %d (%.1f%%)\n', length(nfr), 100 * length(nfr) / length(gal)))
    cat(sprintf('Mono (180-247): %d (%.1f%%)\n', length(mono), 100 * length(mono) / length(gal)))
    cat(sprintf('Di (315-473): %d (%.1f%%)\n', length(di), 100 * length(di) / length(gal)))

    # Tn5 shift correction is required before splitting by fragment class
    cat('Applying Tn5 shift correction...\n')
    gal_shifted <- shiftGAlignmentsList(GAlignmentsList(gal))

    # Get TSS regions (protein-coding only avoids ncRNA noise)
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
    tss <- promoters(txs, upstream=upstream, downstream=downstream)

    # Split fragments into NFR / mono / di and compute aggregate signal at TSS
    cat('Computing TSS-aligned fragment-class signal...\n')
    objs <- splitGAlignmentsByCut(gal_shifted, txs=txs, genome=BSgenome.Hsapiens.UCSC.hg38)
    sigs <- featureAlignedSignal(cvglist=objs, feature.gr=tss,
                                 upstream=upstream, downstream=downstream)

    # V-plot at TSS (fragment size vs position; classic V/W indicates positioning is recoverable)
    cat('Generating V-plot at TSS...\n')
    pdf(sprintf('%s_vplot.pdf', output_prefix), 8, 6)
    vPlot(gal_shifted, tss, genome=BSgenome.Hsapiens.UCSC.hg38,
          upstream=upstream, downstream=downstream)
    dev.off()

    # Nucleosome positioning heatmap (NFR + mono signal aligned to TSS)
    pdf(sprintf('%s_heatmap.pdf', output_prefix), 8, 10)
    featureAlignedHeatmap(sigs, tss, upstream=upstream, downstream=downstream)
    dev.off()

    # Export NFR and mono BAMs for downstream tools
    rtracklayer::export(objs$NucleosomeFree %||% objs$NussomeFree,
                        sprintf('%s_nfr.bam', output_prefix))
    rtracklayer::export(objs$mononucleosome,
                        sprintf('%s_mono.bam', output_prefix))

    # Summary
    summary <- data.frame(
        sample = output_prefix,
        total_pairs = length(gal),
        nfr_count = length(nfr),
        mono_count = length(mono),
        di_count = length(di),
        nfr_to_mono_ratio = round(length(nfr) / length(mono), 2)
    )
    write.csv(summary, sprintf('%s_summary.csv', output_prefix), row.names=FALSE)
    cat('\nSummary:\n'); print(summary)
    invisible(summary)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) analyze_nucleosomes(args[1],
                                          if (length(args) > 1) args[2] else 'nucleosome')
