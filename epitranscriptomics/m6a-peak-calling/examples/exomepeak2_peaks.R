# Reference: exomePeak2 1.14+ (Bioconductor 3.18+), GenomicFeatures 1.54+, BSgenome.Hsapiens.UCSC.hg38 1.4+, rtracklayer 1.62+ | Verify with packageVersion('exomePeak2'); ?exomePeak2 if installed releases differ.
# m6A peak calling with exomePeak2 — transcript-aware, GC-bias-corrected, paired IP/input GLM.
# Output: BED12 + RDS + per-peak fold-change / FDR in save_dir/experiment_name/.

library(exomePeak2)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)

ip_bams <- c('aligned/IP_rep1_Aligned.sortedByCoord.out.bam',
             'aligned/IP_rep2_Aligned.sortedByCoord.out.bam',
             'aligned/IP_rep3_Aligned.sortedByCoord.out.bam')

input_bams <- c('aligned/Input_rep1_Aligned.sortedByCoord.out.bam',
                'aligned/Input_rep2_Aligned.sortedByCoord.out.bam',
                'aligned/Input_rep3_Aligned.sortedByCoord.out.bam')

txdb <- makeTxDbFromGFF('refs/annotation.gtf', format='gtf')

result <- exomePeak2(
    bam_ip          = ip_bams,
    bam_input       = input_bams,
    txdb            = txdb,
    genome          = BSgenome.Hsapiens.UCSC.hg38,
    paired_end      = TRUE,
    library_type    = 'unstranded',
    save_dir        = 'exomepeak2_output',
    experiment_name = 'm6a_run1'
)

peaks <- result
cat('exomePeak2 peaks called:', length(peaks), '\n')

# Flag TSS-proximal peaks as m6A-or-m6Am ambiguous (within 50 nt of transcript 5' end).
TSS_AMBIGUITY_NT <- 50

tx_5p <- promoters(txdb, upstream=0, downstream=TSS_AMBIGUITY_NT)

tss_proximal <- overlapsAny(peaks, tx_5p, ignore.strand=FALSE)
peaks$m6a_or_m6am_ambiguous <- tss_proximal

internal_peaks <- peaks[!tss_proximal]
ambiguous_peaks <- peaks[tss_proximal]

cat('Internal peaks (>=', TSS_AMBIGUITY_NT, 'nt from TSS):', length(internal_peaks), '\n')
cat('5UTR-proximal peaks (m6A/m6Am ambiguous):', length(ambiguous_peaks), '\n')

# Export internal-m6A peaks for downstream METTL3-biology analysis.
export(internal_peaks, 'exomepeak2_output/m6a_run1/peaks_internal.bed')
export(ambiguous_peaks, 'exomepeak2_output/m6a_run1/peaks_5utr_ambiguous.bed')
saveRDS(result, 'exomepeak2_output/m6a_run1/exomepeak2_result.rds')
