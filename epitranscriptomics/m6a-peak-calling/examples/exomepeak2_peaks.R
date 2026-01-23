# m6A peak calling with exomePeak2

library(exomePeak2)

ip_bams <- c('IP_rep1.bam', 'IP_rep2.bam')
input_bams <- c('Input_rep1.bam', 'Input_rep2.bam')
gtf_file <- 'genes.gtf'

# Peak calling
# Uses negative binomial model to identify enriched regions
result <- exomePeak2(
    bam_ip = ip_bams,
    bam_input = input_bams,
    gff = gtf_file,
    genome = 'hg38',
    paired_end = TRUE,
    # p < 0.05: Standard threshold for peak calling
    p_cutoff = 0.05,
    # log2FC > 1: Require at least 2-fold enrichment
    log2FC_cutoff = 1
)

# View peak summary
print(result)

# Get peaks as GRanges
peaks <- granges(result)
print(paste('Found', length(peaks), 'm6A peaks'))

# Export to BED
exportResults(result, format = 'BED', file = 'm6a_peaks.bed')

# Export full results with statistics
exportResults(result, format = 'CSV', file = 'm6a_peaks.csv')

# Peak annotation summary
anno <- as.data.frame(result)
table(anno$feature)  # Distribution across UTRs, CDS, etc
