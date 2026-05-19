# Reference: ChIPseeker 1.40+, GenomicFeatures 1.56+, TxDb.Hsapiens.UCSC.hg38.knownGene 3.20+ | Verify API if version differs
# Annotate CLIP-seq peaks with CLIP-appropriate priority hierarchy.
# Default ChIPseeker tssRegion c(-3000, 3000) over-extends for CLIP - 6 kb "promoter" window catches
# deep 5' UTR peaks and labels them spuriously. Tighten to c(-100, 100).
# Use level='transcript' to preserve isoform context (splicing factors need this).

library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('peaks.stringent.bed')

# CLIP-appropriate annotation
# Priority hierarchy explicit; promoter window tight
anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    level = 'transcript',
    tssRegion = c(-100, 100),
    genomicAnnotationPriority = c(
        'Promoter', '5UTR', '3UTR', 'Exon',
        'Intron', 'Downstream', 'Intergenic'
    )
)

cat('\nGlobal region distribution\n')
print(anno)

# Region pie chart
png('clip_region_pie.png', width=800, height=800, res=150)
plotAnnoPie(anno)
dev.off()

# Distance to TSS (= distance to mRNA 5' end for CLIP)
png('clip_dist_to_tss.png', width=1000, height=600, res=150)
plotDistToTSS(anno, title='CLIP peak distance to mRNA TSS')
dev.off()

# Per-peak annotation table
anno_df <- as.data.frame(anno)
write.table(anno_df, 'clip_peaks_annotated.tsv', sep='\t', quote=FALSE, row.names=FALSE)

cat('\nFraction by region:\n')
region_frac <- table(anno_df$annotation) / nrow(anno_df)
print(round(region_frac, 3))

# Diagnostic: RBP-class expected distributions
# HuR / ELAVL1 -> > 50% 3' UTR
# PTBP1 -> > 50% intron, with peak < 500 nt of splice junctions
# EIF3 / RPS proteins -> 5' UTR / CDS / snoRNA
# FASTKD2 -> chrM
# MATR3 -> intron with > 15% repeat (Alu/LINE)
cat('\nIf >> 30% in Promoter category, tighten tssRegion.\n')
cat('If >> 50% Intergenic, verify TxDb covers all annotated transcripts (GENCODE v38+).\n')
cat('If chrM peaks missing, TxDb may lack mt-chromosome; add custom mt-mRNA BED.\n')

# Repeat-element axis (separate from region)
# Run separately:
#   bedtools intersect -s -wa -u -a peaks.stringent.bed -b RepeatMasker.bed > peaks.repeat.bed
#   bedtools intersect -s -wa -v -a peaks.stringent.bed -b RepeatMasker.bed > peaks.norepeat.bed
# Cross-tabulate against region axis for full picture.

cat('\nDone. Region pie: clip_region_pie.png; per-peak table: clip_peaks_annotated.tsv\n')
cat('For splicing factors, additionally run RBP-Maps (Yeo lab) with cassette-exon BED.\n')
