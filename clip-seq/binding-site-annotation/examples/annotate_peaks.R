# Annotate CLIP-seq peaks

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('peaks.bed')

anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-1000, 1000)
)

# Summary
print(anno)

# Pie chart
plotAnnoPie(anno)

# Save annotations
write.csv(as.data.frame(anno), 'peaks_annotated.csv')
