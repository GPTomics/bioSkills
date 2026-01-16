library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('sample_peaks.narrowPeak')

promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

tagMatrix <- getTagMatrix(peaks, windows = promoter)

tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = 'red')

plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = 'Distance from TSS (bp)',
            ylab = 'Peak Count Frequency', conf = 0.95)

peak_files <- list(
    H3K4me3 = 'H3K4me3_peaks.narrowPeak',
    H3K27ac = 'H3K27ac_peaks.narrowPeak'
)

tagMatrixList <- lapply(peak_files, function(f) {
    peaks <- readPeakFile(f)
    getTagMatrix(peaks, windows = promoter)
})

plotAvgProf(tagMatrixList, xlim = c(-3000, 3000), xlab = 'Distance from TSS (bp)',
            ylab = 'Peak Count Frequency', facet = 'row')
