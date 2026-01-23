# Metagene visualization with Guitar

library(Guitar)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load m6A peaks
peaks <- import('m6a_peaks.bed')

# Convert to GRanges if needed
if (!is(peaks, 'GRanges')) {
    peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
}

# Standard Guitar metagene plot
# Shows m6A distribution across 5'UTR, CDS, 3'UTR
GuitarPlot(
    peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    saveToPDFprefix = 'm6a_metagene',
    enableCI = TRUE  # Confidence intervals
)

# Custom with specific features
# Focus on regions around stop codon where m6A typically enriched
guitar_coords <- GuitarCoords(
    TxDb.Hsapiens.UCSC.hg38.knownGene,
    noBins = 100
)

# Compare multiple peak sets
if (file.exists('m6a_ctrl.bed') && file.exists('m6a_treat.bed')) {
    ctrl_peaks <- import('m6a_ctrl.bed')
    treat_peaks <- import('m6a_treat.bed')

    GuitarPlot(
        list(Control = ctrl_peaks, Treatment = treat_peaks),
        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        saveToPDFprefix = 'm6a_comparison'
    )
}
