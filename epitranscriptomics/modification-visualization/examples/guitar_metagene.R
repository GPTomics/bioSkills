# Reference: Guitar 2.18+ (Bioconductor), TxDb.Hsapiens.UCSC.hg38.knownGene 3.18+, rtracklayer 1.62+ | Verify with packageVersion('Guitar'); ?GuitarPlot if installed releases differ.
# THE canonical m6A metagene plot — Guitar transcript-feature-scaled, expecting stop-codon-proximal enrichment.
# Stop-codon enrichment is the biological QC anchor (Dominissini 2012 *Nature* 485:201; Meyer 2012 *Cell* 149:1635).

library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Single-condition metagene: confirm canonical stop-codon-proximal enrichment.
# Older Guitar versions: argument is `txdb=`; newer: `txTxdb=`. Verify with ?GuitarPlot.
GuitarPlot(
    txTxdb            = txdb,
    stBedFiles        = list('exomepeak2_output/m6a_run1/peaks.bed'),
    miscOutFilePrefix = 'figures/m6a_metagene',
    enableCI          = FALSE,
    saveToPDFprefix   = 'figures/m6a_metagene'
)

# Multi-condition comparison (WT vs KO) as overlaid metagenes.
if (file.exists('exomepeak2_output/wt/peaks.bed') &&
    file.exists('exomepeak2_output/ko/peaks.bed')) {

    GuitarPlot(
        txTxdb            = txdb,
        stBedFiles        = list(
            WT = 'exomepeak2_output/wt/peaks.bed',
            KO = 'exomepeak2_output/ko/peaks.bed'
        ),
        miscOutFilePrefix = 'figures/m6a_wt_vs_ko_metagene',
        enableCI          = TRUE,
        saveToPDFprefix   = 'figures/m6a_wt_vs_ko_metagene'
    )
}

# Visual sanity check: if the resulting metagene PDF does NOT show density rising toward and
# peaking near the stop codon in the 3'UTR, do NOT proceed to downstream visualisation.
# Re-inspect merip-preprocessing IP enrichment QC; investigate antibody / protocol failure.

cat('Wrote figures/m6a_metagene*.pdf.\n')
cat('Expected pattern: density rising toward and peaking near stop codon in 3UTR-proximal CDS.\n')
cat('If pattern absent: IP failure suspect; re-inspect plotFingerprint in merip-preprocessing.\n')
