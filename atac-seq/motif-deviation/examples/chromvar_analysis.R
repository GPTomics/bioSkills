#!/usr/bin/env Rscript
# Reference: chromVAR 1.24+, motifmatchr 1.24+, JASPAR2024 0.99+, limma 3.58+ | Verify API if version differs
# chromVAR bulk workflow: peak counts -> matched-background z-scores -> variability + limma differential.

suppressPackageStartupMessages({
    library(chromVAR); library(motifmatchr); library(BSgenome.Hsapiens.UCSC.hg38)
    library(JASPAR2024); library(TFBSTools); library(SummarizedExperiment)
    library(GenomicRanges); library(limma); library(pheatmap); library(RSQLite)
})

# Load consensus peaks + counts
peaks <- read.table('peaks.bed', col.names=c('chr', 'start', 'end'))
peak_ranges <- GRanges(seqnames=peaks$chr, ranges=IRanges(peaks$start + 1, peaks$end))

counts_matrix <- as.matrix(read.delim('counts.tsv', row.names=1))

se <- SummarizedExperiment(assays=list(counts=counts_matrix), rowRanges=peak_ranges)
colData(se)$condition <- factor(c('control','control','control','treated','treated','treated'),
                                levels=c('control','treated'))

# Add GC content for matched-bias background sampling
se <- addGCBias(se, genome=BSgenome.Hsapiens.UCSC.hg38)

# chromVAR thresholds: 1500 reads/sample, FRiP >= 0.15; 10 reads/peak in >= 10% samples
se <- filterSamples(se, min_depth=1500, min_in_peaks=0.15, shiny=FALSE)
se <- filterPeaks(se, non_overlapping=TRUE, min_count=10, n_samples_frac=0.1)
cat(sprintf('After filtering: %d peaks, %d samples\n', nrow(se), ncol(se)))
stopifnot(nrow(se) >= 5000)             # chromVAR requires >= 5000 peaks for stable inference

# JASPAR 2024 vertebrate CORE motifs (~1900)
# TFBSTools issue #39: getMatrixSet does not dispatch on JASPAR2024 directly.
# Use the SQLite handle as the workaround.
jaspar2024 <- JASPAR2024::JASPAR2024()
sq <- dbConnect(SQLite(), db(jaspar2024))
pfm <- getMatrixSet(sq, opts=list(collection='CORE', tax_group='vertebrates'))
motif_ix <- matchMotifs(pfm, se, genome=BSgenome.Hsapiens.UCSC.hg38, p.cutoff=5e-5)

# Background peaks: 50 iterations is the default; do not reduce below 30
bg <- getBackgroundPeaks(object=se, niterations=50)
dev <- computeDeviations(object=se, annotations=motif_ix, background_peaks=bg)
zscores <- deviations(dev)
variability <- computeVariability(dev)
variability <- variability[order(-variability$variability), ]
cat('\nTop 10 variable motifs:\n'); print(head(variability, 10))

# Differential motif activity via limma on z-scores
# Use adj.P.Val (limma's BH FDR), not 'FDR' which limma does not return
groups <- colData(se)$condition
design <- model.matrix(~groups)
fit <- eBayes(lmFit(zscores, design))
diff_motifs <- topTable(fit, coef=2, number=Inf, sort.by='P')
sig <- diff_motifs[diff_motifs$adj.P.Val < 0.05 & abs(diff_motifs$logFC) >= 0.5, ]
cat(sprintf('\nDifferential motifs (FDR < 0.05, |logFC| >= 0.5): %d\n', nrow(sig)))
print(head(sig, 10))

# Plot variability (signal vs noise)
pdf('chromvar_variability.pdf', 8, 6)
plotVariability(variability, use_plotly=FALSE)
dev.off()

# Heatmap top 50 variable motifs, scaled across samples
top50 <- head(rownames(variability), 50)
sample_info <- data.frame(Condition=colData(se)$condition, row.names=colnames(zscores))
pdf('chromvar_heatmap.pdf', 10, 12)
pheatmap(zscores[top50, ], annotation_col=sample_info, scale='row',
         clustering_method='ward.D2', fontsize_row=8)
dev.off()

# Save outputs
write.csv(as.data.frame(zscores), 'chromvar_deviations.csv')
write.csv(variability, 'chromvar_variability.csv')
write.csv(diff_motifs, 'chromvar_differential.csv')
cat('\nResults saved.\n')
