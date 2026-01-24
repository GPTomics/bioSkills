library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(SummarizedExperiment)
library(pheatmap)
library(ggplot2)

# For test data, use chromVAR example dataset:
# data(example_counts, package = 'chromVAR')

# Load peak counts from files
peaks <- read.table('peaks.bed', col.names = c('chr', 'start', 'end'))
peak_ranges <- GRanges(seqnames = peaks$chr, ranges = IRanges(peaks$start, peaks$end))

counts <- read.table('counts.txt', header = TRUE, row.names = 1)
counts_matrix <- as.matrix(counts)

fragment_counts <- SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowRanges = peak_ranges
)

# Add sample metadata
colData(fragment_counts)$condition <- c('Control', 'Control', 'Treated', 'Treated')

# Add GC bias correction
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)

# min_depth=1500: Minimum reads per sample. Lower for shallow sequencing.
# min_in_peaks=0.15: FRiP threshold. Typical range 0.1-0.3 for ATAC-seq.
fragment_counts <- filterSamples(fragment_counts, min_depth = 1500, min_in_peaks = 0.15)

# min_count=10: Peak read threshold across samples.
fragment_counts <- filterPeaks(fragment_counts, non_overlapping = TRUE, min_count = 10)

cat('Filtered peaks:', nrow(fragment_counts), '\n')
cat('Samples:', ncol(fragment_counts), '\n')

# Get JASPAR motifs
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = 'CORE', tax_group = 'vertebrates'))
cat('Motifs:', length(pfm), '\n')

# Match motifs to peaks
# p.cutoff=5e-5: Match threshold. Use 1e-4 for more matches, 1e-5 for stringent.
motif_ix <- matchMotifs(pfm, fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 5e-5)

# Compute deviations
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix)
deviation_scores <- deviations(dev)

# Compute variability
variability <- computeVariability(dev)
variability <- variability[order(-variability$variability), ]

cat('\nTop 10 variable motifs:\n')
print(head(variability, 10))

# Variability plot
pdf('chromvar_variability.pdf', width = 8, height = 6)
plotVariability(variability, use_plotly = FALSE)
dev.off()

# Heatmap of top variable motifs
# n_top=50: Adjust based on how many TFs you want to display.
n_top <- 50
top_motifs <- head(rownames(variability), n_top)
top_dev <- deviation_scores[top_motifs, ]

sample_info <- data.frame(
    Condition = colData(fragment_counts)$condition,
    row.names = colnames(top_dev)
)

pdf('chromvar_heatmap.pdf', width = 10, height = 12)
pheatmap(top_dev, annotation_col = sample_info, scale = 'row',
         clustering_method = 'ward.D2', fontsize_row = 8)
dev.off()

# Save results
write.csv(as.data.frame(deviation_scores), 'chromvar_deviations.csv')
write.csv(variability, 'chromvar_variability.csv')

cat('\nResults saved to chromvar_deviations.csv and chromvar_variability.csv\n')
