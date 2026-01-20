# ChIP-seq peak annotation with ChIPseeker

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)

# Configuration
peak_file <- 'chipseq_results/peaks/experiment_peaks.narrowPeak'
output_dir <- 'chipseq_results/annotation'
dir.create(output_dir, showWarnings = FALSE)

# Load annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peaks
peaks <- readPeakFile(peak_file)
cat('Total peaks:', length(peaks), '\n')

# Annotate peaks
peak_anno <- annotatePeak(peaks, TxDb = txdb, annoDb = 'org.Hs.eg.db',
                          tssRegion = c(-3000, 3000), verbose = FALSE)

# Visualizations
pdf(file.path(output_dir, 'annotation_plots.pdf'), width = 10, height = 8)

plotAnnoPie(peak_anno)
plotAnnoBar(peak_anno)
plotDistToTSS(peak_anno, title = 'Distribution of peaks relative to TSS')

# Coverage plot around TSS
peakHeatmap(peaks, TxDb = txdb, upstream = 3000, downstream = 3000)

dev.off()

# Export annotated peaks
anno_df <- as.data.frame(peak_anno)
write.csv(anno_df, file.path(output_dir, 'annotated_peaks.csv'), row.names = FALSE)

# Summary statistics
cat('\n=== Annotation Summary ===\n')
anno_table <- table(anno_df$annotation)
for (i in seq_along(anno_table)) {
    cat(sprintf('%s: %d (%.1f%%)\n', names(anno_table)[i], anno_table[i],
                100 * anno_table[i] / sum(anno_table)))
}

# Extract genes with promoter peaks
promoter_peaks <- anno_df[grepl('Promoter', anno_df$annotation),]
promoter_genes <- unique(promoter_peaks$SYMBOL)
promoter_genes <- promoter_genes[!is.na(promoter_genes)]

cat('\nGenes with peaks in promoter:', length(promoter_genes), '\n')
write.table(promoter_genes, file.path(output_dir, 'promoter_genes.txt'),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Extract all target genes
all_genes <- unique(anno_df$SYMBOL)
all_genes <- all_genes[!is.na(all_genes)]
write.table(all_genes, file.path(output_dir, 'all_target_genes.txt'),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat('All target genes:', length(all_genes), '\n')
cat('\nResults saved to:', output_dir, '\n')
