# Differential accessibility analysis with DiffBind

library(DiffBind)
library(ggplot2)

# Configuration
output_dir <- 'atac_results/differential'
dir.create(output_dir, showWarnings = FALSE)

# Create sample sheet
samples <- data.frame(
    SampleID = c('control_1', 'control_2', 'treated_1', 'treated_2'),
    Condition = c('control', 'control', 'treated', 'treated'),
    Replicate = c(1, 2, 1, 2),
    bamReads = c('atac_results/aligned/control_1.shifted.bam',
                 'atac_results/aligned/control_2.shifted.bam',
                 'atac_results/aligned/treated_1.shifted.bam',
                 'atac_results/aligned/treated_2.shifted.bam'),
    Peaks = c('atac_results/peaks/control_1_peaks.narrowPeak',
              'atac_results/peaks/control_2_peaks.narrowPeak',
              'atac_results/peaks/treated_1_peaks.narrowPeak',
              'atac_results/peaks/treated_2_peaks.narrowPeak'),
    PeakCaller = rep('narrow', 4)
)

# Create DBA object
cat('Loading samples...\n')
dba_obj <- dba(sampleSheet = samples)
print(dba_obj)

# Count reads
cat('Counting reads in peaks...\n')
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps = TRUE)

# Normalize
cat('Normalizing...\n')
dba_obj <- dba.normalize(dba_obj)

# PCA and correlation
pdf(file.path(output_dir, 'sample_correlation.pdf'), width = 8, height = 6)
dba.plotHeatmap(dba_obj, correlations = TRUE)
dev.off()

pdf(file.path(output_dir, 'pca.pdf'), width = 6, height = 5)
dba.plotPCA(dba_obj, label = DBA_ID)
dev.off()

# Set up contrast
cat('Setting up contrast...\n')
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, minMembers = 2)

# Differential analysis
cat('Running differential analysis...\n')
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# Summary
cat('\nDifferential analysis summary:\n')
dba.show(dba_obj, bContrasts = TRUE)

# Get report
report <- dba.report(dba_obj, th = 0.05)
cat('\nDifferential peaks (FDR < 0.05):', length(report), '\n')

# Save report
report_df <- as.data.frame(report)
write.csv(report_df, file.path(output_dir, 'differential_peaks.csv'), row.names = FALSE)

# Separate up and down
up_peaks <- report_df[report_df$Fold > 0,]
down_peaks <- report_df[report_df$Fold < 0,]
cat('More accessible in treated:', nrow(up_peaks), '\n')
cat('Less accessible in treated:', nrow(down_peaks), '\n')

# Plots
pdf(file.path(output_dir, 'ma_plot.pdf'), width = 7, height = 6)
dba.plotMA(dba_obj)
dev.off()

pdf(file.path(output_dir, 'volcano_plot.pdf'), width = 7, height = 6)
dba.plotVolcano(dba_obj)
dev.off()

pdf(file.path(output_dir, 'venn.pdf'), width = 6, height = 5)
dba.plotVenn(dba_obj, contrast = 1, bDB = TRUE)
dev.off()

# Export BED files
export.bed <- function(gr, filename) {
    df <- data.frame(
        chr = seqnames(gr),
        start = start(gr) - 1,
        end = end(gr),
        name = paste0('peak_', seq_along(gr)),
        score = abs(gr$Fold) * 100
    )
    write.table(df, filename, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
}

export.bed(report[report$Fold > 0], file.path(output_dir, 'up_peaks.bed'))
export.bed(report[report$Fold < 0], file.path(output_dir, 'down_peaks.bed'))

cat('\nResults saved to:', output_dir, '\n')
