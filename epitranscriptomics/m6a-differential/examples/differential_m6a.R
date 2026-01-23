# Differential m6A analysis with exomePeak2

library(exomePeak2)
library(ggplot2)

# Sample setup
# 2 conditions (ctrl, treat) x 2 replicates x 2 fractions (IP, Input) = 8 BAMs
ip_bams <- c('ctrl_IP1.bam', 'ctrl_IP2.bam', 'treat_IP1.bam', 'treat_IP2.bam')
input_bams <- c('ctrl_Input1.bam', 'ctrl_Input2.bam', 'treat_Input1.bam', 'treat_Input2.bam')

# Experimental design
design <- data.frame(
    condition = factor(c('ctrl', 'ctrl', 'treat', 'treat'), levels = c('ctrl', 'treat'))
)

# Differential peak analysis
result <- exomePeak2(
    bam_ip = ip_bams,
    bam_input = input_bams,
    gff = 'genes.gtf',
    genome = 'hg38',
    paired_end = TRUE,
    experiment_design = design
)

# Extract differential results
diff_sites <- results(result, contrast = c('condition', 'treat', 'ctrl'))

# Apply thresholds
# padj < 0.05: FDR-corrected significance
# |log2FC| > 1: At least 2-fold change in methylation
sig_sites <- diff_sites[diff_sites$padj < 0.05 & abs(diff_sites$log2FoldChange) > 1, ]

print(paste('Total tested sites:', nrow(diff_sites)))
print(paste('Significant differential sites:', nrow(sig_sites)))
print(paste('Hypermethylated in treatment:', sum(sig_sites$log2FoldChange > 0)))
print(paste('Hypomethylated in treatment:', sum(sig_sites$log2FoldChange < 0)))

# Volcano plot
diff_df <- as.data.frame(diff_sites)
diff_df$significant <- diff_df$padj < 0.05 & abs(diff_df$log2FoldChange) > 1

pdf('differential_m6a_volcano.pdf', width = 8, height = 6)
ggplot(diff_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.5) +
    scale_color_manual(values = c('grey', 'red')) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'blue') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'blue') +
    labs(x = 'log2 Fold Change (treat/ctrl)', y = '-log10(adjusted p-value)',
         title = 'Differential m6A Methylation') +
    theme_minimal()
dev.off()

# Export results
write.csv(as.data.frame(sig_sites), 'differential_m6a_significant.csv', row.names = FALSE)
