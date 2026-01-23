# Differential expression analysis of miRNAs using DESeq2

library(DESeq2)

# Load count matrix (from miRge3 or miRDeep2)
# Rows = miRNAs, Columns = samples
counts <- read.csv('miR.Counts.csv', row.names = 1)

# Create sample metadata
# Adjust according to your experimental design
coldata <- data.frame(
    sample = colnames(counts),
    condition = factor(c('control', 'control', 'control',
                        'treated', 'treated', 'treated')),
    row.names = colnames(counts)
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = ~ condition
)

# Filter low-expressed miRNAs
# Threshold: 10 reads total across samples
# Lower than typical mRNA threshold because miRNAs are smaller/fewer
MIN_COUNT <- 10
keep <- rowSums(counts(dds)) >= MIN_COUNT
dds <- dds[keep, ]
cat('Kept', sum(keep), 'of', length(keep), 'miRNAs\n')

# Run DESeq2
dds <- DESeq(dds)

# Get results with apeglm shrinkage
# apeglm provides better log2FC estimates for low-count genes
library(apeglm)
res <- lfcShrink(
    dds,
    coef = 'condition_treated_vs_control',
    type = 'apeglm'
)

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Filter significant miRNAs
# Standard thresholds:
# - padj < 0.05: FDR-corrected significance
# - |log2FC| > 1: 2-fold change (biologically meaningful)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1

sig <- subset(res, padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)

cat('Significant miRNAs:', nrow(sig), '\n')
cat('Upregulated:', sum(sig$log2FoldChange > 0, na.rm = TRUE), '\n')
cat('Downregulated:', sum(sig$log2FoldChange < 0, na.rm = TRUE), '\n')

# Top significant miRNAs
cat('\nTop 10 DE miRNAs:\n')
print(head(sig, 10))

# Visualization
library(ggplot2)
library(EnhancedVolcano)

# Volcano plot
EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = PADJ_CUTOFF,
    FCcutoff = LFC_CUTOFF,
    title = 'Differential miRNA Expression',
    subtitle = paste0('padj < ', PADJ_CUTOFF, ', |log2FC| > ', LFC_CUTOFF)
)
ggsave('volcano_plot.pdf', width = 10, height = 8)

# Heatmap of significant miRNAs
if (nrow(sig) > 0) {
    library(pheatmap)

    # Variance-stabilized counts for visualization
    vsd <- vst(dds, blind = FALSE)

    # Select significant miRNAs
    sig_mirnas <- rownames(sig)
    mat <- assay(vsd)[sig_mirnas, , drop = FALSE]

    # Z-score scale rows for better visualization
    mat_scaled <- t(scale(t(mat)))

    pheatmap(
        mat_scaled,
        annotation_col = coldata['condition'],
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = nrow(mat) < 50,
        main = 'Significant DE miRNAs'
    )
}

# Export results
res_df <- as.data.frame(res)
res_df$miRNA <- rownames(res_df)
write.csv(res_df, 'DE_miRNAs_all.csv', row.names = FALSE)
write.csv(as.data.frame(sig), 'DE_miRNAs_significant.csv')

cat('\nResults saved to DE_miRNAs_*.csv\n')
