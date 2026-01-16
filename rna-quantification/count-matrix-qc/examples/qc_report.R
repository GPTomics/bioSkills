library(DESeq2)
library(pheatmap)
library(ggplot2)

counts <- read.csv('count_matrix.csv', row.names = 1)
coldata <- read.csv('sample_info.csv', row.names = 1)
coldata$condition <- factor(coldata$condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

cat('=== RNA-seq QC Report ===\n\n')
cat('Samples:', ncol(dds), '\n')
cat('Genes:', nrow(dds), '\n\n')

cat('Library sizes:\n')
lib_sizes <- colSums(counts(dds))
print(summary(lib_sizes))
cat('\n')

cat('Genes detected per sample:\n')
genes_detected <- colSums(counts(dds) > 0)
print(summary(genes_detected))
cat('\n')

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat('Genes after filtering:', nrow(dds), '\n\n')

vsd <- vst(dds, blind = TRUE)

png('qc_sample_correlation.png', width = 800, height = 800)
pheatmap(cor(assay(vsd)), annotation_col = coldata['condition'], main = 'Sample Correlation')
dev.off()

pca_data <- plotPCA(vsd, intgroup = 'condition', returnData = TRUE)
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    theme_minimal() +
    ggtitle('PCA of Samples')
ggsave('qc_pca.png', pca_plot, width = 8, height = 6)

png('qc_library_sizes.png', width = 800, height = 400)
barplot(lib_sizes / 1e6, las = 2, ylab = 'Million Reads', main = 'Library Sizes')
dev.off()

cat('QC plots saved: qc_sample_correlation.png, qc_pca.png, qc_library_sizes.png\n')
