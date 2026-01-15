# Create PCA plot from DESeq2 vst data

library(DESeq2)
library(ggplot2)

# Simulate data (in practice, use your own dds object)
set.seed(42)
n_genes <- 1000
n_samples <- 12

counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('gene', 1:n_genes),
                                paste0('sample', 1:n_samples)))

coldata <- data.frame(
    condition = factor(rep(c('control', 'treated'), 6)),
    batch = factor(rep(c('batch1', 'batch2', 'batch3'), each = 4)),
    row.names = colnames(counts)
)

# Create DESeq object and run vst
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ batch + condition)
vsd <- vst(dds, blind = FALSE)

# Get PCA data
pca_data <- plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, 'percentVar'))

# Create custom PCA plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
    geom_point(size = 4) +
    xlab(paste0('PC1: ', percentVar[1], '% variance')) +
    ylab(paste0('PC2: ', percentVar[2], '% variance')) +
    ggtitle('PCA of Samples') +
    theme_bw() +
    theme(legend.position = 'right')

print(p)
ggsave('pca_plot.png', p, width = 8, height = 6, dpi = 300)
