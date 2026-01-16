# Basic DESeq2 workflow for differential expression analysis

library(DESeq2)
library(apeglm)

# Simulate example data
set.seed(42)
n_genes <- 1000
n_samples <- 6

counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('gene', 1:n_genes),
                                paste0('sample', 1:n_samples)))

coldata <- data.frame(
    condition = factor(rep(c('control', 'treated'), each = 3)),
    row.names = colnames(counts)
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat('Genes after filtering:', nrow(dds), '\n')

# Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Get results with shrinkage
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# Summary
summary(res)

# Get significant genes
sig_genes <- subset(res, padj < 0.05)
cat('\nSignificant genes (padj < 0.05):', nrow(sig_genes), '\n')
