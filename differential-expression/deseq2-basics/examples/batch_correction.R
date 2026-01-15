# DESeq2 with batch effect correction

library(DESeq2)
library(apeglm)

# Example with batch effects
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

# Create DESeqDataSet with batch in design
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ batch + condition)

# Pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2
dds <- DESeq(dds)

# Check available coefficients
cat('Available coefficients:\n')
print(resultsNames(dds))

# Get results for condition effect (controlling for batch)
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

summary(res)
