# DESeq2 with multiple conditions

library(DESeq2)

# Example with three conditions
set.seed(42)
n_genes <- 1000
n_samples <- 9

counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('gene', 1:n_genes),
                                paste0('sample', 1:n_samples)))

coldata <- data.frame(
    condition = factor(rep(c('control', 'treatmentA', 'treatmentB'), each = 3)),
    row.names = colnames(counts)
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)

# Pre-filter and set reference
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2
dds <- DESeq(dds)

# Check available coefficients
cat('Available coefficients:\n')
print(resultsNames(dds))

# Compare treatmentA vs control
res_A <- results(dds, name = 'condition_treatmentA_vs_control')
cat('\nTreatmentA vs Control:\n')
summary(res_A)

# Compare treatmentB vs control
res_B <- results(dds, name = 'condition_treatmentB_vs_control')
cat('\nTreatmentB vs Control:\n')
summary(res_B)

# Compare treatmentA vs treatmentB (using contrast)
res_AB <- results(dds, contrast = c('condition', 'treatmentA', 'treatmentB'))
cat('\nTreatmentA vs TreatmentB:\n')
summary(res_AB)
