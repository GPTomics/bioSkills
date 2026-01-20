library(CATALYST)
library(diffcyt)

# Load clustered data
sce <- readRDS('sce_clustered.rds')
cat('Loaded', ncol(sce), 'cells from', length(unique(sce$sample_id)), 'samples\n')

# Check conditions
print(table(sce$condition))

# Create design matrix
design <- createDesignMatrix(ei(sce), cols_design = 'condition')

# Create contrast (Treatment vs Control)
contrast <- createContrast(c(0, 1))

# Differential abundance
cat('\nRunning differential abundance analysis...\n')
res_DA <- testDA_edgeR(sce, design, contrast, cluster_id = 'meta20')

# Results
da_results <- as.data.frame(rowData(res_DA))
da_results <- da_results[order(da_results$p_adj), ]

cat('\nSignificant clusters (FDR < 0.05):\n')
sig <- da_results[da_results$p_adj < 0.05, ]
print(sig[, c('cluster_id', 'logFC', 'p_val', 'p_adj')])

# Save
write.csv(da_results, 'da_results.csv', row.names = FALSE)
cat('\nResults saved to da_results.csv\n')
