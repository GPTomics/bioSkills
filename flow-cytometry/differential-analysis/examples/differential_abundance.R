# Reference: diffcyt 1.22+, CATALYST 1.26+, edgeR 4.0+, limma 3.58+ | Verify API if version differs
library(CATALYST)
library(diffcyt)
library(SummarizedExperiment)

# Load clustered data
sce <- readRDS('sce_clustered.rds')
cat('Loaded', ncol(sce), 'cells from', length(unique(sce$sample_id)), 'samples\n')

# Check conditions
print(table(sce$condition))

# Create design matrix
design <- createDesignMatrix(ei(sce), cols_design = 'condition')

# Create contrast (Treatment vs Control)
contrast <- createContrast(c(0, 1))

# Differential abundance via the CATALYST-integrated diffcyt wrapper
# (sample is the unit: diffcyt aggregates cells to per-sample-per-cluster counts)
cat('\nRunning differential abundance analysis...\n')
res_DA <- diffcyt(sce, clustering_to_use = 'meta20',
                  analysis_type = 'DA', method_DA = 'diffcyt-DA-edgeR',
                  design = design, contrast = contrast)

# Results (rowData of the inner result object; p_adj is BH across clusters)
da_results <- as.data.frame(rowData(res_DA$res))
da_results <- da_results[order(da_results$p_adj), ]

cat('\nSignificant clusters (FDR < 0.05):\n')
sig <- da_results[da_results$p_adj < 0.05, ]
print(sig[, c('cluster_id', 'logFC', 'p_val', 'p_adj')])

# Save
write.csv(da_results, 'da_results.csv', row.names = FALSE)
cat('\nResults saved to da_results.csv\n')
