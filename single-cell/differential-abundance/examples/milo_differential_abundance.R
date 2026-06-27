# Reference: miloR 2.0+ | Verify API if version differs
# Cluster-free differential abundance on kNN-graph neighborhoods.
# Expects a SingleCellExperiment 'sce' with an integrated 'PCA' reducedDim,
# a 'sample' column (biological replicate) and a 'condition' column in colData.
library(miloR)
library(SingleCellExperiment)
library(dplyr)

milo <- Milo(sce)

# k and d set neighborhood size; prop samples representative neighborhoods.
# k=30 balances power against resolution; refined sampling reduces redundant nhoods.
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = 'PCA')
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = 'PCA')

milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = 'sample')
milo <- calcNhoodDistance(milo, d = 30, reduced.dim = 'PCA')

design <- distinct(as.data.frame(colData(milo))[, c('sample', 'condition')])
rownames(design) <- design$sample

# GLM on neighborhood counts; SpatialFDR corrects for overlapping neighborhoods.
da <- testNhoods(milo, design = ~ condition, design.df = design, reduced.dim = 'PCA')
da <- annotateNhoods(milo, da, coldata_col = 'cell_type')

# SpatialFDR < 0.1 is the reportable significance, not raw PValue.
sig <- da[da$SpatialFDR < 0.1, ]
print(table(sign(sig$logFC), sig$cell_type))
print(head(da[order(da$SpatialFDR), c('Nhood', 'logFC', 'SpatialFDR', 'cell_type')], 20))
