# Cluster single-cell data with Seurat

library(Seurat)

seurat_obj <- readRDS('preprocessed.rds')
cat('Loaded', ncol(seurat_obj), 'cells\n')

seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
pdf('elbow_plot.pdf')
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
cat('Found', length(unique(Idents(seurat_obj))), 'clusters\n')
print(table(Idents(seurat_obj)))

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
pdf('umap_clusters.pdf')
DimPlot(seurat_obj, reduction = 'umap', label = TRUE)
dev.off()

saveRDS(seurat_obj, file = 'clustered.rds')
cat('Saved clustered data\n')
