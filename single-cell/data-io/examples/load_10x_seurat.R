# Reference: Cell Ranger 8.0+, anndata 0.10+, numpy 1.26+, pandas 2.2+, scanpy 1.10+ | Verify API if version differs
# Load 10X Genomics data with Seurat

library(Seurat)

counts <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
seurat_obj <- CreateSeuratObject(counts = counts, project = 'PBMC', min.cells = 3, min.features = 200)

print(seurat_obj)
cat('Cells:', ncol(seurat_obj), '\n')
cat('Genes:', nrow(seurat_obj), '\n')

saveRDS(seurat_obj, file = 'pbmc.rds')
cat('Saved to pbmc.rds\n')
