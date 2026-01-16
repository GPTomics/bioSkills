# Preprocess single-cell data with Seurat

library(Seurat)

counts <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
cat('Raw:', ncol(seurat_obj), 'cells,', nrow(seurat_obj), 'genes\n')

seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')

seurat_obj <- subset(seurat_obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
cat('Filtered:', ncol(seurat_obj), 'cells\n')

seurat_obj <- SCTransform(seurat_obj, vars.to.regress = 'percent.mt', verbose = FALSE)
cat('HVGs:', length(VariableFeatures(seurat_obj)), '\n')

saveRDS(seurat_obj, file = 'preprocessed.rds')
cat('Saved preprocessed data\n')
