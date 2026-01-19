library(Seurat)
library(DoubletFinder)

args <- commandArgs(trailingOnly = TRUE)
input_path <- ifelse(length(args) > 0, args[1], 'filtered_feature_bc_matrix/')

counts <- Read10X(data.dir = input_path)
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

cat('Loaded', ncol(seurat_obj), 'cells\n')

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

sweep.res <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
cat('Optimal pK:', optimal_pk, '\n')

nExp_poi <- round(0.06 * ncol(seurat_obj))
cat('Expected doublets:', nExp_poi, '\n')

seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, pK = optimal_pk,
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

df_col <- grep('DF.classifications', colnames(seurat_obj@meta.data), value = TRUE)
n_doublets <- sum(seurat_obj@meta.data[[df_col]] == 'Doublet')
cat('Detected', n_doublets, 'doublets\n')

pdf('doublet_umap.pdf', width = 8, height = 6)
print(DimPlot(seurat_obj, group.by = df_col))
dev.off()

seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[df_col]] == 'Singlet'])
cat('Kept', ncol(seurat_obj), 'singlets\n')

saveRDS(seurat_obj, 'singlets.rds')
