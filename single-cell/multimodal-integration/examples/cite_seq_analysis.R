library(Seurat)

data <- Read10X('filtered_feature_bc_matrix/')
rna_counts <- data$`Gene Expression`
adt_counts <- data$`Antibody Capture`

obj <- CreateSeuratObject(counts = rna_counts, assay = 'RNA', min.cells = 3)
obj[['ADT']] <- CreateAssayObject(counts = adt_counts[, colnames(obj)])

obj <- PercentageFeatureSet(obj, pattern = '^MT-', col.name = 'percent.mt')
obj <- subset(obj, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

obj <- NormalizeData(obj, assay = 'RNA')
obj <- FindVariableFeatures(obj, assay = 'RNA')
obj <- ScaleData(obj, assay = 'RNA')
obj <- RunPCA(obj, assay = 'RNA', reduction.name = 'pca')

obj <- NormalizeData(obj, assay = 'ADT', normalization.method = 'CLR', margin = 2)
obj <- ScaleData(obj, assay = 'ADT')
obj <- RunPCA(obj, assay = 'ADT', reduction.name = 'apca',
              features = rownames(obj[['ADT']]), npcs = min(18, nrow(obj[['ADT']]) - 1))

obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'),
                                dims.list = list(1:30, 1:min(18, ncol(obj[['apca']]))))
obj <- FindClusters(obj, graph.name = 'wsnn', resolution = 0.5)
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')

pdf('cite_seq_wnn_umap.pdf', width = 10, height = 8)
print(DimPlot(obj, reduction = 'wnn.umap', label = TRUE))
dev.off()

saveRDS(obj, 'cite_seq_analyzed.rds')
cat('Analysis complete. Saved to cite_seq_analyzed.rds\n')
