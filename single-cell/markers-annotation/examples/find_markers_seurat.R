# Find marker genes with Seurat

library(Seurat)
library(dplyr)

seurat_obj <- readRDS('clustered.rds')

all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
print(top_markers)

write.csv(all_markers, file = 'all_markers.csv', row.names = FALSE)

pbmc_markers <- c('CD3D', 'CD8A', 'MS4A1', 'CD14', 'FCGR3A', 'NKG7')
pdf('dotplot_markers.pdf', width = 10, height = 6)
DotPlot(seurat_obj, features = pbmc_markers) + RotatedAxis()
dev.off()

new_cluster_ids <- c('0' = 'T cells', '1' = 'Monocytes', '2' = 'B cells')
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
seurat_obj$cell_type <- Idents(seurat_obj)

pdf('umap_celltypes.pdf')
DimPlot(seurat_obj, reduction = 'umap', label = TRUE)
dev.off()

saveRDS(seurat_obj, file = 'annotated.rds')
cat('Saved annotated data\n')
