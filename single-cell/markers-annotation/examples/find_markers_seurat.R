# Reference: Seurat 5.0+ | Verify API if version differs
# Find cluster marker genes with Seurat. Install presto so Wilcoxon is fast;
# without it Seurat silently falls back to slow base-R. Seurat v5 lowered the
# defaults to logfc.threshold=0.1 / min.pct=0.01, so an explicit higher threshold
# plus a pct.1-pct.2 specificity filter restores a marker-grade shortlist.

library(Seurat)
library(dplyr)

seurat_obj <- readRDS('clustered.rds')

all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)

top_markers <- all_markers %>%
    filter(p_val_adj < 0.05, avg_log2FC > 1, (pct.1 - pct.2) > 0.2) %>%
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
