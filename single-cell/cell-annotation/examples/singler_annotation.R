library(SingleR)
library(celldex)
library(Seurat)
library(SingleCellExperiment)
library(pheatmap)

seurat_obj <- readRDS('seurat_processed.rds')
sce <- as.SingleCellExperiment(seurat_obj)

ref_hpca <- celldex::HumanPrimaryCellAtlasData()
ref_monaco <- celldex::MonacoImmuneData()

pred_main <- SingleR(test = sce, ref = ref_hpca, labels = ref_hpca$label.main)
pred_fine <- SingleR(test = sce, ref = ref_hpca, labels = ref_hpca$label.fine)
pred_immune <- SingleR(test = sce, ref = ref_monaco, labels = ref_monaco$label.main)

seurat_obj$SingleR_main <- pred_main$labels
seurat_obj$SingleR_fine <- pred_fine$labels
seurat_obj$SingleR_immune <- pred_immune$labels
seurat_obj$SingleR_main_pruned <- pred_main$pruned.labels

pdf('singler_diagnostics.pdf', width = 12, height = 8)
plotScoreHeatmap(pred_main)
plotDeltaDistribution(pred_main, ncol = 4)
dev.off()

pdf('singler_umap.pdf', width = 14, height = 5)
p1 <- DimPlot(seurat_obj, group.by = 'SingleR_main', label = TRUE, repel = TRUE) +
    ggtitle('Main Labels') + NoLegend()
p2 <- DimPlot(seurat_obj, group.by = 'SingleR_immune', label = TRUE, repel = TRUE) +
    ggtitle('Immune Labels') + NoLegend()
print(p1 | p2)
dev.off()

label_counts <- table(seurat_obj$SingleR_main)
write.csv(as.data.frame(label_counts), 'singler_label_counts.csv', row.names = FALSE)

saveRDS(seurat_obj, 'seurat_annotated.rds')
