# Complete single-cell RNA-seq workflow with Seurat

library(Seurat)
library(scDblFinder)
library(ggplot2)
library(dplyr)

# Public scRNA-seq datasets:
# - 10x Genomics: https://www.10xgenomics.com/datasets (PBMC 3k, 10k)
# - SeuratData package: InstallData('pbmc3k') for built-in tutorial data
# - CELLxGENE: https://cellxgene.cziscience.com (curated annotated datasets)
# - GEO: GSE149173 (COVID-19 PBMC), GSE126030 (multi-tissue)
# - Human Cell Atlas: https://data.humancellatlas.org

# Configuration
data_dir <- 'filtered_feature_bc_matrix'
output_dir <- 'scrnaseq_results'
project_name <- 'my_scrnaseq'

dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, 'plots'), showWarnings = FALSE)

# === Step 1: Load Data ===
cat('Loading data...\n')
counts <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = counts, project = project_name,
                                  min.cells = 3, min.features = 200)
cat('Initial cells:', ncol(seurat_obj), 'Genes:', nrow(seurat_obj), '\n')

# === Step 2: QC Metrics ===
cat('Calculating QC metrics...\n')
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
seurat_obj[['percent.ribo']] <- PercentageFeatureSet(seurat_obj, pattern = '^RP[SL]')

# QC violin plots
pdf(file.path(output_dir, 'plots', 'qc_violin.pdf'), width = 12, height = 4)
VlnPlot(seurat_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
dev.off()

# QC scatter plots
pdf(file.path(output_dir, 'plots', 'qc_scatter.pdf'), width = 10, height = 4)
p1 <- FeatureScatter(seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
p2 <- FeatureScatter(seurat_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
p1 + p2
dev.off()

# Filter cells
seurat_obj <- subset(seurat_obj,
                     nFeature_RNA > 200 &
                     nFeature_RNA < 5000 &
                     percent.mt < 20 &
                     nCount_RNA > 500)
cat('After QC filtering:', ncol(seurat_obj), 'cells\n')

# === Step 3: Doublet Detection ===
cat('Detecting doublets...\n')
sce <- as.SingleCellExperiment(seurat_obj)
set.seed(42)
sce <- scDblFinder(sce)

seurat_obj$doublet_class <- sce$scDblFinder.class
seurat_obj$doublet_score <- sce$scDblFinder.score

doublet_rate <- sum(seurat_obj$doublet_class == 'doublet') / ncol(seurat_obj)
cat('Doublet rate:', round(doublet_rate * 100, 1), '%\n')

seurat_obj <- subset(seurat_obj, doublet_class == 'singlet')
cat('After doublet removal:', ncol(seurat_obj), 'cells\n')

# === Step 4: Normalization ===
cat('Running SCTransform...\n')
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = 'percent.mt', verbose = FALSE)

# === Step 5: Dimensionality Reduction ===
cat('Running PCA...\n')
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

# Elbow plot
pdf(file.path(output_dir, 'plots', 'elbow.pdf'), width = 6, height = 4)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

# UMAP
n_pcs <- 30
cat('Running UMAP with', n_pcs, 'PCs...\n')
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)

# === Step 6: Clustering ===
cat('Finding clusters...\n')
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0), verbose = FALSE)

# UMAP plots at different resolutions
pdf(file.path(output_dir, 'plots', 'umap_resolutions.pdf'), width = 15, height = 10)
p1 <- DimPlot(seurat_obj, group.by = 'SCT_snn_res.0.2', label = TRUE) + ggtitle('res=0.2')
p2 <- DimPlot(seurat_obj, group.by = 'SCT_snn_res.0.4', label = TRUE) + ggtitle('res=0.4')
p3 <- DimPlot(seurat_obj, group.by = 'SCT_snn_res.0.6', label = TRUE) + ggtitle('res=0.6')
p4 <- DimPlot(seurat_obj, group.by = 'SCT_snn_res.0.8', label = TRUE) + ggtitle('res=0.8')
(p1 | p2) / (p3 | p4)
dev.off()

# Set default resolution
Idents(seurat_obj) <- 'SCT_snn_res.0.5'
seurat_obj$seurat_clusters <- Idents(seurat_obj)

# === Step 7: Find Markers ===
cat('Finding marker genes...\n')
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)

# Top markers per cluster
top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

write.csv(markers, file.path(output_dir, 'all_markers.csv'), row.names = FALSE)
write.csv(top_markers, file.path(output_dir, 'top10_markers.csv'), row.names = FALSE)

# Heatmap
pdf(file.path(output_dir, 'plots', 'marker_heatmap.pdf'), width = 12, height = 10)
top5 <- markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
DoHeatmap(seurat_obj, features = top5$gene) + NoLegend()
dev.off()

# === Step 8: Save Results ===
saveRDS(seurat_obj, file.path(output_dir, 'seurat_object.rds'))

# Summary
cat('\n=== Pipeline Complete ===\n')
cat('Final cells:', ncol(seurat_obj), '\n')
cat('Clusters:', length(unique(Idents(seurat_obj))), '\n')
cat('Markers found:', nrow(markers), '\n')
cat('Results saved to:', output_dir, '\n')
cat('\nTop 3 markers per cluster:\n')
print(markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC) %>%
      select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj))
