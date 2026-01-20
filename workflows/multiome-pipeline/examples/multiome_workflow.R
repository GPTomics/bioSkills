# Complete 10X Multiome workflow with Seurat and Signac

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)

# Configuration
data_dir <- 'cellranger_multiome_output'
output_dir <- 'multiome_results'
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, 'plots'), showWarnings = FALSE)

# === Step 1: Load Data ===
cat('=== Step 1: Loading Data ===\n')

# Read multiome data
counts <- Read10X_h5(file.path(data_dir, 'filtered_feature_bc_matrix.h5'))
frags_path <- file.path(data_dir, 'atac_fragments.tsv.gz')

# Create Seurat object with RNA
seurat_obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = 'RNA',
    min.cells = 3,
    min.features = 200
)

# Get annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

# Add ATAC assay
seurat_obj[['ATAC']] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(':', '-'),
    fragments = frags_path,
    annotation = annotations
)

cat('Initial cells:', ncol(seurat_obj), '\n')
cat('RNA features:', nrow(seurat_obj[['RNA']]), '\n')
cat('ATAC peaks:', nrow(seurat_obj[['ATAC']]), '\n')

# === Step 2: RNA QC ===
cat('\n=== Step 2: RNA QC ===\n')

seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
seurat_obj[['percent.ribo']] <- PercentageFeatureSet(seurat_obj, pattern = '^RP[SL]')

# QC plots
pdf(file.path(output_dir, 'plots', 'rna_qc.pdf'), width = 12, height = 4)
VlnPlot(seurat_obj, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), pt.size = 0, ncol = 3)
dev.off()

# Filter RNA
seurat_obj <- subset(seurat_obj,
    nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nFeature_RNA > 500 &
    percent.mt < 20
)
cat('After RNA QC:', ncol(seurat_obj), 'cells\n')

# === Step 3: ATAC QC ===
cat('\n=== Step 3: ATAC QC ===\n')

DefaultAssay(seurat_obj) <- 'ATAC'
seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj)

# QC plots
pdf(file.path(output_dir, 'plots', 'atac_qc.pdf'), width = 12, height = 4)
VlnPlot(seurat_obj, features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'), pt.size = 0, ncol = 3)
dev.off()

# Filter ATAC
seurat_obj <- subset(seurat_obj,
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    TSS.enrichment > 2 &
    nucleosome_signal < 4
)
cat('After ATAC QC:', ncol(seurat_obj), 'cells\n')

# === Step 4: RNA Processing ===
cat('\n=== Step 4: RNA Processing ===\n')

DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = 'percent.mt', verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

# === Step 5: ATAC Processing ===
cat('\n=== Step 5: ATAC Processing ===\n')

DefaultAssay(seurat_obj) <- 'ATAC'
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- RunSVD(seurat_obj, n = 50)

# Check depth correlation
pdf(file.path(output_dir, 'plots', 'lsi_depth_cor.pdf'), width = 8, height = 4)
DepthCor(seurat_obj)
dev.off()

# === Step 6: WNN Integration ===
cat('\n=== Step 6: WNN Integration ===\n')

seurat_obj <- FindMultiModalNeighbors(
    seurat_obj,
    reduction.list = list('pca', 'lsi'),
    dims.list = list(1:30, 2:30),
    modality.weight.name = 'RNA.weight'
)

seurat_obj <- RunUMAP(seurat_obj, nn.name = 'weighted.nn',
    reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')

seurat_obj <- FindClusters(seurat_obj, graph.name = 'wsnn',
    algorithm = 3, resolution = 0.5, verbose = FALSE)

cat('Clusters:', length(unique(seurat_obj$seurat_clusters)), '\n')

# === Step 7: Visualization ===
cat('\n=== Step 7: Visualization ===\n')

# Compare embeddings
p1 <- DimPlot(seurat_obj, reduction = 'pca', label = TRUE) + ggtitle('RNA (PCA)')
p2 <- DimPlot(seurat_obj, reduction = 'lsi', label = TRUE) + ggtitle('ATAC (LSI)')
p3 <- DimPlot(seurat_obj, reduction = 'wnn.umap', label = TRUE) + ggtitle('WNN (Joint)')

pdf(file.path(output_dir, 'plots', 'embeddings_comparison.pdf'), width = 15, height = 5)
print(p1 + p2 + p3)
dev.off()

# Modality weights
pdf(file.path(output_dir, 'plots', 'modality_weights.pdf'), width = 8, height = 4)
VlnPlot(seurat_obj, features = 'RNA.weight', group.by = 'seurat_clusters', pt.size = 0)
dev.off()

# === Step 8: Markers ===
cat('\n=== Step 8: Finding Markers ===\n')

# RNA markers
DefaultAssay(seurat_obj) <- 'SCT'
rna_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25,
    logfc.threshold = 0.25, verbose = FALSE)
write.csv(rna_markers, file.path(output_dir, 'rna_markers.csv'), row.names = FALSE)

# ATAC markers
DefaultAssay(seurat_obj) <- 'ATAC'
atac_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.05,
    test.use = 'LR', latent.vars = 'nCount_ATAC', verbose = FALSE)
write.csv(atac_markers, file.path(output_dir, 'atac_markers.csv'), row.names = FALSE)

# === Step 9: Save ===
cat('\n=== Step 9: Saving Results ===\n')
saveRDS(seurat_obj, file.path(output_dir, 'multiome_analyzed.rds'))

cat('\n=== Analysis Complete ===\n')
cat('Results saved to:', output_dir, '\n')
cat('  - Seurat object: multiome_analyzed.rds\n')
cat('  - RNA markers: rna_markers.csv\n')
cat('  - ATAC markers: atac_markers.csv\n')
cat('  - Plots: plots/\n')
