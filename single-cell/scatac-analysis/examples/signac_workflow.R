# scATAC-seq analysis with Signac
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load 10X data
counts <- Read10X_h5('filtered_peak_bc_matrix.h5')
metadata <- read.csv('singlecell.csv', header = TRUE, row.names = 1)

# Create Seurat object with chromatin assay
chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(':', '-'),
    genome = 'hg38',
    fragments = 'fragments.tsv.gz',
    min.cells = 10,
    min.features = 200
)
obj <- CreateSeuratObject(counts = chrom_assay, assay = 'peaks', meta.data = metadata)

# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(obj) <- annotations

# QC metrics
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj, fast = FALSE)
obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments

# QC plots
pdf('qc_plots.pdf', width = 12, height = 4)
VlnPlot(obj, features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal'), pt.size = 0.1, ncol = 4)
dev.off()

# Filter cells
obj <- subset(obj,
    peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

cat('Cells after QC:', ncol(obj), '\n')

# Normalization and LSI
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)

# Check depth correlation (skip first component if correlated)
depth_cor <- DepthCor(obj)
dims_use <- if(abs(depth_cor[1]) > 0.5) 2:30 else 1:30

# Clustering
obj <- RunUMAP(obj, reduction = 'lsi', dims = dims_use)
obj <- FindNeighbors(obj, reduction = 'lsi', dims = dims_use)
obj <- FindClusters(obj, algorithm = 3, resolution = 0.5)

# UMAP plot
pdf('umap.pdf', width = 8, height = 6)
DimPlot(obj, label = TRUE) + NoLegend()
dev.off()

# Gene activity scores
gene_activities <- GeneActivity(obj)
obj[['RNA']] <- CreateAssayObject(counts = gene_activities)
obj <- NormalizeData(obj, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = median(obj$nCount_RNA))

# Plot marker gene activities
DefaultAssay(obj) <- 'RNA'
pdf('gene_activity.pdf', width = 10, height = 8)
FeaturePlot(obj, features = c('CD34', 'MS4A1', 'CD3D', 'CD14', 'GATA1', 'PAX5'), pt.size = 0.1, max.cutoff = 'q95', ncol = 3)
dev.off()

# Save object
saveRDS(obj, 'scatac_processed.rds')
cat('Done. Saved to scatac_processed.rds\n')
