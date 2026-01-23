# Perturb-seq analysis with Seurat Mixscape

library(Seurat)

# Load data (assumes guide_ID in metadata)
seurat <- Read10X('filtered_feature_bc_matrix/')
seurat <- CreateSeuratObject(counts = seurat, project = 'perturb_seq')

# Add guide info
guide_calls <- read.csv('guide_calls.csv', row.names = 1)
seurat <- AddMetaData(seurat, guide_calls)

# Standard preprocessing
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:20)

# Mixscape: classify cells as perturbed or escaped
# Escaped cells received guide but show no perturbation effect
seurat <- CalcPerturbSig(
    seurat,
    assay = 'RNA',
    layer = 'data',
    gd.class = 'guide_ID',  # Guide identity column
    nt.cell.class = 'NT',   # Non-targeting label
    new.class.name = 'mixscape_class'
)

# Run Mixscape classification
seurat <- RunMixscape(
    seurat,
    assay = 'RNA',
    layer = 'scale.data',
    labels = 'gene',        # Target gene column
    nt.class.name = 'NT',
    min.de.genes = 5,       # Min DE genes to classify as perturbed
    iter.num = 10
)

# View classification results
table(seurat$mixscape_class)

# Perturbation vs NT differential expression for specific gene
de_results <- FindMarkers(
    seurat,
    ident.1 = 'GENE_KO',    # Target gene knockout
    ident.2 = 'NT',
    group.by = 'guide_ID'
)

# Save results
write.csv(de_results, 'de_GENE_vs_NT.csv')

# Visualize
DimPlot(seurat, group.by = 'mixscape_class')
ggsave('mixscape_classification.pdf')
