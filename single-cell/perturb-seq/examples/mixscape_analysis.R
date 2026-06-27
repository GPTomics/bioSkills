# Reference: Seurat 5.0+ | Verify API if version differs
# Perturb-seq escaper removal with Seurat Mixscape

library(Seurat)

# Assumes guide_ID and gene columns in metadata (gene = target, 'NT' = non-targeting)
seurat <- Read10X('filtered_feature_bc_matrix/')
seurat <- CreateSeuratObject(counts = seurat, project = 'perturb_seq')
guide_calls <- read.csv('guide_calls.csv', row.names = 1)
seurat <- AddMetaData(seurat, guide_calls)

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

# Local perturbation signature: subtract each cell's NT neighbors to cancel shared cell-state structure
seurat <- CalcPerturbSig(
    seurat,
    assay = 'RNA',
    slot = 'data',
    gd.class = 'guide_ID',
    nt.cell.class = 'NT',
    num.neighbors = 20,
    reduction = 'pca',
    ndims = 15,
    new.assay.name = 'PRTB'
)

# The slot bug: CalcPerturbSig writes PRTB 'data' but RunMixscape reads PRTB 'scale.data' -- scale it first
DefaultAssay(seurat) <- 'PRTB'
seurat <- ScaleData(seurat, assay = 'PRTB')

# Per-target 2-component mixture classifies KO vs NP (non-perturbed); NP cells are escapers to remove before DE
seurat <- RunMixscape(
    seurat,
    assay = 'PRTB',
    slot = 'scale.data',
    labels = 'gene',
    nt.class.name = 'NT',
    min.de.genes = 5,
    iter.num = 10,
    de.assay = 'RNA',
    prtb.type = 'KO'
)

# Report the perturbed fraction; an all-NP target is confounded with low guide efficiency, not proof of no function
table(seurat$mixscape_class.global)
