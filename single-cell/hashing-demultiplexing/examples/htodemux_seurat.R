# Reference: Seurat 5.0+ | Verify API if version differs
# Hashtag demultiplexing with Seurat HTODemux: assign cells to samples and call cross-sample doublets.
# Inputs are illustrative; replace gex_counts and hto_counts with real matrices (cells must match).

library(Seurat)

# gex_counts: gene x cell sparse matrix; hto_counts: HTO x cell matrix with the SAME cell barcodes
obj <- CreateSeuratObject(counts = gex_counts)
obj[['HTO']] <- CreateAssay5Object(counts = hto_counts[, colnames(obj)])

# CLR margin=2 normalizes each tag across cells, correcting per-tag capture-efficiency differences
obj <- NormalizeData(obj, assay = 'HTO', normalization.method = 'CLR', margin = 2)

# positive.quantile=0.99: a cell is positive for a tag above the 0.99 quantile of that tag's
# inferred negative distribution; one positive = Singlet, two+ = Doublet, none = Negative
obj <- HTODemux(obj, assay = 'HTO', positive.quantile = 0.99)

global_counts <- table(obj$HTO_classification.global)   # Singlet / Doublet / Negative
sample_counts <- table(obj$hash.ID)                      # per-sample singlets + Doublet + Negative

# Cross-sample doublet fraction calibrates the expected total doublet rate; sanity-check it
doublet_rate <- global_counts['Doublet'] / sum(global_counts)

singlets <- subset(obj, subset = HTO_classification.global == 'Singlet')

# Ridge plots per tag reveal whether staining is bimodal (clean) or smeared (rescue with demuxmix)
RidgePlot(obj, assay = 'HTO', features = rownames(obj[['HTO']]), ncol = 2)

global_counts
doublet_rate
singlets
