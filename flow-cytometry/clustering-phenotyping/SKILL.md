---
name: bio-flow-cytometry-clustering-phenotyping
description: Unsupervised clustering and cell type identification for flow/mass cytometry. Covers FlowSOM, Phenograph, and CATALYST workflows. Use when discovering cell populations in high-dimensional cytometry data without predefined gates.
tool_type: r
primary_tool: CATALYST
---

# Clustering and Phenotyping

## FlowSOM Clustering

```r
library(FlowSOM)

# Prepare data
expr <- exprs(fcs)
marker_cols <- grep('CD|HLA', colnames(fcs), value = TRUE)

# Build SOM
fsom <- FlowSOM(fcs,
                colsToUse = marker_cols,
                xdim = 10, ydim = 10,
                nClus = 20,
                seed = 42)

# Get cluster assignments
clusters <- GetMetaclusters(fsom)

# Add to flowFrame
exprs(fcs) <- cbind(exprs(fcs), cluster = clusters)
```

## CATALYST Workflow (Full Pipeline)

```r
library(CATALYST)
library(SingleCellExperiment)

# Create SCE from flowSet
sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 5)

# Clustering
sce <- cluster(sce,
               features = 'type',  # Use 'type' markers from panel
               xdim = 10, ydim = 10,
               maxK = 20,
               seed = 42)

# View cluster assignments
table(cluster_ids(sce, 'meta20'))
```

## Phenograph Clustering

```r
library(Rphenograph)

# Extract expression matrix
expr <- assay(sce, 'exprs')

# Run Phenograph
pheno_result <- Rphenograph(t(expr[rowData(sce)$marker_class == 'type', ]), k = 30)

# Get clusters
sce$phenograph <- factor(membership(pheno_result[[2]]))
```

## Dimensionality Reduction

```r
# UMAP
sce <- runDR(sce, dr = 'UMAP', features = 'type')

# tSNE
sce <- runDR(sce, dr = 'TSNE', features = 'type')

# Plot
plotDR(sce, 'UMAP', color_by = 'meta20')
```

## Cluster Annotation

```r
# Heatmap of marker expression by cluster
plotExprHeatmap(sce, features = 'type',
                by = 'cluster_id', k = 'meta20',
                scale = 'first', row_anno = FALSE)

# Manual annotation
cluster_annotation <- c(
    '1' = 'CD4 T cells',
    '2' = 'CD8 T cells',
    '3' = 'B cells',
    '4' = 'NK cells',
    '5' = 'Monocytes'
)

sce$cell_type <- cluster_annotation[as.character(cluster_ids(sce, 'meta20'))]
```

## Cluster Merging

```r
# Merge similar clusters
merging_table <- data.frame(
    original = 1:20,
    merged = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
               6, 6, 7, 7, 8, 8, 9, 9, 10, 10)
)

sce <- mergeClusters(sce, k = 'meta20', table = merging_table, id = 'merged')
```

## Abundance Analysis (per sample)

```r
# Cluster frequencies per sample
abundances <- table(cluster_ids(sce, 'meta20'), sce$sample_id)
freq <- prop.table(abundances, margin = 2)

# Plot
plotAbundances(sce, k = 'meta20', by = 'cluster_id', group_by = 'condition')
```

## Marker Expression Summary

```r
# Median expression per cluster
plotClusterExprs(sce, k = 'meta20', features = 'type')

# Expression by cluster and condition
plotPbExprs(sce, k = 'meta20', features = 'type', facet_by = 'cluster_id')
```

## Export Results

```r
# Add cluster info to metadata
colData(sce)$cluster <- cluster_ids(sce, 'meta20')

# Export to CSV
results <- as.data.frame(colData(sce))
write.csv(results, 'clustering_results.csv', row.names = FALSE)

# Save SCE
saveRDS(sce, 'sce_clustered.rds')
```

## Choosing Number of Clusters

```r
# Delta area plot
plotNRS(sce, features = 'type')

# Or visual inspection of heatmap at different K
plotExprHeatmap(sce, features = 'type', by = 'cluster_id', k = 'meta10')
plotExprHeatmap(sce, features = 'type', by = 'cluster_id', k = 'meta20')
```

## Batch Integration

```r
# If batch effects present
library(batchelor)

sce <- runDR(sce, dr = 'UMAP', features = 'type')

# Check for batch effects
plotDR(sce, 'UMAP', color_by = 'batch')

# MNN correction if needed
sce_corrected <- fastMNN(sce, batch = sce$batch)
```

## Related Skills

- gating-analysis - Manual alternative
- differential-analysis - Compare clusters between conditions
- single-cell/clustering - Similar concepts for scRNA-seq
