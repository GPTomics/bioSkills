#!/usr/bin/env Rscript
# Batch integration with Harmony

library(Seurat)
library(harmony)
library(ggplot2)

# Example: integrate multiple 10X samples
# Assumes sample1.rds, sample2.rds, sample3.rds exist as Seurat objects

integrate_with_harmony <- function(sample_files, output_prefix = 'integrated') {
    cat('Loading samples...\n')
    samples <- lapply(sample_files, readRDS)
    names(samples) <- gsub('\\.rds$', '', basename(sample_files))

    cat('Merging samples...\n')
    merged <- merge(samples[[1]], y = samples[-1], add.cell.ids = names(samples))

    cat('Standard preprocessing...\n')
    merged <- NormalizeData(merged)
    merged <- FindVariableFeatures(merged, selection.method = 'vst', nfeatures = 2000)
    merged <- ScaleData(merged)
    merged <- RunPCA(merged, npcs = 50)

    cat('Running Harmony integration...\n')
    merged <- RunHarmony(merged, group.by.vars = 'orig.ident', dims.use = 1:30)

    cat('Post-integration analysis...\n')
    merged <- RunUMAP(merged, reduction = 'harmony', dims = 1:30)
    merged <- FindNeighbors(merged, reduction = 'harmony', dims = 1:30)
    merged <- FindClusters(merged, resolution = 0.5)

    # Visualizations
    p1 <- DimPlot(merged, reduction = 'umap', group.by = 'orig.ident') +
        ggtitle('By Sample')
    p2 <- DimPlot(merged, reduction = 'umap', label = TRUE) +
        ggtitle('By Cluster')

    pdf(paste0(output_prefix, '_umap.pdf'), width = 14, height = 6)
    print(p1 + p2)
    dev.off()

    saveRDS(merged, paste0(output_prefix, '.rds'))
    cat('Saved:', paste0(output_prefix, '.rds\n'))

    return(merged)
}

# Run if executed directly
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    sample_files <- args
    integrate_with_harmony(sample_files)
} else {
    cat('Usage: Rscript harmony_integration.R sample1.rds sample2.rds sample3.rds\n')
}
