---
name: bio-differential-expression-batch-correction
description: Remove batch effects from RNA-seq data using ComBat, ComBat-Seq, limma removeBatchEffect, and SVA for unknown batch variables. Use when correcting batch effects in expression data.
tool_type: r
primary_tool: sva
---

# Batch Effect Correction

## ComBat-Seq (Count Data)

```r
library(sva)

# counts: raw count matrix (genes x samples)
# batch: vector of batch labels
# group: vector of biological condition (optional, to preserve)

corrected_counts <- ComBat_seq(counts = as.matrix(counts),
                                batch = batch,
                                group = condition,
                                full_mod = TRUE)

# Result is batch-corrected count matrix
# Use for visualization, clustering, but NOT for DE (use design formula instead)
```

## ComBat (Normalized Data)

```r
library(sva)

# For normalized expression (log-transformed, TPM, etc.)
# NOT for raw counts

# Create model matrix
mod <- model.matrix(~ condition, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)

# Run ComBat
corrected_expr <- ComBat(dat = as.matrix(normalized_expr),
                          batch = metadata$batch,
                          mod = mod,
                          par.prior = TRUE)
```

## limma removeBatchEffect

```r
library(limma)

# For visualization/clustering only
# Preserves group differences while removing batch

design <- model.matrix(~ condition, data = metadata)
corrected_expr <- removeBatchEffect(normalized_expr,
                                     batch = metadata$batch,
                                     design = design)

# For PCA, heatmaps, etc.
```

## DESeq2 Design Formula (Recommended for DE)

```r
library(DESeq2)

# Include batch in design formula - preferred for DE analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ batch + condition)

# Batch is modeled, not removed
# DE results are adjusted for batch
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'treatment', 'control'))
```

## Surrogate Variable Analysis (SVA)

```r
library(sva)

# When batch is unknown, estimate surrogate variables
mod <- model.matrix(~ condition, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)

# Estimate number of surrogate variables
n_sv <- num.sv(normalized_expr, mod, method = 'leek')

# Estimate surrogate variables
svobj <- sva(normalized_expr, mod, mod0, n.sv = n_sv)

# Add SVs to design for DE
design_with_sv <- cbind(mod, svobj$sv)
```

## SVA with DESeq2

```r
library(DESeq2)
library(sva)

# Normalize for SV estimation
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Estimate SVs
mod <- model.matrix(~ condition, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)
svobj <- sva(norm_counts, mod, mod0)

# Add SVs to colData
for (i in seq_len(ncol(svobj$sv))) {
    colData(dds)[[paste0('SV', i)]] <- svobj$sv[, i]
}

# Update design
sv_formula <- as.formula(paste('~', paste(paste0('SV', 1:ncol(svobj$sv)), collapse = ' + '), '+ condition'))
design(dds) <- sv_formula

# Run DESeq2
dds <- DESeq(dds)
```

## Visualize Batch Effects

```r
library(ggplot2)

# PCA before correction
pca_before <- prcomp(t(normalized_expr), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_before$x[, 1], PC2 = pca_before$x[, 2],
                     batch = metadata$batch, condition = metadata$condition)

p1 <- ggplot(pca_df, aes(PC1, PC2, color = batch, shape = condition)) +
    geom_point(size = 3) + ggtitle('Before Correction')

# PCA after correction
pca_after <- prcomp(t(corrected_expr), scale. = TRUE)
pca_df_after <- data.frame(PC1 = pca_after$x[, 1], PC2 = pca_after$x[, 2],
                           batch = metadata$batch, condition = metadata$condition)

p2 <- ggplot(pca_df_after, aes(PC1, PC2, color = batch, shape = condition)) +
    geom_point(size = 3) + ggtitle('After Correction')

library(patchwork)
p1 + p2
```

## Quantify Batch Effect

```r
# PVCA - Principal Variance Component Analysis
library(pvca)

# Proportion of variance explained by batch vs condition
pvcaObj <- pvcaBatchAssess(normalized_expr, metadata, threshold = 0.6,
                            theInteractionTerms = c('batch', 'condition'))

# Or manual approach
pca <- prcomp(t(normalized_expr), scale. = TRUE)
variance_explained <- summary(pca)$importance[2, 1:5]

# Correlation of PCs with batch
cor(pca$x[, 1], as.numeric(as.factor(metadata$batch)))
```

## Harmony (Single-Cell Integration)

```r
library(harmony)
library(Seurat)

# For single-cell data with multiple batches
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = 'batch', reduction = 'pca',
                          dims.use = 1:30)

# Use harmony reduction for downstream
seurat_obj <- RunUMAP(seurat_obj, reduction = 'harmony', dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = 'harmony', dims = 1:30)
```

## When NOT to Correct

```r
# DON'T use batch-corrected values for:
# - Differential expression (use design formula instead)
# - Count-based methods expecting raw/normalized counts

# DO use batch-corrected values for:
# - Visualization (PCA, UMAP, heatmaps)
# - Clustering
# - Machine learning features
# - Cross-study comparisons
```

## Related Skills

- differential-expression/deseq2-basics - DE with batch in design
- single-cell/clustering - Integration methods
- expression-matrix/matrix-operations - Data transformation
