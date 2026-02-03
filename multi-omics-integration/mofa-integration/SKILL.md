---
name: bio-multi-omics-mofa-integration
description: Multi-Omics Factor Analysis (MOFA2) for unsupervised integration of multiple data modalities. Identifies shared and view-specific sources of variation. Use when integrating RNA-seq, proteomics, methylation, or other omics to discover latent factors driving biological variation across modalities.
tool_type: r
primary_tool: MOFA2
---

# MOFA2 Integration

## Prepare Multi-Omics Data

```r
library(MOFA2)
library(MultiAssayExperiment)

# Load individual omics matrices (samples x features)
rna <- as.matrix(read.csv('rnaseq_matrix.csv', row.names = 1))
protein <- as.matrix(read.csv('proteomics_matrix.csv', row.names = 1))
methylation <- as.matrix(read.csv('methylation_matrix.csv', row.names = 1))

# Ensure consistent sample names across views
common_samples <- Reduce(intersect, list(rownames(rna), rownames(protein), rownames(methylation)))
rna <- rna[common_samples, ]
protein <- protein[common_samples, ]
methylation <- methylation[common_samples, ]

# Transpose to features x samples (MOFA format)
data_list <- list(
    RNA = t(rna),
    Protein = t(protein),
    Methylation = t(methylation)
)
```

## Create and Train MOFA Model

```r
# Create MOFA object
mofa <- create_mofa(data_list)

# View data overview
plot_data_overview(mofa)

# Set model options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15  # Number of factors to learn

# Set training options
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- 'slow'
train_opts$seed <- 42

# Prepare and train
mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa, outfile = 'mofa_model.hdf5')
```

## Analyze Factors

```r
# Variance explained per factor per view
plot_variance_explained(mofa, max_r2 = 15)
plot_variance_explained(mofa, plot_total = TRUE)[[2]]

# Factor values (sample scores)
factors <- get_factors(mofa, as.data.frame = TRUE)

# Factor weights (feature loadings)
weights <- get_weights(mofa, as.data.frame = TRUE)
```

## Visualize Results

```r
# Scatter plot of factor values
plot_factor(mofa, factor = 1, color_by = 'Group')

# Heatmap of factor weights
plot_weights(mofa, view = 'RNA', factor = 1, nfeatures = 20)

# Top features per factor
plot_top_weights(mofa, view = 'RNA', factor = 1, nfeatures = 10)

# Correlation between factors
plot_factor_cor(mofa)
```

## Factor Interpretation

```r
# Extract top features for pathway analysis
top_rna_factor1 <- get_weights(mofa, views = 'RNA', factors = 1, as.data.frame = TRUE)
top_rna_factor1 <- top_rna_factor1[order(abs(top_rna_factor1$value), decreasing = TRUE), ]
gene_list <- head(top_rna_factor1$feature, 100)

# Gene set enrichment on factor weights
library(MOFA2)
enrichment <- run_enrichment(mofa, feature.sets = msigdb_genesets,
                              view = 'RNA', factors = 1:5)
plot_enrichment(enrichment, factor = 1, max.pathways = 15)
```

## Add Sample Metadata

```r
# Load sample annotations
metadata <- read.csv('sample_metadata.csv', row.names = 1)

# Add to MOFA object
samples_metadata(mofa) <- metadata[samples_names(mofa)[[1]], ]

# Color by metadata
plot_factor(mofa, factor = 1, color_by = 'Condition')
plot_factors(mofa, factors = 1:3, color_by = 'Condition')
```

## Multi-Group MOFA

```r
# For batch/group-structured data
data_list_grouped <- list(
    group1 = list(RNA = rna_g1, Protein = prot_g1),
    group2 = list(RNA = rna_g2, Protein = prot_g2)
)

mofa_grouped <- create_mofa(data_list_grouped)
mofa_grouped <- prepare_mofa(mofa_grouped)
mofa_grouped <- run_mofa(mofa_grouped)

# Compare factor activity across groups
plot_factor(mofa_grouped, factor = 1, group_by = 'group', color_by = 'group')
```

## MOFA+ for Single-Cell

```r
# For single-cell multi-omics (CITE-seq, Multiome)
library(Seurat)

# Extract modalities from Seurat object
rna_mat <- GetAssayData(seurat_obj, assay = 'RNA', layer = 'data')
adt_mat <- GetAssayData(seurat_obj, assay = 'ADT', layer = 'data')

# Create MOFA with single-cell settings
mofa_sc <- create_mofa(list(RNA = rna_mat, ADT = adt_mat))
model_opts <- get_default_model_options(mofa_sc)
model_opts$num_factors <- 10

# Use stochastic inference for large datasets
train_opts <- get_default_training_options(mofa_sc)
train_opts$stochastic <- TRUE
```

## Export Results

```r
# Save factor values
factors_df <- get_factors(mofa, as.data.frame = TRUE)
write.csv(factors_df, 'mofa_factors.csv', row.names = FALSE)

# Save weights
weights_df <- get_weights(mofa, as.data.frame = TRUE)
write.csv(weights_df, 'mofa_weights.csv', row.names = FALSE)

# Save variance explained
var_exp <- get_variance_explained(mofa)
write.csv(var_exp$r2_per_factor, 'mofa_variance_explained.csv')
```

## Related Skills

- mixomics-analysis - Supervised multi-omics integration
- data-harmonization - Preprocess data before MOFA
- pathway-analysis/go-enrichment - Interpret MOFA factors
- single-cell/multimodal-integration - Single-cell multi-omics
