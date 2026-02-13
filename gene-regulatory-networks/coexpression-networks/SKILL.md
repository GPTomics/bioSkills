---
name: bio-gene-regulatory-networks-coexpression-networks
description: Build weighted gene co-expression networks to identify modules of co-regulated genes and relate them to phenotypes using WGCNA and CEMiTool. Detects hub genes and module-trait relationships from bulk or single-cell expression data. Use when finding co-expression modules, identifying hub genes, or relating gene networks to clinical or experimental variables.
tool_type: r
primary_tool: WGCNA
---

## Version Compatibility

Reference examples tested with: WGCNA 1.72+, CEMiTool 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion("<pkg>")` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Co-expression Networks

**"Find co-expression modules and hub genes from my RNA-seq data"** → Build a weighted gene co-expression network, detect modules of co-regulated genes via hierarchical clustering, and correlate modules with sample traits to identify hub genes.
- R: `WGCNA::blockwiseModules()` for network construction and module detection
- R: `CEMiTool::cemitool()` for automated co-expression analysis

Build weighted gene co-expression networks to identify modules of co-regulated genes and relate them to sample traits.

## WGCNA Workflow

### Required Libraries

```r
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
```

### Input Preparation

```r
# Expression matrix: genes as columns, samples as rows (WGCNA convention)
expr_data <- read.csv('normalized_counts.csv', row.names = 1)
expr_data <- t(expr_data)  # transpose if genes are rows

# Filter low-variance genes
gene_vars <- apply(expr_data, 2, var)
expr_data <- expr_data[, gene_vars > quantile(gene_vars, 0.25)]

# Check for outlier samples
sample_tree <- hclust(dist(expr_data), method = 'average')
plot(sample_tree, main = 'Sample dendrogram')
```

### Soft-Thresholding Power Selection

```r
powers <- c(1:20)
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5)

# Plot scale-free topology fit
# R^2 > 0.85 is the default target; acceptable range is 0.8-0.9
# Higher R^2 means better scale-free fit; if no power reaches 0.85,
# use the power where R^2 plateaus (but should be >= 0.80)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology R^2',
     main = 'Scale independence')
abline(h = 0.85, col = 'red')

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity',
     main = 'Mean connectivity')

soft_power <- sft$powerEstimate
```

### Network Construction and Module Detection

```r
net <- blockwiseModules(
    expr_data, power = soft_power,
    TOMType = 'unsigned', minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE, saveTOMFileBase = 'TOM',
    verbose = 3
)

module_colors <- labels2colors(net$colors)
table(module_colors)

plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                    'Module colors', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

### Module Eigengenes and Trait Relationships

```r
# Module eigengenes (first PC of each module)
MEs <- net$MEs
MEs <- orderMEs(MEs)

# Trait data (samples as rows, traits as columns)
traits <- read.csv('sample_traits.csv', row.names = 1)

# Correlate eigengenes with traits
module_trait_cor <- cor(MEs, traits, use = 'p')
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_data))

# Heatmap of module-trait correlations
textMatrix <- paste(signif(module_trait_cor, 2), '\n(',
                    signif(module_trait_pval, 1), ')', sep = '')
dim(textMatrix) <- dim(module_trait_cor)

labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(traits), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.5,
               main = 'Module-trait relationships')
```

### Hub Gene Identification

**Goal:** Identify the most influential genes within a co-expression module by combining module membership with trait association.

**Approach:** Correlate each gene's expression with its module eigengene (module membership) and with the trait of interest (gene significance), then select genes exceeding both thresholds; additionally rank by intramodular connectivity (kWithin).

```r
# Gene module membership (correlation with module eigengene)
module_of_interest <- 'turquoise'
module_genes <- colnames(expr_data)[module_colors == module_of_interest]

gene_module_membership <- cor(expr_data, MEs, use = 'p')
gene_trait_significance <- cor(expr_data, traits$phenotype, use = 'p')

# Hub genes: high module membership AND high trait significance
hub_threshold_mm <- 0.8   # module membership cutoff
hub_threshold_gs <- 0.2   # gene significance cutoff
hub_genes <- module_genes[
    abs(gene_module_membership[module_genes, paste0('ME', module_of_interest)]) > hub_threshold_mm &
    abs(gene_trait_significance[module_genes, 1]) > hub_threshold_gs
]

# Intramodular connectivity (kWithin)
connectivity <- intramodularConnectivity(
    adjacency(expr_data, power = soft_power),
    module_colors
)
top_hubs <- connectivity[module_genes, ] %>%
    dplyr::arrange(desc(kWithin)) %>%
    head(20)
```

### Export for Cytoscape

**Goal:** Export a co-expression module as an edge/node list for interactive network visualization in Cytoscape.

**Approach:** Recompute the TOM (Topological Overlap Matrix) for the module of interest, threshold weak connections, and export as Cytoscape-compatible edge and node files.

```r
# Export top connections for network visualization
TOM <- TOMsimilarityFromExpr(expr_data, power = soft_power)
dimnames(TOM) <- list(colnames(expr_data), colnames(expr_data))

module_genes <- colnames(expr_data)[module_colors == module_of_interest]
module_TOM <- TOM[module_genes, module_genes]

cyt <- exportNetworkToCytoscape(module_TOM, edgeFile = 'edges.txt',
                                 nodeFile = 'nodes.txt',
                                 weighted = TRUE, threshold = 0.02)
```

## CEMiTool (Automated Analysis)

```r
library(CEMiTool)

# CEMiTool automates soft-threshold selection and module detection
expr_matrix <- read.csv('normalized_counts.csv', row.names = 1)
sample_annot <- read.csv('sample_annotation.csv')

cem <- cemitool(expr_matrix, sample_annot, filter = TRUE, plot = TRUE, verbose = TRUE)

# Results
nmodules(cem)
module_genes <- module_genes(cem)
generate_report(cem, directory = 'cemitool_report')

# Enrichment with gene sets (optional)
gene_sets <- read_gmt('pathways.gmt')
cem <- mod_ora(cem, gene_sets)
cem <- plot_ora(cem)
```

## hdWGCNA (Single-Cell)

```r
library(hdWGCNA)
library(Seurat)

seurat_obj <- readRDS('clustered.rds')
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = 'fraction',
                             fraction = 0.05, wgcna_name = 'hdwgcna')

# Metacells to reduce sparsity (groups of 25-50 cells)
seurat_obj <- MetacellsByGroups(seurat_obj, group.by = c('seurat_clusters'),
                                 k = 25, max_shared = 10, ident.group = 'seurat_clusters')
seurat_obj <- NormalizeMetacells(seurat_obj)

# Standard WGCNA on metacells
seurat_obj <- SetDatExpr(seurat_obj, group.by = 'seurat_clusters', group_name = 'all')
seurat_obj <- TestSoftPowers(seurat_obj)
seurat_obj <- ConstructNetwork(seurat_obj, soft_power = 6, setDatExpr = FALSE)
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

# Visualize on UMAP
seurat_obj <- RunModuleUMAP(seurat_obj, n_hubs = 10, n_neighbors = 15, min_dist = 0.1)
ModuleFeaturePlot(seurat_obj, features = 'hMEs')
```

## PyWGCNA (Python Alternative)

```python
import PyWGCNA

pywgcna = PyWGCNA.WGCNA(name='analysis', species='homo sapiens',
                         geneExp='normalized_counts.csv', outputPath='pywgcna_output/')
pywgcna.preprocess()
pywgcna.findModules()
pywgcna.updateSampleInfo(path='sample_traits.csv')
pywgcna.analyseWGCNA()
```

## Statistical Considerations

| Consideration | Guideline | Rationale |
|---------------|-----------|-----------|
| Minimum samples | 20 recommended, 15 absolute floor | Correlation stability; below 15, module detection is unreliable |
| Scale-free R^2 | > 0.85 default (0.80-0.90 acceptable) | Ensures biological network topology |
| Min module size | 30 genes | Smaller modules are often noise |
| Merge cut height | 0.25 | Merges modules with >75% eigengene correlation |
| Gene filtering | Top 5000-10000 most variable | Reduces noise, speeds computation |

## Related Skills

- scenic-regulons - TF-centric regulon inference from scRNA-seq
- differential-networks - Compare networks between conditions
- differential-expression/deseq2-basics - DE analysis to prioritize network genes
- differential-expression/batch-correction - Remove batch effects before network construction
- temporal-genomics/temporal-grn - Dynamic GRN inference from bulk time-series data
