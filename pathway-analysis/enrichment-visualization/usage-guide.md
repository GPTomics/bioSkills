# Enrichment Visualization Usage Guide

The enrichplot package provides visualization functions for clusterProfiler results.

## Prerequisites

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'enrichplot'))
```

## Quick Reference

| Plot Type | Function | Description |
|-----------|----------|-------------|
| Dot plot | `dotplot()` | Points sized by count, colored by p-value |
| Bar plot | `barplot()` | Bars showing count or gene ratio |
| Network | `cnetplot()` | Gene-concept relationships |
| Map | `emapplot()` | Term similarity clusters |
| Tree | `treeplot()` | Hierarchical term grouping |
| Upset | `upsetplot()` | Overlapping genes |
| GSEA | `gseaplot2()` | Running enrichment score |
| Ridge | `ridgeplot()` | Fold change distribution |
| Heatmap | `heatplot()` | Gene-concept matrix |

## Choosing a Visualization

### For Over-Representation Results

**Starting point:** `dotplot()` - best overview
**Show relationships:** `cnetplot()` - genes and terms
**Cluster terms:** `emapplot()` - similar terms grouped
**Compare groups:** `dotplot()` with compareCluster result

### For GSEA Results

**Starting point:** `ridgeplot()` - all gene sets
**Detailed view:** `gseaplot2()` - running score for specific sets
**Overview:** `dotplot()` - works for GSEA too

## Dot Plot (Most Common)

```r
library(enrichplot)

dotplot(ego, showCategory = 20)
```

**What it shows:**
- X-axis: Gene ratio (proportion of your genes in term)
- Y-axis: Term description
- Dot size: Number of genes
- Dot color: Adjusted p-value

**Customization:**
```r
dotplot(ego, showCategory = 15, font.size = 10) +
    ggtitle('GO Enrichment') +
    theme(axis.text.y = element_text(size = 9))
```

## Gene-Concept Network

Shows which genes contribute to which terms:

```r
cnetplot(ego)

# Color by fold change
cnetplot(ego, foldChange = gene_list)

# Circular layout
cnetplot(ego, circular = TRUE)
```

## Enrichment Map

Clusters similar terms based on shared genes:

```r
# Must compute similarity first
ego_sim <- pairwise_termsim(ego)
emapplot(ego_sim)

# Group clusters
emapplot(ego_sim, group_category = TRUE)
```

## GSEA Plots

### Running Score Plot

```r
# By index
gseaplot2(gse, geneSetID = 1)

# Multiple gene sets
gseaplot2(gse, geneSetID = 1:3)

# By term name
gseaplot2(gse, geneSetID = 'HALLMARK_INFLAMMATORY_RESPONSE')
```

### Ridge Plot

Shows fold change distribution for each gene set:

```r
ridgeplot(gse, showCategory = 15)
```

## Saving Plots

### PDF (Best for Publication)

```r
pdf('figure.pdf', width = 10, height = 8)
dotplot(ego, showCategory = 20)
dev.off()
```

### PNG (Good for Presentations)

```r
png('figure.png', width = 1000, height = 800, res = 150)
dotplot(ego, showCategory = 20)
dev.off()
```

### Using ggsave

```r
p <- dotplot(ego)
ggsave('figure.pdf', p, width = 10, height = 8)
ggsave('figure.png', p, width = 10, height = 8, dpi = 150)
```

## Common Customizations

### Change Number of Terms

```r
dotplot(ego, showCategory = 30)  # Show more
dotplot(ego, showCategory = 10)  # Show fewer
```

### Adjust Text Size

```r
dotplot(ego) + theme(axis.text.y = element_text(size = 8))
```

### Change Colors

```r
dotplot(ego) + scale_color_viridis_c()
dotplot(ego) + scale_color_gradient(low = 'blue', high = 'red')
```

### Add Title

```r
dotplot(ego) + ggtitle('GO Biological Process Enrichment')
```

## Troubleshooting

**Terms truncated:**
- Increase plot width in pdf/png
- Use `theme(axis.text.y = element_text(size = 8))`

**emapplot error:**
- Run `pairwise_termsim()` first

**cnetplot too crowded:**
- Use `showCategory = 5` to limit terms
- Try `circular = TRUE` layout
