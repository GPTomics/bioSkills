# Time-Series DE Usage Guide

Identify genes with significant temporal expression patterns.

## Prerequisites

```r
BiocManager::install(c('limma', 'edgeR', 'splines', 'maSigPro'))
```

## Approach Comparison

| Method | Best For |
|--------|----------|
| limma + splines | Smooth temporal patterns |
| maSigPro | Multiple conditions over time |
| ImpulseDE2 | Impulse-like responses |
| DESeq2 LRT | Discrete time comparisons |

## Quick Start: limma with Splines

```r
library(limma)
library(edgeR)
library(splines)

# Load data
counts <- read.csv('counts.csv', row.names=1)
time <- c(0, 2, 4, 8, 12, 24)  # Timepoints
group <- factor(rep(c('ctrl', 'treat'), each=3))

# Create design with natural splines
design <- model.matrix(~ group * ns(time, df=3))

# voom transformation
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
v <- voom(dge, design)

# Fit model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Test for time effect
results <- topTable(fit, coef=grep('time', colnames(design)), number=Inf)
```

## maSigPro (Multi-condition Time Series)

```r
library(maSigPro)

# Design matrix
time <- rep(c(0, 2, 6, 12, 24), each=3)
group <- rep(c('ctrl', 'treat'), each=15)
replicates <- rep(1:3, 10)

edesign <- cbind(Time=time, Replicate=replicates,
                 ctrl=(group=='ctrl'), treat=(group=='treat'))

# Run maSigPro
design <- make.design.matrix(edesign, degree=2)
fit <- p.vector(counts, design, Q=0.05)
tstep <- T.fit(fit, step.method='backward')
sigs <- get.siggenes(tstep, rsq=0.6, vars='groups')
```

## DESeq2 Likelihood Ratio Test

```r
library(DESeq2)

# Full model with time
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ condition + time)

# Reduced model without time
dds <- DESeq(dds, test='LRT', reduced = ~ condition)
results <- results(dds)
```

## Clustering Temporal Patterns

```r
library(maSigPro)

# Cluster significant genes
see <- see.genes(sigs$sig.genes, show.fit=TRUE, k=6)

# Or use custom clustering
library(pheatmap)
sig_genes <- results[results$adj.P.Val < 0.05, ]
pheatmap(sig_genes, cluster_rows=TRUE, scale='row')
```

## Visualization

```r
library(ggplot2)

# Plot gene expression over time
plot_gene <- function(gene, counts, time, group) {
    df <- data.frame(
        time = time,
        expr = as.numeric(counts[gene, ]),
        group = group
    )
    ggplot(df, aes(x=time, y=expr, color=group)) +
        geom_point() +
        geom_smooth(method='loess') +
        labs(title=gene, y='Expression')
}
```

## Tips

- Choose spline degrees based on number of timepoints
- More timepoints = higher df possible
- Consider biological replicates at each timepoint
- Test for group:time interaction for condition-specific dynamics

## See Also

- [limma user guide](https://bioconductor.org/packages/limma/)
- [maSigPro vignette](https://bioconductor.org/packages/maSigPro/)
