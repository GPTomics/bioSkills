# Basic DESeq2 workflow for differential expression analysis

library(DESeq2)
library(apeglm)

# --- ALTERNATIVE: Use real Bioconductor datasets ---
# The 'pasilla' or 'airway' packages provide realistic RNA-seq count data:
#
# library(pasilla)
# data('pasillaGenes')
# counts <- counts(pasillaGenes)
# coldata <- pData(pasillaGenes)[, c('condition', 'type')]
#
# Or with airway (Himes et al. 2014, dexamethasone-treated airway cells):
# library(airway)
# data('airway')
# counts <- assay(airway)
# coldata <- colData(airway)[, c('cell', 'dex')]

# Simulate example data with biologically realistic patterns
set.seed(42)
n_genes <- 1000
n_samples <- 6

base_counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                      nrow = n_genes,
                      dimnames = list(paste0('gene', 1:n_genes),
                                     paste0('sample', 1:n_samples)))

# Add differential expression signal: ~10% genes upregulated, ~10% downregulated
de_up <- 1:100
de_down <- 101:200
fold_change <- 2

counts <- base_counts
counts[de_up, 4:6] <- counts[de_up, 4:6] * fold_change
counts[de_down, 4:6] <- pmax(1, round(counts[de_down, 4:6] / fold_change))

coldata <- data.frame(
    condition = factor(rep(c('control', 'treated'), each = 3)),
    row.names = colnames(counts)
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat('Genes after filtering:', nrow(dds), '\n')

# Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Get results with shrinkage
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

# Summary
summary(res)

# Get significant genes
sig_genes <- subset(res, padj < 0.05)
cat('\nSignificant genes (padj < 0.05):', nrow(sig_genes), '\n')
