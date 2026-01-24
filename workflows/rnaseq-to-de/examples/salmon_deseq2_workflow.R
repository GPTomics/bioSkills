# Complete RNA-seq workflow: Salmon + tximport + DESeq2

library(tximport)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)

# Public RNA-seq datasets:
# - Bioconductor: airway package (GSE52778, dexamethasone treatment)
# - recount3: https://rna.recount.bio (precomputed counts for GEO/SRA)
# - ENCODE: https://www.encodeproject.org (cell line RNA-seq)
# - GTEx: https://gtexportal.org (tissue expression)
# - GEO: GSE81089, GSE130963 (treatment studies with Salmon-compatible FASTQ)

# Configuration
samples <- c('control_1', 'control_2', 'control_3', 'treated_1', 'treated_2', 'treated_3')
conditions <- factor(c(rep('control', 3), rep('treated', 3)))
quant_dir <- 'salmon_quants'

# Create tx2gene mapping (example for Ensembl)
# In practice, generate this from your GTF:
# awk -F'\t' '$3=="transcript" {match($9, /gene_id "([^"]+)".*transcript_id "([^"]+)"/, a); print a[2]"\t"a[1]}' genes.gtf > tx2gene.tsv
tx2gene <- read.table('tx2gene.tsv', header = FALSE, col.names = c('TXNAME', 'GENEID'))

# Import Salmon quantifications
files <- file.path(quant_dir, samples, 'quant.sf')
names(files) <- samples
stopifnot(all(file.exists(files)))

# ignoreTxVersion=TRUE strips version suffixes (.1, .2) from Ensembl IDs to match tx2gene
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Create sample metadata
coldata <- data.frame(condition = conditions, row.names = samples)

# Create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat('Genes after filtering:', nrow(dds), '\n')

# Set reference level
dds$condition <- relevel(dds$condition, ref = 'control')

# Run DESeq2
dds <- DESeq(dds)

# Check available coefficients
resultsNames(dds)

# Get shrunken results
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
summary(res)

# QC: PCA plot
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = 'condition', returnData = TRUE)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = 'PCA of samples')
ggsave('pca_plot.pdf', width = 6, height = 5)

# QC: Sample distances
sample_dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dists), main = 'Sample distances')

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c('grey60', 'red3')) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey40') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'grey40') +
    theme_minimal() +
    labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 P-value')
ggsave('volcano_plot.pdf', width = 7, height = 6)

# Heatmap of top DE genes
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
if (length(sig_genes) > 0) {
    top_n <- min(50, length(sig_genes))
    top_genes <- head(sig_genes[order(res[sig_genes, 'padj'])], top_n)
    pheatmap(assay(vsd)[top_genes,], scale = 'row', show_rownames = (top_n <= 30),
             main = paste('Top', top_n, 'DE genes'))
}

# Export results
write.csv(as.data.frame(res), 'deseq2_all_results.csv')
write.csv(as.data.frame(res[sig_genes,]), 'deseq2_significant_genes.csv')

# Summary statistics
cat('\n=== Summary ===\n')
cat('Total genes tested:', nrow(res), '\n')
cat('Significant (padj < 0.05, |LFC| > 1):', length(sig_genes), '\n')
cat('Upregulated:', sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE), '\n')
cat('Downregulated:', sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE), '\n')
