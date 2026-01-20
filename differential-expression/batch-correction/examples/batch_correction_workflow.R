library(DESeq2)
library(sva)
library(limma)
library(ggplot2)
library(patchwork)

counts <- read.csv('counts.csv', row.names = 1)
metadata <- read.csv('metadata.csv', row.names = 1)

stopifnot(all(colnames(counts) == rownames(metadata)))
stopifnot('batch' %in% colnames(metadata))
stopifnot('condition' %in% colnames(metadata))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
norm_expr <- assay(vsd)

pca <- prcomp(t(norm_expr), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     batch = metadata$batch, condition = metadata$condition)

p_before <- ggplot(pca_df, aes(PC1, PC2, color = batch, shape = condition)) +
    geom_point(size = 3) + theme_bw() + ggtitle('Before Batch Correction')

design(dds) <- ~ batch + condition
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'treatment', 'control'))
write.csv(as.data.frame(res), 'de_results_batch_adjusted.csv')

corrected_counts <- ComBat_seq(counts = as.matrix(counts),
                                batch = metadata$batch,
                                group = metadata$condition)
write.csv(corrected_counts, 'combat_seq_corrected_counts.csv')

design_matrix <- model.matrix(~ condition, data = metadata)
corrected_expr <- removeBatchEffect(norm_expr,
                                     batch = metadata$batch,
                                     design = design_matrix)

pca_after <- prcomp(t(corrected_expr), scale. = TRUE)
pca_df_after <- data.frame(PC1 = pca_after$x[, 1], PC2 = pca_after$x[, 2],
                           batch = metadata$batch, condition = metadata$condition)

p_after <- ggplot(pca_df_after, aes(PC1, PC2, color = batch, shape = condition)) +
    geom_point(size = 3) + theme_bw() + ggtitle('After Batch Correction')

pdf('batch_correction_pca.pdf', width = 12, height = 5)
print(p_before + p_after)
dev.off()

batch_cor_pc1 <- cor(pca$x[, 1], as.numeric(as.factor(metadata$batch)))
batch_cor_pc1_after <- cor(pca_after$x[, 1], as.numeric(as.factor(metadata$batch)))

cat('Batch correlation with PC1 (before):', round(batch_cor_pc1, 3), '\n')
cat('Batch correlation with PC1 (after):', round(batch_cor_pc1_after, 3), '\n')
