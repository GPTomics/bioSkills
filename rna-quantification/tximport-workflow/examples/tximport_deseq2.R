library(tximport)
library(DESeq2)

sample_names <- c('ctrl1', 'ctrl2', 'ctrl3', 'treat1', 'treat2', 'treat3')
files <- file.path(paste0(sample_names, '_quant'), 'quant.sf')
names(files) <- sample_names

tx2gene <- read.csv('tx2gene.csv')

txi <- tximport(files, type = 'salmon', tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

coldata <- data.frame(
    condition = factor(rep(c('control', 'treated'), each = 3)),
    row.names = sample_names
)

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c('condition', 'treated', 'control'))
res_sig <- subset(res, padj < 0.05)

write.csv(as.data.frame(res), 'deseq2_results.csv')
write.csv(as.data.frame(res_sig), 'deseq2_significant.csv')

cat('Total genes tested:', nrow(res), '\n')
cat('Significant genes (padj < 0.05):', nrow(res_sig), '\n')
