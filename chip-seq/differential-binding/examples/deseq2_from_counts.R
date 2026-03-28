# Reference: DESeq2 1.42+, apeglm 1.24+ | Verify API if version differs
library(DESeq2)
library(apeglm)

counts <- read.delim('counts.tsv', row.names = 1, check.names = FALSE)
coldata <- data.frame(
    condition = factor(c(rep('treated', 3), rep('control', 3))),
    row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# ChIP-seq peaks are already enriched regions; filter less aggressively than RNA-seq
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = 'control')
dds <- DESeq(dds)

# alpha matches intended significance threshold for optimal independent filtering
res <- results(dds, alpha = 0.05)
summary(res)

res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, 'TRUE', 'FALSE')

out <- res_df[, c('peak_id', 'log2FoldChange', 'pvalue', 'padj', 'significant')]
colnames(out)[colnames(out) == 'log2FoldChange'] <- 'log2fc'
write.table(out, file = 'differential.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

# LFC shrinkage for ranking and visualization (does not change padj)
resLFC <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
