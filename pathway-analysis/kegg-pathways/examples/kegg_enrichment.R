# KEGG pathway enrichment analysis

library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')

sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
cat('Converted', length(gene_list), 'genes to Entrez IDs\n')

kk <- enrichKEGG(gene = gene_list, organism = 'hsa', pvalueCutoff = 0.05, pAdjustMethod = 'BH')

kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
print(head(kk_readable))

results_df <- as.data.frame(kk_readable)
write.csv(results_df, 'kegg_enrichment.csv', row.names = FALSE)

cat('Found', nrow(results_df), 'enriched KEGG pathways\n')
