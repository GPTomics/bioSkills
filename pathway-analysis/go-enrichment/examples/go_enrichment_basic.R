# Basic GO enrichment analysis

library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')

sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID

ego_bp <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2)

ego_bp <- setReadable(ego_bp, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

print(head(ego_bp))

results_df <- as.data.frame(ego_bp)
write.csv(results_df, 'go_bp_enrichment.csv', row.names = FALSE)

cat('Found', nrow(results_df), 'enriched GO terms\n')
