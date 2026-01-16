library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1, 'gene_symbol']

gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID

wp_result <- enrichWP(
    gene = entrez_ids,
    organism = 'Homo sapiens',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)

wp_readable <- setReadable(wp_result, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
results_df <- as.data.frame(wp_readable)
results_df

dotplot(wp_readable, showCategory = 15, title = 'WikiPathways Enrichment')

wp_result <- pairwise_termsim(wp_result)
emapplot(wp_result)

write.csv(results_df, 'wikipathways_enrichment_results.csv', row.names = FALSE)
