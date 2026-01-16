library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

de_results <- read.csv('de_results.csv')

gene_ids <- bitr(de_results$gene_symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
de_merged <- merge(de_results, gene_ids, by.x = 'gene_symbol', by.y = 'SYMBOL')

gene_list <- de_merged$log2FoldChange
names(gene_list) <- de_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_result <- gsePathway(
    geneList = gene_list,
    organism = 'human',
    pvalueCutoff = 0.1,
    pAdjustMethod = 'BH',
    verbose = FALSE
)

results_df <- as.data.frame(gsea_result)
results_df

gseaplot2(gsea_result, geneSetID = 1:3, title = 'Reactome GSEA')

ridgeplot(gsea_result, showCategory = 15)

write.csv(results_df, 'reactome_gsea_results.csv', row.names = FALSE)
