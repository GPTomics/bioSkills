library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1, 'gene_symbol']

gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID

pathway_result <- enrichPathway(
    gene = entrez_ids,
    organism = 'human',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    readable = TRUE
)

results_df <- as.data.frame(pathway_result)
results_df

dotplot(pathway_result, showCategory = 15, title = 'Reactome Pathway Enrichment')

pathway_result <- pairwise_termsim(pathway_result)
emapplot(pathway_result)

write.csv(results_df, 'reactome_enrichment_results.csv', row.names = FALSE)
