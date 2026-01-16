# GO enrichment for all three ontologies

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results %>% filter(padj < 0.05, abs(log2FoldChange) > 1) %>% pull(gene_id)

gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

all_genes <- de_results$gene_id
universe_ids <- bitr(all_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego_all <- enrichGO(gene = gene_ids$ENTREZID, universe = universe_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)

results <- as.data.frame(ego_all)
for (ont in c('BP', 'MF', 'CC')) {
    ont_results <- results[results$ONTOLOGY == ont, ]
    cat(ont, ':', nrow(ont_results), 'terms\n')
}

ego_simplified <- simplify(ego_all, cutoff = 0.7, by = 'p.adjust', select_fun = min)
cat('After simplification:', nrow(as.data.frame(ego_simplified)), 'terms\n')

write.csv(as.data.frame(ego_simplified), 'go_all_simplified.csv', row.names = FALSE)
