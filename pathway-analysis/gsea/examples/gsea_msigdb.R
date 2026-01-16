# GSEA with MSigDB Hallmark gene sets

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

de_results <- read.csv('de_results.csv')

gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$gene_id
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

gene_ids <- bitr(names(gene_list), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list_entrez <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID[match(names(gene_list_entrez), gene_ids$SYMBOL)]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

hallmarks <- msigdbr(species = 'Homo sapiens', category = 'H')
hallmarks_t2g <- hallmarks[, c('gs_name', 'entrez_gene')]

gse_hallmark <- GSEA(geneList = gene_list_entrez, TERM2GENE = hallmarks_t2g, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)

results <- as.data.frame(gse_hallmark)
cat('Enriched Hallmarks:', nrow(results), '\n')

up_hallmarks <- results[results$NES > 0, c('ID', 'NES', 'p.adjust')]
down_hallmarks <- results[results$NES < 0, c('ID', 'NES', 'p.adjust')]

cat('\nUpregulated hallmarks:\n')
print(head(up_hallmarks[order(-up_hallmarks$NES), ], 5))

cat('\nDownregulated hallmarks:\n')
print(head(down_hallmarks[order(down_hallmarks$NES), ], 5))

write.csv(results, 'gsea_hallmark_results.csv', row.names = FALSE)
