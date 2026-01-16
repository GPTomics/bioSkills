# GSEA with Gene Ontology

library(clusterProfiler)
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')

gene_list <- -log10(de_results$pvalue) * sign(de_results$log2FoldChange)
names(gene_list) <- de_results$gene_id
gene_list <- gene_list[!is.na(gene_list) & is.finite(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

gene_ids <- bitr(names(gene_list), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

gene_list_entrez <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID[match(names(gene_list_entrez), gene_ids$SYMBOL)]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

cat('Ranked gene list:', length(gene_list_entrez), 'genes\n')

gse_bp <- gseGO(geneList = gene_list_entrez, OrgDb = org.Hs.eg.db, ont = 'BP', minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)

gse_bp <- setReadable(gse_bp, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

results <- as.data.frame(gse_bp)
cat('Enriched GO terms:', nrow(results), '\n')
cat('Positive NES (upregulated):', sum(results$NES > 0), '\n')
cat('Negative NES (downregulated):', sum(results$NES < 0), '\n')

write.csv(results, 'gsea_go_results.csv', row.names = FALSE)
print(head(gse_bp[, c('Description', 'NES', 'p.adjust')]))
