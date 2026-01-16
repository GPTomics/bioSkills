# KEGG module and pathway comparison

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

de_results <- read.csv('de_results.csv')

up_genes <- de_results %>% filter(padj < 0.05, log2FoldChange > 1) %>% pull(gene_id)
down_genes <- de_results %>% filter(padj < 0.05, log2FoldChange < -1) %>% pull(gene_id)

up_entrez <- bitr(up_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID
down_entrez <- bitr(down_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID

gene_lists <- list(upregulated = up_entrez, downregulated = down_entrez)

ck <- compareCluster(geneClusters = gene_lists, fun = 'enrichKEGG', organism = 'hsa', pvalueCutoff = 0.05)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

pdf('kegg_compare_clusters.pdf', width = 10, height = 8)
dotplot(ck, showCategory = 10) + ggtitle('KEGG Pathway Enrichment: Up vs Down')
dev.off()

all_genes <- c(up_entrez, down_entrez)
mkk <- enrichMKEGG(gene = all_genes, organism = 'hsa', pvalueCutoff = 0.05)

cat('KEGG Modules enriched:', nrow(as.data.frame(mkk)), '\n')
print(head(mkk))
