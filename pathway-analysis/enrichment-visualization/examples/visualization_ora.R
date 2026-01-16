# Visualize over-representation analysis results

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)

pdf('go_dotplot.pdf', width = 10, height = 8)
dotplot(ego, showCategory = 20) + ggtitle('GO Biological Process Enrichment')
dev.off()

pdf('go_barplot.pdf', width = 10, height = 6)
barplot(ego, showCategory = 15)
dev.off()

gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$gene_id

pdf('go_cnetplot.pdf', width = 12, height = 10)
cnetplot(ego, foldChange = gene_list, showCategory = 5, cex_label_gene = 0.7)
dev.off()

ego_sim <- pairwise_termsim(ego)

pdf('go_emapplot.pdf', width = 10, height = 10)
emapplot(ego_sim, showCategory = 30)
dev.off()

pdf('go_treeplot.pdf', width = 12, height = 8)
treeplot(ego_sim, showCategory = 20)
dev.off()

cat('Saved visualization plots\n')
