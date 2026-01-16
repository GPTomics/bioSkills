# Visualize GSEA results

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

de_results <- read.csv('de_results.csv')
gene_list <- de_results$log2FoldChange
names(gene_list) <- de_results$gene_id
gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)

gene_ids <- bitr(names(gene_list), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list_entrez <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID[match(names(gene_list_entrez), gene_ids$SYMBOL)]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

gse <- gseGO(geneList = gene_list_entrez, OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, verbose = FALSE)
gse <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

pdf('gsea_dotplot.pdf', width = 10, height = 8)
dotplot(gse, showCategory = 20, split = '.sign') + facet_grid(~.sign)
dev.off()

pdf('gsea_ridgeplot.pdf', width = 10, height = 10)
ridgeplot(gse, showCategory = 20) + theme(axis.text.y = element_text(size = 8))
dev.off()

results <- as.data.frame(gse)
if (nrow(results) > 0) {
    top_up <- which(results$NES > 0)[1]
    top_down <- which(results$NES < 0)[1]

    if (!is.na(top_up)) {
        pdf('gsea_plot_up.pdf', width = 8, height = 6)
        print(gseaplot2(gse, geneSetID = top_up, title = results$Description[top_up]))
        dev.off()
    }

    if (!is.na(top_down)) {
        pdf('gsea_plot_down.pdf', width = 8, height = 6)
        print(gseaplot2(gse, geneSetID = top_down, title = results$Description[top_down]))
        dev.off()
    }
}

cat('Saved GSEA visualization plots\n')
