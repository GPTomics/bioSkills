# Reference: clusterProfiler 4.18+, rWikiPathways 1.26+ | Verify API if version differs
# NOTE: enrichWP/gseWP download the current/ WikiPathways GMT over the network at run time
# (needs internet). current/ is the latest MONTHLY release, so identical code returns different
# results across months. For a reproducible run, pin a dated GMT (see wikipathways_explore.R).
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

pvalue_cut <- 0.05   # filters on p.adjust (BH) by default; standard FDR gate
qvalue_cut <- 0.2    # clusterProfiler enricher default secondary q-value gate
min_gs     <- 10     # drop tiny WP pathways that overfit; enricher default
max_gs     <- 500    # drop overly broad sets that always enrich; enricher default

de_results <- read.csv('de_results.csv')
sig_symbols <- de_results[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1, 'gene_symbol']

# WP GMT is Entrez-keyed; symbols/Ensembl overlap nothing with no error -> convert first
sig <- bitr(sig_symbols, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID
all_entrez <- bitr(de_results$gene_symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID

# universe = tested genes; default NULL = all-WP-genes background inflates significance
wp <- enrichWP(gene = sig, organism = 'Homo sapiens', universe = all_entrez,
               pvalueCutoff = pvalue_cut, pAdjustMethod = 'BH',
               minGSSize = min_gs, maxGSSize = max_gs, qvalueCutoff = qvalue_cut)

wp <- setReadable(wp, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
results_df <- as.data.frame(wp)
results_df

# GSEA path: named Entrez vector sorted decreasing; no universe (FCS uses the whole list)
# WP GMT is Entrez-keyed, so map symbols and carry log2FC by match (never names(x)<-bitr()$ENTREZID)
gsea_map <- bitr(de_results$gene_symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gl <- setNames(de_results$log2FoldChange[match(gsea_map$SYMBOL, de_results$gene_symbol)], gsea_map$ENTREZID)
gl <- sort(gl[!duplicated(names(gl))], decreasing = TRUE)
set.seed(123)   # fix permutation reproducibility; p-values drift across runs otherwise
wp_gsea <- gseWP(geneList = gl, organism = 'Homo sapiens',
                 pvalueCutoff = pvalue_cut, pAdjustMethod = 'BH',
                 minGSSize = min_gs, maxGSSize = max_gs)
as.data.frame(wp_gsea)

# plots are owned by enrichment-visualization; route to a tempdir device so nothing lands in the CWD
wp <- pairwise_termsim(wp)   # required before emapplot
pdf(file.path(tempdir(), 'wikipathways_plots.pdf'))
print(dotplot(wp, showCategory = 15, title = 'WikiPathways Enrichment'))
print(emapplot(wp))
dev.off()

write.csv(results_df, file.path(tempdir(), 'wikipathways_enrichment_results.csv'), row.names = FALSE)
