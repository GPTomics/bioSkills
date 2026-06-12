# Reference: clusterProfiler 4.18.4+, org.Hs.eg.db 3.22+ | Verify API if version differs
# GO ORA across all three ontologies, simplified PER ONTOLOGY.
# simplify() operates on one ontology (GOSemSim similarity is defined within a single DAG),
# so an ont='ALL' object must be split into BP/MF/CC and each simplified separately.

library(clusterProfiler)
library(org.Hs.eg.db)

simplify_cut <- 0.7   # semantic-similarity redundancy cutoff; lower keeps more terms

de_results <- read.csv('de_results.csv')

sig_genes  <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]
all_tested <- de_results$gene_id[!is.na(de_results$pvalue)]   # universe = tested genes, NOT all rows, NOT the genome

fg_map <- bitr(sig_genes,  fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
bg_map <- bitr(all_tested, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list    <- unique(fg_map$ENTREZID)   # deduplicate one-to-many maps before counting
universe_ids <- unique(bg_map$ENTREZID)

simplified <- lapply(c('BP', 'MF', 'CC'), function(ont) {
    ego <- enrichGO(gene = gene_list, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID',
                    ont = ont, pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego <- simplify(ego, cutoff = simplify_cut, by = 'p.adjust', select_fun = min, measure = 'Wang')
    df <- as.data.frame(ego)
    cat(ont, ':', nrow(df), 'terms after simplify\n')
    df
})

combined <- do.call(rbind, simplified)
write.csv(combined, file.path(tempdir(), 'go_all_simplified.csv'), row.names = FALSE)
