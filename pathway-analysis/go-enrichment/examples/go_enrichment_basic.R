# Reference: clusterProfiler 4.18.4+, org.Hs.eg.db 3.22+ | Verify API if version differs
# Basic GO over-representation analysis: foreground + matched universe, ont set explicitly.
# Self-contained: draws a small foreground and a tested-gene universe from org.Hs.eg.db so it runs offline.

library(clusterProfiler)
library(org.Hs.eg.db)

pvalue_cutoff <- 0.05   # filters p.adjust (despite the name); standard FDR gate
qvalue_cutoff <- 0.2    # clusterProfiler default secondary q-value gate
min_gs        <- 10     # drop tiny over-fitting gene sets; enrichGO default
max_gs        <- 500    # drop overly broad always-enriched sets; enrichGO default

all_entrez   <- keys(org.Hs.eg.db, keytype = 'ENTREZID')
universe_ids <- head(all_entrez, 3000)   # the genes TESTED in this assay, NOT the whole genome
gene_list    <- head(universe_ids, 200)  # foreground = the flagged hits

# enrichGO source default is ont='MF'; always set ont explicitly.
ego_bp <- enrichGO(gene = gene_list, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID',
                   ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
                   minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE)

results_df <- as.data.frame(ego_bp)
print(head(results_df))
cat('Found', nrow(results_df), 'enriched GO BP terms\n')

write.csv(results_df, file.path(tempdir(), 'go_bp_enrichment.csv'), row.names = FALSE)
