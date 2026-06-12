# Reference: clusterProfiler 4.18.4+ | Verify API if version differs
# KEGG over-representation (enrichKEGG) + module enrichment + reproducible gson pinning.
# NOTE: enrichKEGG/enrichMKEGG/gson_KEGG query the LIVE KEGG REST API (https://rest.kegg.jp/).
#       They need internet and are NOT reproducible across KEGG releases unless pinned.
#       Pin with gson_KEGG() and record the access date (below).

library(clusterProfiler)   # supplies gson_KEGG() and the enrich/GSEA functions
library(org.Hs.eg.db)
library(gson)              # supplies the GSON class + write.gson/read.gson (NOT gson_KEGG)

pvalue_cutoff <- 0.05   # filters on p.adjust by default in enrichResult; standard FDR gate
qvalue_cutoff <- 0.2    # clusterProfiler default; secondary q-value gate
min_gs_size   <- 10     # drop tiny sets that overfit; enrichKEGG default
max_gs_size   <- 500    # drop overly broad sets that always enrich; enrichKEGG default
lfc_cut       <- 1      # |log2FC| gate when selecting DE genes

de <- read.csv('de_results.csv')

# KEGG needs Entrez (keyType='ncbi-geneid'); ENSEMBL/SYMBOL into enrichKEGG return zero hits
sig_symbols <- de$gene[de$padj < pvalue_cutoff & abs(de$log2FoldChange) > lfc_cut]
sig_entrez  <- bitr(sig_symbols, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID

# universe = genes that could have been called DE (non-NA test statistic), same ID type
universe <- bitr(de$gene[!is.na(de$pvalue)], fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID
cat('Converted', length(sig_entrez), 'DE genes and', length(universe), 'universe genes to Entrez\n')

kk <- enrichKEGG(gene = sig_entrez, organism = 'hsa', keyType = 'ncbi-geneid', universe = universe,
                 pvalueCutoff = pvalue_cutoff, pAdjustMethod = 'BH',
                 minGSSize = min_gs_size, maxGSSize = max_gs_size, qvalueCutoff = qvalue_cutoff)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')   # eukaryotes only; no OrgDb -> keep raw IDs
cat('Found', nrow(as.data.frame(kk)), 'enriched KEGG pathways\n')

# KEGG modules (M-numbers): higher resolution, lower power, sparser coverage
mkk <- enrichMKEGG(gene = sig_entrez, organism = 'hsa', keyType = 'ncbi-geneid', universe = universe,
                   pvalueCutoff = pvalue_cutoff)
cat('Found', nrow(as.data.frame(mkk)), 'enriched KEGG modules\n')

# Reproducible pinning: snapshot the current KEGG release and run against it offline.
# use_internal_data=TRUE does NOT pin current KEGG (it loads the deprecated 2012 KEGG.db) -- use gson instead.
k <- gson_KEGG('hsa')
k@accessed_date <- as.character(Sys.Date())   # the GSON accessed_date slot survives write/read; a base attr() does not
snapshot <- file.path(tempdir(), 'kegg_hsa.gson')
write.gson(k, snapshot)
k <- read.gson(snapshot)

kk_pinned <- enricher(sig_entrez, gson = k, universe = universe,
                      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff)
cat('Pinned KEGG ORA against snapshot from', k@accessed_date, '-', nrow(as.data.frame(kk_pinned)), 'pathways\n')

write.csv(as.data.frame(kk), file.path(tempdir(), 'kegg_enrichment.csv'), row.names = FALSE)
