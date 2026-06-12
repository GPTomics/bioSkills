# Reference: ReactomePA 1.54+, clusterProfiler 4.18+ | Verify API if version differs
# reactome.db is a LOCAL Bioconductor annotation, so this script runs offline (no live database).
# enrichPathway has NO keyType argument: gene MUST be a character vector of ENTREZ ids -> bitr first.
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

sig_symbols <- c('CDK1','CCNB1','CCNB2','CDC20','BUB1','MAD2L1','PLK1','AURKA','AURKB','CDC25C',
                 'CCNA2','CDK2','E2F1','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','ORC1')
all_symbols <- keys(org.Hs.eg.db, keytype = 'SYMBOL')   # stand-in for the genes actually measured

sig_entrez <- bitr(sig_symbols, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID
universe   <- bitr(all_symbols, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID

pvalue_cutoff <- 0.05   # filters on p.adjust (BH) by default; standard FDR gate
qvalue_cutoff <- 0.2    # enrichPathway default secondary q-value gate
min_gs_size   <- 10     # drop tiny leaf pathways (2-3 genes) that inflate false positives
max_gs_size   <- 500    # drop huge top-level pathways that always enrich and are uninformative

ora <- enrichPathway(gene = sig_entrez, organism = 'human', universe = universe,
                     pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
                     minGSSize = min_gs_size, maxGSSize = max_gs_size, readable = TRUE)

results_df <- as.data.frame(ora)
results_df[, c('ID', 'Description', 'GeneRatio', 'BgRatio', 'FoldEnrichment', 'p.adjust', 'Count')]

# The top hits are often parent/child of one signal (e.g. cell-cycle checkpoints): report the
# deepest significant node and treat ancestors as context. ReactomePA has no simplify() equivalent.

# viewPathway takes the pathway NAME (Description), NOT the R-HSA id, and draws a LOCAL ggraph plot.
# Route it to a tempdir device so the default Rplots.pdf is never written to the working directory.
top_name <- results_df$Description[1]
pdf(file.path(tempdir(), 'reactome_viewpathway.pdf'))
print(viewPathway(top_name, organism = 'human', readable = TRUE))
dev.off()

# The interactive web diagram needs the R-HSA id (not viewPathway):
web_url <- paste0('https://reactome.org/PathwayBrowser/#/', results_df$ID[1])
web_url

write.csv(results_df, file.path(tempdir(), 'reactome_ora_results.csv'), row.names = FALSE)
