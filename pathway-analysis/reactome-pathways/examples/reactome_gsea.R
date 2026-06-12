# Reference: ReactomePA 1.54+, clusterProfiler 4.18+ | Verify API if version differs
# reactome.db is LOCAL, so this runs offline. gsePathway needs a NAMED numeric vector with ENTREZ
# names, sorted DECREASING; there is no keyType argument, so map non-Entrez ids first.
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# stand-in ranked statistic: a per-gene value (t-stat / signed -log10 p / shrunken log2FC) for every gene.
de <- data.frame(symbol = keys(org.Hs.eg.db, keytype = 'SYMBOL'), stringsAsFactors = FALSE)
set.seed(123)   # only to make this self-contained demo deterministic
de$stat <- rnorm(nrow(de))

mapped <- bitr(de$symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
de <- merge(de, mapped, by.x = 'symbol', by.y = 'SYMBOL')

gene_list <- de$stat
names(gene_list) <- de$ENTREZID                 # names MUST be ENTREZ
gene_list <- sort(gene_list, decreasing = TRUE)

pvalue_cutoff <- 0.05   # filters on p.adjust (BH) by default

set.seed(123)           # gsePathway permutes; fix the seed so p-values reproduce across runs
gse <- gsePathway(geneList = gene_list, organism = 'human',
                  pvalueCutoff = pvalue_cutoff, pAdjustMethod = 'BH', verbose = FALSE)

results_df <- as.data.frame(gse)
results_df

# gseaResult plots are owned by enrichment-visualization; shown here for a Reactome GSEA result.
# Route to a tempdir device so nothing lands in the working directory (ridgeplot needs ggridges).
if (nrow(results_df) > 0) {
    pdf(file.path(tempdir(), 'reactome_gsea_plots.pdf'))
    print(gseaplot2(gse, geneSetID = 1:min(3, nrow(results_df)), title = 'Reactome GSEA'))
    print(ridgeplot(gse, showCategory = 15))
    dev.off()
}

write.csv(results_df, file.path(tempdir(), 'reactome_gsea_results.csv'), row.names = FALSE)
