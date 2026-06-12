# Reference: enrichplot 1.30+, clusterProfiler 4.18+, ggplot2 3.5+ | Verify API if version differs
# Visualize a GSEA gseaResult preserving DIRECTION (NES sign).
# Note: there is NO barplot method for gseaResult by design - a bar from zero cannot carry a signed NES.
# Note: ridgeplot() needs the ggridges package (enrichplot Suggests-only); install.packages('ggridges') if missing.
# Offline-runnable: needs only org.Hs.eg.db locally.

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

set.seed(123)   # fix the GSEA permutation seed so p-values are reproducible

# A named numeric vector sorted DECREASING by the ranking metric (here a synthetic log2FC).
# Replace with a real ranking from differential-expression/de-results.
genes <- keys(org.Hs.eg.db, keytype = 'ENTREZID')[1:3000]
gene_list <- sort(setNames(rnorm(length(genes), sd = 1.5), genes), decreasing = TRUE)

gse <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = 'BP',
             minGSSize = 10, maxGSSize = 500,   # drop tiny overfit and overly broad always-enriched sets
             pvalueCutoff = 0.25, verbose = FALSE)
gse <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

show_n <- 20

out <- file.path(tempdir(), 'gsea_visualization.pdf')
pdf(out, width = 10, height = 9)

if (nrow(as.data.frame(gse)) > 0) {
    print(dotplot(gse, x = 'NES', showCategory = show_n, color = 'p.adjust') + ggtitle('GSEA: signed NES'))
    print(ridgeplot(gse, showCategory = show_n) + theme(axis.text.y = element_text(size = 8)))
    print(gseaplot2(gse, geneSetID = 1, title = as.data.frame(gse)$Description[1]))
} else {
    cat('No enriched sets in this synthetic run; rerun with real ranked data.\n')
}

dev.off()
cat('Wrote', out, '\n')
