# Reference: SPIA 2.50+, graphite 1.56+, clusterProfiler 4.18.4+ | Verify API if version differs
# Signed-topology perturbation analysis (the third pathway-analysis generation).
# SPIA propagates DE fold-changes through KEGG's signed signaling wiring (KGML), combining
# over-representation (pNDE) with perturbation (pPERT) into a global pG.
# NOTE: SPIA and graphite query / depend on the live KEGG release (graphite ships versioned
#       topology data, SPIA's bundled hsaSPIA is an OLDER snapshot). Pin the release before
#       publishing. SPIA is SIGNALING-ONLY: it is undefined for metabolic maps (e.g. glycolysis).

library(SPIA)
library(graphite)
library(clusterProfiler)
library(org.Hs.eg.db)

n_boot   <- 2000   # SPIA default; bootstrap replicates for the pPERT null
padj_cut <- 0.05   # DESeq2 adjusted-p gate for selecting DE genes

de <- read.csv('de_results.csv')

# SPIA needs a NAMED vector of log2 fold-changes (DE genes only) plus the universe, in Entrez space.
# bitr drops unmapped symbols and can map many-to-one, so it must be MERGED back by symbol -- never
# assigned as names directly (that recycles and silently attaches the wrong fold-change to each Entrez).
sig <- de[de$padj < padj_cut, ]
sig_map <- bitr(sig$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
de_vec  <- setNames(sig$log2FoldChange[match(sig_map$SYMBOL, sig$gene)], sig_map$ENTREZID)
de_vec  <- de_vec[!duplicated(names(de_vec))]

# universe = all measured genes (same ID space); SPIA aborts if >1% of DE IDs are absent from it
universe <- bitr(de$gene[!is.na(de$pvalue)], fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID

# Direct SPIA against KEGG (organism code; signaling maps only)
res <- spia(de = de_vec, all = universe, organism = 'hsa', nB = n_boot, plots = FALSE)
# output cols: Name, ID, pSize, NDE, pNDE, tA, pPERT, pG, pGFdr, pGFWER, Status, KEGGLINK
cat('SPIA scored', nrow(res), 'pathways;', sum(res$pGFdr < 0.05), 'significant after FDR\n')

# graphite route: harmonizes node IDs, resolves complexes/families, removes compounds,
# and runs SPIA over the cleaned graphs (also works on Reactome topology)
db <- pathways('hsapiens', 'kegg')
db <- convertIdentifiers(db, 'ENTREZID')
spia_set <- file.path(tempdir(), 'kegg_hsa_spia')
prepareSPIA(db, spia_set)
gr <- runSPIA(de = de_vec, all = universe, spia_set)

write.csv(res, file.path(tempdir(), 'spia_results.csv'), row.names = FALSE)
