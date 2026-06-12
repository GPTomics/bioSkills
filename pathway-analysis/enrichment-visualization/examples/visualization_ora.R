# Reference: enrichplot 1.30+, clusterProfiler 4.18+, ggplot2 3.5+ | Verify API if version differs
# Visualize an ORA enrichResult, choosing whether to SHOW or DELETE gene-set redundancy.
# Offline-runnable: needs only org.Hs.eg.db locally (no live KEGG/WikiPathways calls).

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# A small fixed Entrez gene list standing in for a DE hit list; replace with real input.
entrez <- c('1029', '1019', '1021', '890', '983', '991', '993', '4085', '4174', '5111',
            '1869', '1871', '7027', '4171', '4172', '4175', '4176', '4998', '5424', '5425')

ego <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = 'BP',
                pvalueCutoff = 0.05,   # filters on p.adjust by default; standard FDR gate
                qvalueCutoff = 0.2,    # clusterProfiler default secondary q-value gate
                readable = TRUE)       # map Entrez to symbols for plot labels

show_n <- 20   # dotplot top-N window; report the total significant count in the caption

# Collapse GO DAG redundancy BEFORE a flat dotplot (simplify lives conceptually in go-enrichment).
ego_simple <- simplify(ego, cutoff = 0.7, by = 'p.adjust', select_fun = min)   # 0.7 semantic-similarity cutoff; lower keeps more

p_dot <- dotplot(ego_simple, showCategory = show_n) + ggtitle('GO BP (redundancy collapsed)')
p_fold <- dotplot(ego_simple, x = 'FoldEnrichment', showCategory = show_n) + ggtitle('GO BP ordered by fold enrichment')

# SHOW the redundancy as structure: term-similarity matrix is mandatory before emapplot/treeplot.
ego_ts <- pairwise_termsim(ego)        # JC (Jaccard on gene overlap), the default; any gene-set type
p_emap <- emapplot(ego_ts, showCategory = 30)                   # edges = overlap >= min_edge (0.2 default)
p_tree <- treeplot(ego_ts, showCategory = 20, nCluster = 5)     # deterministic Ward clusters

out <- file.path(tempdir(), 'ora_visualization.pdf')
pdf(out, width = 11, height = 9)
print(p_dot); print(p_fold); print(p_emap); print(p_tree)
dev.off()

cat('Wrote', out, '\n')
