# Reference: clusterProfiler 4.10+, org.Hs.eg.db 3.18+, ReactomePA 1.46+, enrichplot 1.22+ | Verify API if version differs
# Orchestrate DE results -> ORA and GSEA -> redundancy-collapsed visualization.
# Decide the generation FIRST: a ranking for all genes -> GSEA; a pre-selected list -> ORA.
# NETWORK NOTE: enrichKEGG queries the LIVE KEGG REST API (needs internet, not reproducible across
# releases). enrichGO/gseGO and Reactome read local annotation and are reproducible given the release.

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

de_results_file <- 'deseq2_results.csv'
output_dir <- file.path(tempdir(), 'pathway_results')   # tempdir() so the repo stays clean
padj_cutoff <- 0.05
lfc_cutoff <- 1

dir.create(file.path(output_dir, 'plots'), showWarnings = FALSE, recursive = TRUE)

# === Stage 1: prepare the list, the ranked vector, and the universe ===
res <- read.csv(de_results_file, row.names = 1)

sig_genes <- rownames(subset(res, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff))
cat('Significant genes:', length(sig_genes), '\n')

# Background universe = genes that entered the DE test (testable), NOT the genome
universe_genes <- rownames(res[!is.na(res$pvalue), ])

# Ranked vector of ALL genes: prefer the Wald stat over a bare log2FC.
# NOTE: lfcShrink(apeglm/ashr) tables have NO `stat` column -- rank from the unshrunken
# results(dds)$stat (or limma `t` / edgeR sign(logFC)*-log10(PValue)).
ranked <- res$stat
names(ranked) <- rownames(res)
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)

# === Stage 2: convert IDs (ENTREZ is the safe downstream lingua franca) ===
sig_entrez <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
bg_entrez <- bitr(universe_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
conv_rate <- nrow(sig_entrez) / length(sig_genes)
cat('ID conversion rate:', round(conv_rate, 3), '\n')   # <0.85 flags a wrong ID type/organism

ranked_map <- bitr(names(ranked), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
ranked_list <- ranked[ranked_map$SYMBOL]
names(ranked_list) <- ranked_map$ENTREZID
ranked_list <- ranked_list[!duplicated(names(ranked_list))]   # dedup or GSEA biases the score
ranked_list <- sort(ranked_list, decreasing = TRUE)           # re-sort: bitr remap can reorder rows

# === Stage 3a: ORA branch (list + universe) ===
cat('\n=== GO ORA ===\n')
go_bp <- enrichGO(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)
go_bp <- simplify(go_bp, cutoff = 0.7, by = 'p.adjust')   # collapse DAG redundancy; one ontology at a time
cat('GO BP terms (simplified):', nrow(as.data.frame(go_bp)), '\n')

cat('\n=== KEGG ORA (LIVE DB) ===\n')
kegg <- enrichKEGG(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
cat('KEGG pathways:', nrow(as.data.frame(kegg)), '\n')

cat('\n=== Reactome ORA (local DB) ===\n')
reactome <- enrichPathway(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, organism = 'human', pvalueCutoff = 0.05, readable = TRUE)
cat('Reactome pathways:', nrow(as.data.frame(reactome)), '\n')

# === Stage 3b: GSEA branch (named decreasing vector of all genes) ===
cat('\n=== GSEA ===\n')
set.seed(123)   # permutation reproducibility; without it p-values drift across runs
gsea_go <- gseGO(ranked_list, OrgDb = org.Hs.eg.db, ont = 'BP',
                 minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
cat('GSEA GO terms:', nrow(as.data.frame(gsea_go)), '\n')

# === Stage 4: collapse redundancy, then visualize ===
cat('\n=== Visualization ===\n')
if (nrow(as.data.frame(go_bp)) > 0) {
    p <- dotplot(go_bp, showCategory = 20) + ggtitle('GO Biological Process Enrichment')
    ggsave(file.path(output_dir, 'plots', 'go_bp_dotplot.pdf'), p, width = 10, height = 8)
}
if (nrow(as.data.frame(go_bp)) > 5) {
    go_bp_sim <- pairwise_termsim(go_bp)   # required before emapplot
    p <- emapplot(go_bp_sim, showCategory = 30) + ggtitle('GO Term Similarity Network')
    ggsave(file.path(output_dir, 'plots', 'go_network.pdf'), p, width = 10, height = 10)
}
if (nrow(as.data.frame(gsea_go)) > 0) {
    p <- gseaplot2(gsea_go, geneSetID = 1:min(3, nrow(as.data.frame(gsea_go))), pvalue_table = TRUE)
    ggsave(file.path(output_dir, 'plots', 'gsea_plot.pdf'), p, width = 10, height = 8)
}

# === Export (provenance matters: record tool + DB version/date, metric, universe) ===
write.csv(as.data.frame(go_bp), file.path(output_dir, 'go_bp_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(kegg), file.path(output_dir, 'kegg_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(reactome), file.path(output_dir, 'reactome_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(gsea_go), file.path(output_dir, 'gsea_go_results.csv'), row.names = FALSE)

cat('\nResults written under:', output_dir, '\n')
