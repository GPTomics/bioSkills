# Complete pathway enrichment analysis from DE results

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Configuration
de_results_file <- 'deseq2_results.csv'
output_dir <- 'pathway_results'
padj_cutoff <- 0.05
lfc_cutoff <- 1

dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, 'plots'), showWarnings = FALSE)

# === Load and Prepare Data ===
cat('Loading DE results...\n')
res <- read.csv(de_results_file, row.names = 1)

# Significant genes
sig_up <- rownames(subset(res, padj < padj_cutoff & log2FoldChange > lfc_cutoff))
sig_down <- rownames(subset(res, padj < padj_cutoff & log2FoldChange < -lfc_cutoff))
sig_all <- c(sig_up, sig_down)

cat('Upregulated genes:', length(sig_up), '\n')
cat('Downregulated genes:', length(sig_down), '\n')

# Convert to Entrez IDs
convert_ids <- function(genes) {
    converted <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
    return(converted$ENTREZID)
}

entrez_up <- convert_ids(sig_up)
entrez_down <- convert_ids(sig_down)
entrez_all <- convert_ids(sig_all)

# Ranked list for GSEA
ranked <- res$log2FoldChange
names(ranked) <- rownames(res)
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)
ranked_conv <- bitr(names(ranked), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
ranked_list <- ranked[ranked_conv$SYMBOL]
names(ranked_list) <- ranked_conv$ENTREZID

# === GO Enrichment ===
cat('\n=== GO Enrichment ===\n')

go_bp_all <- enrichGO(entrez_all, OrgDb = org.Hs.eg.db, ont = 'BP',
                      pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)
go_bp_simple <- simplify(go_bp_all, cutoff = 0.7)
cat('GO BP terms (simplified):', nrow(as.data.frame(go_bp_simple)), '\n')

go_mf <- enrichGO(entrez_all, OrgDb = org.Hs.eg.db, ont = 'MF',
                  pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)
cat('GO MF terms:', nrow(as.data.frame(go_mf)), '\n')

go_cc <- enrichGO(entrez_all, OrgDb = org.Hs.eg.db, ont = 'CC',
                  pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)
cat('GO CC terms:', nrow(as.data.frame(go_cc)), '\n')

# === KEGG Enrichment ===
cat('\n=== KEGG Enrichment ===\n')
kegg <- enrichKEGG(entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
cat('KEGG pathways:', nrow(as.data.frame(kegg)), '\n')

# === Reactome Enrichment ===
cat('\n=== Reactome Enrichment ===\n')
reactome <- enrichPathway(entrez_all, organism = 'human', pvalueCutoff = 0.05, readable = TRUE)
cat('Reactome pathways:', nrow(as.data.frame(reactome)), '\n')

# === GSEA ===
cat('\n=== GSEA ===\n')
gsea_go <- gseGO(ranked_list, OrgDb = org.Hs.eg.db, ont = 'BP',
                 minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
cat('GSEA GO terms:', nrow(as.data.frame(gsea_go)), '\n')

gsea_kegg <- gseKEGG(ranked_list, organism = 'hsa',
                     minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
cat('GSEA KEGG pathways:', nrow(as.data.frame(gsea_kegg)), '\n')

# === Visualizations ===
cat('\n=== Creating Visualizations ===\n')

# GO dot plot
if (nrow(as.data.frame(go_bp_simple)) > 0) {
    p <- dotplot(go_bp_simple, showCategory = 20) +
        ggtitle('GO Biological Process Enrichment')
    ggsave(file.path(output_dir, 'plots', 'go_bp_dotplot.pdf'), p, width = 10, height = 8)
}

# KEGG bar plot
if (nrow(as.data.frame(kegg)) > 0) {
    p <- barplot(kegg, showCategory = 15) + ggtitle('KEGG Pathway Enrichment')
    ggsave(file.path(output_dir, 'plots', 'kegg_barplot.pdf'), p, width = 9, height = 6)
}

# Reactome dot plot
if (nrow(as.data.frame(reactome)) > 0) {
    p <- dotplot(reactome, showCategory = 15) + ggtitle('Reactome Pathway Enrichment')
    ggsave(file.path(output_dir, 'plots', 'reactome_dotplot.pdf'), p, width = 10, height = 8)
}

# Enrichment map
if (nrow(as.data.frame(go_bp_simple)) > 5) {
    go_bp_sim <- pairwise_termsim(go_bp_simple)
    p <- emapplot(go_bp_sim, showCategory = 30) + ggtitle('GO Term Network')
    ggsave(file.path(output_dir, 'plots', 'go_network.pdf'), p, width = 10, height = 10)
}

# GSEA plot
if (nrow(as.data.frame(gsea_go)) > 0) {
    p <- gseaplot2(gsea_go, geneSetID = 1:min(3, nrow(as.data.frame(gsea_go))), pvalue_table = TRUE)
    ggsave(file.path(output_dir, 'plots', 'gsea_plot.pdf'), p, width = 10, height = 8)
}

# === Export Results ===
cat('\n=== Exporting Results ===\n')
write.csv(as.data.frame(go_bp_simple), file.path(output_dir, 'go_bp_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(go_mf), file.path(output_dir, 'go_mf_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(go_cc), file.path(output_dir, 'go_cc_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(kegg), file.path(output_dir, 'kegg_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(reactome), file.path(output_dir, 'reactome_enrichment.csv'), row.names = FALSE)
write.csv(as.data.frame(gsea_go), file.path(output_dir, 'gsea_go_results.csv'), row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), file.path(output_dir, 'gsea_kegg_results.csv'), row.names = FALSE)

# Summary table
cat('\n=== Analysis Complete ===\n')
cat('Results saved to:', output_dir, '\n')
cat('\nSummary:\n')
cat('  GO BP (simplified):', nrow(as.data.frame(go_bp_simple)), 'terms\n')
cat('  GO MF:', nrow(as.data.frame(go_mf)), 'terms\n')
cat('  GO CC:', nrow(as.data.frame(go_cc)), 'terms\n')
cat('  KEGG:', nrow(as.data.frame(kegg)), 'pathways\n')
cat('  Reactome:', nrow(as.data.frame(reactome)), 'pathways\n')
cat('  GSEA GO:', nrow(as.data.frame(gsea_go)), 'terms\n')
cat('  GSEA KEGG:', nrow(as.data.frame(gsea_kegg)), 'pathways\n')
