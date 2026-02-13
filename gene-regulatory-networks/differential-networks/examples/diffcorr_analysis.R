# Reference: matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs
# Differential co-expression analysis with DiffCorr

library(DiffCorr)
library(igraph)

expr_all <- read.csv('normalized_counts.csv', row.names = 1)
sample_info <- read.csv('sample_info.csv', row.names = 1)

expr_ctrl <- t(expr_all[, sample_info$condition == 'control'])
expr_disease <- t(expr_all[, sample_info$condition == 'disease'])
cat('Control:', nrow(expr_ctrl), 'samples. Disease:', nrow(expr_disease), 'samples.\n')

# Top 3000 variable genes to reduce multiple testing burden
gene_vars <- apply(expr_all, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
expr_ctrl <- expr_ctrl[, top_genes]
expr_disease <- expr_disease[, top_genes]

result <- comp.2.cc.fdr(
    output.file = 'diffcorr_results.txt',
    data1 = expr_ctrl,
    data2 = expr_disease,
    threshold = 0.05
)

diffcorr <- read.delim('diffcorr_results.txt')
cat('Total tested pairs:', nrow(diffcorr), '\n')

# Classify differential correlations
# threshold 0.3: minimum absolute correlation to consider an edge present
classify_edge <- function(cor1, cor2, threshold = 0.3) {
    if (abs(cor1) < threshold & abs(cor2) >= threshold) return('gained')
    if (abs(cor1) >= threshold & abs(cor2) < threshold) return('lost')
    if (cor1 > threshold & cor2 < -threshold) return('reversed')
    if (cor1 < -threshold & cor2 > threshold) return('reversed')
    return('unchanged')
}

diffcorr$edge_type <- mapply(classify_edge, diffcorr$cor1, diffcorr$cor2)
significant <- diffcorr[diffcorr$p.adj < 0.05, ]
rewired <- significant[significant$edge_type != 'unchanged', ]

cat('\nEdge classification (significant only):\n')
print(table(significant$edge_type))
cat('\nTotal rewired edges:', nrow(rewired), '\n')

# Rewired hub genes (most differential connections)
gene_rewiring <- c(as.character(rewired$gene1), as.character(rewired$gene2))
rewiring_counts <- sort(table(gene_rewiring), decreasing = TRUE)
cat('\nTop 20 rewired hub genes:\n')
print(head(rewiring_counts, 20))

# Build differential network for visualization
edges <- rewired[, c('gene1', 'gene2', 'edge_type')]
g <- graph_from_data_frame(edges, directed = FALSE)

color_map <- c(gained = '#2ca02c', lost = '#d62728', reversed = '#9467bd')
E(g)$color <- color_map[E(g)$edge_type]

pdf('differential_network.pdf', width = 12, height = 12)
plot(g, vertex.size = 3, vertex.label.cex = 0.5, edge.width = 0.5,
     main = 'Differential co-expression network')
legend('topleft', legend = names(color_map), col = color_map, lwd = 2)
dev.off()

write.csv(rewired, 'rewired_edges.csv', row.names = FALSE)
cat('Saved rewired_edges.csv and differential_network.pdf\n')
