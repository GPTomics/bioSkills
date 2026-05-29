# Reference: DiffCorr 0.4.1+, igraph 1.5+ | Verify API if version differs
# Differential co-expression analysis with DiffCorr

library(DiffCorr)
library(igraph)

expr_all <- read.csv('normalized_counts.csv', row.names = 1)
sample_info <- read.csv('sample_info.csv', row.names = 1)

# comp.2.cc.fdr expects genes in ROWS, samples in COLUMNS (it correlates rows = genes).
expr_ctrl <- expr_all[, sample_info$condition == 'control']
expr_disease <- expr_all[, sample_info$condition == 'disease']
cat('Control:', ncol(expr_ctrl), 'samples. Disease:', ncol(expr_disease), 'samples.\n')

# Top 3000 variable genes to reduce multiple testing burden
gene_vars <- apply(expr_all, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
expr_ctrl <- expr_ctrl[top_genes, ]
expr_disease <- expr_disease[top_genes, ]

# comp.2.cc.fdr returns the data.frame directly (save=TRUE also writes the file).
result <- comp.2.cc.fdr(
    output.file = 'diffcorr_results.txt',
    data1 = expr_ctrl,
    data2 = expr_disease,
    threshold = 0.05,
    save = TRUE
)

# Rename the space-containing output columns to safe names before use.
diffcorr <- result
names(diffcorr) <- c('gene1', 'gene2', 'r1', 'p1', 'r2', 'p2',
                     'p_diff', 'r_diff', 'lfdr1', 'lfdr2', 'lfdr_diff')
cat('Exported significant pairs:', nrow(diffcorr), '\n')

# Classify differential correlations
# threshold 0.3: minimum absolute correlation to consider an edge present
classify_edge <- function(cor1, cor2, threshold = 0.3) {
    if (abs(cor1) < threshold & abs(cor2) >= threshold) return('gained')
    if (abs(cor1) >= threshold & abs(cor2) < threshold) return('lost')
    if (cor1 > threshold & cor2 < -threshold) return('reversed')
    if (cor1 < -threshold & cor2 > threshold) return('reversed')
    return('unchanged')
}

diffcorr$edge_type <- mapply(classify_edge, diffcorr$r1, diffcorr$r2)
# Rows are already lfdr-thresholded by comp.2.cc.fdr; tighten on lfdr_diff if desired.
significant <- diffcorr[diffcorr$lfdr_diff < 0.05, ]
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
