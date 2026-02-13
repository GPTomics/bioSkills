# Reference: WGCNA 1.72+ | Verify API if version differs
# Complete WGCNA co-expression network analysis

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

expr_data <- read.csv('normalized_counts.csv', row.names = 1)
expr_data <- t(expr_data)

gene_vars <- apply(expr_data, 2, var)
# Top 5000 most variable genes to reduce noise and speed computation
expr_data <- expr_data[, order(gene_vars, decreasing = TRUE)[1:5000]]
cat('Expression matrix:', nrow(expr_data), 'samples x', ncol(expr_data), 'genes\n')

sample_tree <- hclust(dist(expr_data), method = 'average')
pdf('sample_dendrogram.pdf', width = 12, height = 6)
plot(sample_tree, main = 'Sample clustering to detect outliers')
dev.off()

powers <- c(1:20)
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5)

pdf('soft_threshold.pdf', width = 10, height = 5)
par(mfrow = c(1, 2))
# R^2 > 0.85: default target for scale-free fit (acceptable range 0.8-0.9)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology R^2',
     main = 'Scale independence')
abline(h = 0.85, col = 'red')
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity',
     main = 'Mean connectivity')
dev.off()

soft_power <- sft$powerEstimate
cat('Selected soft power:', soft_power, '\n')

net <- blockwiseModules(
    expr_data, power = soft_power,
    TOMType = 'unsigned',
    # minModuleSize 30: smaller modules are often noise
    minModuleSize = 30,
    # mergeCutHeight 0.25: merges modules with >75% eigengene correlation
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE, saveTOMFileBase = 'TOM',
    verbose = 3
)

module_colors <- labels2colors(net$colors)
cat('Modules found:', length(unique(module_colors)), '\n')
print(table(module_colors))

pdf('module_dendrogram.pdf', width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                    'Module colors', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs <- net$MEs
MEs <- orderMEs(MEs)

traits <- read.csv('sample_traits.csv', row.names = 1)
module_trait_cor <- cor(MEs, traits, use = 'p')
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_data))

pdf('module_trait_heatmap.pdf', width = 10, height = 8)
textMatrix <- paste(signif(module_trait_cor, 2), '\n(',
                    signif(module_trait_pval, 1), ')', sep = '')
dim(textMatrix) <- dim(module_trait_cor)
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(traits), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.5,
               main = 'Module-trait relationships')
dev.off()

# Hub genes: high module membership AND high trait significance
module_of_interest <- names(which.min(module_trait_pval[, 1]))
module_color <- gsub('ME', '', module_of_interest)
module_genes <- colnames(expr_data)[module_colors == module_color]

gene_mm <- cor(expr_data, MEs, use = 'p')
gene_gs <- cor(expr_data, traits[, 1], use = 'p')

hub_genes <- module_genes[
    abs(gene_mm[module_genes, module_of_interest]) > 0.8 &
    abs(gene_gs[module_genes, 1]) > 0.2
]
cat('Hub genes in', module_color, 'module:', length(hub_genes), '\n')
if (length(hub_genes) > 0) print(head(hub_genes, 20))

connectivity <- intramodularConnectivity(adjacency(expr_data, power = soft_power), module_colors)
top_hubs <- connectivity[module_genes, ]
top_hubs <- top_hubs[order(top_hubs$kWithin, decreasing = TRUE), ]
cat('\nTop 10 hub genes by intramodular connectivity:\n')
print(head(top_hubs, 10))

save(net, MEs, module_colors, module_trait_cor, module_trait_pval,
     file = 'wgcna_results.RData')
cat('Saved WGCNA results\n')
