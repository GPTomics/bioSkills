# Reference: WGCNA 1.72+ | Verify API if version differs
# Complete WGCNA co-expression network analysis (signed network, bicor, kME hubs)

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
# Pick the power on the SIGNED fit so it matches blockwiseModules below.
sft <- pickSoftThreshold(expr_data, powerVector = powers, networkType = 'signed', verbose = 5)

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
    # Signed network keeps activators and repressors in separate modules.
    networkType = 'signed', TOMType = 'signed',
    # bicor down-weights outliers; maxPOutliers caps false outlier flags at modest n.
    corType = 'bicor', maxPOutliers = 0.05,
    # minModuleSize 30: smaller modules are often noise
    minModuleSize = 30,
    # mergeCutHeight 0.25: merges modules with >75% eigengene correlation
    reassignThreshold = 0, mergeCutHeight = 0.25,
    # maxBlockSize >= n_genes keeps one block (no cross-block blindness artifact).
    maxBlockSize = ncol(expr_data) + 1,
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

# Recompute eigengenes from COLOR labels so ME/kME column names match module_colors
# (net$MEs is ME0/ME1... under numericLabels=TRUE and would not match color names).
MEs <- orderMEs(moduleEigengenes(expr_data, module_colors)$eigengenes)

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

# Hubs by module membership (kME): signed, bounded, comparable across modules.
module_of_interest <- names(which.min(module_trait_pval[, 1]))
module_color <- gsub('ME', '', module_of_interest)
module_genes <- colnames(expr_data)[module_colors == module_color]

kME <- signedKME(expr_data, MEs)
hub_ranking <- sort(kME[module_genes, paste0('kME', module_color)], decreasing = TRUE)
cat('Top hub genes in', module_color, 'module by kME:\n')
print(head(hub_ranking, 20))

save(net, MEs, module_colors, module_trait_cor, module_trait_pval,
     file = 'wgcna_results.RData')
cat('Saved WGCNA results\n')
