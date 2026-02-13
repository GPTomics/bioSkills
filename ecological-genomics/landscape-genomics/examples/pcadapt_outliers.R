# Reference: vegan 2.6+ | Verify API if version differs
library(pcadapt)
library(qvalue)

# --- Step 1: Read genotype data ---
# Supports bed/bim/fam (PLINK) or lfmm format
geno <- read.pcadapt('genotypes.bed', type = 'bed')

# --- Step 2: Choose K from screeplot ---
# Run with large K to see eigenvalue decay
# K=20: scan up to 20 PCs to identify the elbow
x_scree <- pcadapt(geno, K = 20)

pdf('pcadapt_screeplot.pdf', width = 7, height = 5)
plot(x_scree, option = 'screeplot', K = 20)
dev.off()

# Proportion of variance explained by each PC
pve <- (x_scree$singular.values^2) / sum(x_scree$singular.values^2)
cat('Proportion of variance explained by first 10 PCs:\n')
for (i in 1:min(10, length(pve))) {
    cat(sprintf('  PC%d: %.3f\n', i, pve[i]))
}

# --- Step 3: Run pcadapt with selected K ---
# K=3: number of principal components capturing population structure
# Determined from screeplot elbow above
best_K <- 3
x <- pcadapt(geno, K = best_K)

# --- Step 4: Diagnostic plots ---
pdf('pcadapt_diagnostics.pdf', width = 10, height = 8)
par(mfrow = c(2, 2))

# QQ-plot: points should follow diagonal except in the upper tail
# Inflation in the body indicates poor K choice or structure correction
plot(x, option = 'qqplot')

# Manhattan plot: -log10(p-values) across loci
plot(x, option = 'manhattan')

# Histogram of p-values: should be uniform with a spike near 0
hist(x$pvalues, breaks = 50, col = 'grey80', border = 'grey40',
     main = 'P-value Distribution', xlab = 'p-value')

# Score plot (PC1 vs PC2): population structure visualization
plot(x, option = 'scores', i = 1, j = 2)
dev.off()

# --- Step 5: Identify outlier loci ---
# q-value < 0.05: 5% false discovery rate (Storey method)
qvals <- qvalue(x$pvalues)$qvalues

# Count outliers at different thresholds
for (threshold in c(0.01, 0.05, 0.10)) {
    n_outliers <- sum(qvals < threshold, na.rm = TRUE)
    cat(sprintf('Outliers at q < %.2f: %d\n', threshold, n_outliers))
}

outliers <- which(qvals < 0.05)
cat('\nOutlier loci (q < 0.05):', length(outliers), '\n')

# --- Step 6: Examine PC loadings for outliers ---
# Loadings indicate which PCs drive each outlier
# Large loading on PC1 = differentiated along primary structure axis
loadings <- x$loadings
outlier_loadings <- loadings[outliers, ]

cat('\nPC contributions for top outliers:\n')
top_outliers <- outliers[order(qvals[outliers])][1:min(10, length(outliers))]
for (locus in top_outliers) {
    pc_loads <- abs(loadings[locus, ])
    dominant_pc <- which.max(pc_loads)
    cat(sprintf('  Locus %d: dominant PC%d (loading = %.4f), q = %.4e\n',
                locus, dominant_pc, loadings[locus, dominant_pc], qvals[locus]))
}

# --- Step 7: Publication Manhattan plot ---
pdf('pcadapt_manhattan_publication.pdf', width = 12, height = 5)
pvals <- x$pvalues
log_pvals <- -log10(pvals)
colors <- ifelse(qvals < 0.05, 'red', 'grey60')

plot(log_pvals, pch = 19, cex = 0.4, col = colors,
     xlab = 'SNP index', ylab = expression(-log[10](italic(p))),
     main = sprintf('pcadapt Outlier Detection (K=%d)', best_K))

# Bonferroni threshold (conservative reference line)
bonf_threshold <- -log10(0.05 / length(pvals))
abline(h = bonf_threshold, col = 'blue', lty = 2)
legend('topright', legend = c(sprintf('Outliers (q<0.05, n=%d)', length(outliers)),
                               'Bonferroni threshold'),
       col = c('red', 'blue'), pch = c(19, NA), lty = c(NA, 2), cex = 0.8)
dev.off()

# --- Step 8: Export results ---
results <- data.frame(
    locus = 1:length(pvals),
    pvalue = pvals,
    qvalue = qvals,
    outlier = qvals < 0.05
)
results <- results[order(results$pvalue), ]

write.csv(results, 'pcadapt_results.csv', row.names = FALSE)
write.csv(results[results$outlier, ], 'pcadapt_outliers.csv', row.names = FALSE)

cat('\nResults written to pcadapt_results.csv and pcadapt_outliers.csv\n')
