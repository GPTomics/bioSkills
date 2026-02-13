# Reference: numpy 1.26+, scanpy 1.10+, scikit-learn 1.4+ | Verify API if version differs
library(Mfuzz)
library(Biobase)

set.seed(42)

# --- Simulate temporal expression data ---
# 8 timepoints (0, 2, 4, 8, 12, 24, 36, 48h): common early-response design
timepoints <- c(0, 2, 4, 8, 12, 24, 36, 48)
n_genes <- 500

expr_mat <- matrix(nrow = n_genes, ncol = length(timepoints))
rownames(expr_mat) <- paste0('gene_', seq_len(n_genes))
colnames(expr_mat) <- paste0('T', timepoints, 'h')

# 4 temporal patterns: early up, late up, transient, down
for (i in seq_len(n_genes)) {
    pattern <- ((i - 1) %% 4) + 1
    base <- runif(1, 6, 10)
    if (pattern == 1) {
        signal <- c(0, 1, 2, 2.5, 2, 1.5, 1, 0.5) * runif(1, 0.8, 1.5)
    } else if (pattern == 2) {
        signal <- c(0, 0, 0.2, 0.5, 1, 2, 2.5, 3) * runif(1, 0.8, 1.5)
    } else if (pattern == 3) {
        signal <- c(0, 0.5, 2, 3, 2, 0.5, 0, 0) * runif(1, 0.8, 1.5)
    } else {
        signal <- c(0, -0.5, -1, -1.5, -2, -2.5, -2, -1.5) * runif(1, 0.8, 1.5)
    }
    # SD = 0.3: moderate noise for simulated log-expression
    expr_mat[i, ] <- base + signal + rnorm(length(timepoints), 0, 0.3)
}

# --- Create ExpressionSet and preprocess ---
eset <- ExpressionSet(assayData = expr_mat)

# min.std=0.5: removes genes with near-zero variance (flat expression across time)
eset <- filter.std(eset, min.std = 0.5)
cat(sprintf('Genes after variance filter: %d\n', nrow(exprs(eset))))

eset <- standardise(eset)

# --- Estimate fuzzifier ---
# mestimate(): data-driven fuzzifier prevents over- or under-fuzzy clusters
m <- mestimate(eset)
cat(sprintf('Estimated fuzzifier m: %.2f\n', m))

# --- Select cluster number ---
# Test k=3 to k=12; evaluate minimum centroid distance
# Below k=3 loses resolution; above k=12 rarely adds biological meaning for 8 timepoints
min_dists <- numeric()
k_range <- 3:12
for (k in k_range) {
    cl_k <- mfuzz(eset, c = k, m = m)
    dists <- as.matrix(dist(cl_k$centers))
    diag(dists) <- Inf
    min_dists <- c(min_dists, min(dists))
}

pdf('mfuzz_cluster_selection.pdf', width = 6, height = 4)
plot(k_range, min_dists, type = 'b', pch = 19, col = 'steelblue',
     xlab = 'Number of clusters (k)', ylab = 'Minimum centroid distance',
     main = 'Cluster number selection')
dev.off()

# --- Run Mfuzz with chosen k ---
# k=4: matches the 4 simulated patterns; in practice, pick from validity plot
k_chosen <- 4
cl <- mfuzz(eset, c = k_chosen, m = m)

cat(sprintf('\nCluster sizes:\n'))
for (i in seq_len(k_chosen)) {
    cat(sprintf('  Cluster %d: %d genes\n', i, sum(cl$cluster == i)))
}

# --- Filter by membership ---
# acore extracts genes with membership >= min.acore in at least one cluster
# 0.5 threshold: standard cutoff; genes below this are equidistant from multiple centroids
core <- acore(eset, cl, min.acore = 0.5)
n_core <- sum(sapply(core, nrow))
cat(sprintf('\nCore genes (membership >= 0.5): %d / %d\n', n_core, nrow(exprs(eset))))

# --- Visualization ---
pdf('mfuzz_clusters.pdf', width = 10, height = 8)
mfuzz.plot2(eset, cl, mfrow = c(2, 2), time.labels = colnames(expr_mat),
            centre = TRUE, x11 = FALSE)
dev.off()

pdf('mfuzz_overlap.pdf', width = 6, height = 6)
o <- overlap(cl)
# thres=0.05: show overlap where >5% of genes share substantial membership
overlap.plot(cl, over = o, thres = 0.05)
dev.off()

# --- Export cluster assignments ---
cluster_df <- data.frame(
    gene = names(cl$cluster),
    cluster = cl$cluster,
    max_membership = apply(cl$membership, 1, max)
)
cluster_df <- cluster_df[order(cluster_df$cluster, -cluster_df$max_membership), ]
write.csv(cluster_df, 'mfuzz_cluster_assignments.csv', row.names = FALSE)

cat('\nOutput files:\n')
cat('  mfuzz_cluster_selection.pdf\n')
cat('  mfuzz_clusters.pdf\n')
cat('  mfuzz_overlap.pdf\n')
cat('  mfuzz_cluster_assignments.csv\n')
