library(SNFtool)

# Generate synthetic multi-omics data
set.seed(42)
n <- 100
p1 <- 200
p2 <- 150

# Three clusters with different signals in each omics
true_clusters <- rep(1:3, length.out = n)

data1 <- matrix(rnorm(n * p1), nrow = n)
data2 <- matrix(rnorm(n * p2), nrow = n)

# Add cluster signal
for (k in 1:3) {
    idx <- which(true_clusters == k)
    data1[idx, 1:20] <- data1[idx, 1:20] + k
    data2[idx, 1:15] <- data2[idx, 1:15] + k
}

rownames(data1) <- rownames(data2) <- paste0('Sample', 1:n)

# SNF pipeline
dist1 <- dist2(data1, data1)
dist2_mat <- dist2(data2, data2)

K <- 20
alpha <- 0.5

aff1 <- affinityMatrix(dist1, K, alpha)
aff2 <- affinityMatrix(dist2_mat, K, alpha)

fused <- SNF(list(aff1, aff2), K = K, t = 20)

# Cluster
num_clusters <- estimateNumberOfClustersGivenGraph(fused, NUMC = 2:5)[[1]]
clusters <- spectralClustering(fused, 3)

# Evaluate
nmi <- calNMI(clusters, true_clusters)
cat('Number of clusters:', 3, '\n')
cat('NMI with true labels:', round(nmi, 3), '\n')

# Compare single vs fused
nmi_single1 <- calNMI(spectralClustering(aff1, 3), true_clusters)
nmi_single2 <- calNMI(spectralClustering(aff2, 3), true_clusters)
cat('NMI data1 only:', round(nmi_single1, 3), '\n')
cat('NMI data2 only:', round(nmi_single2, 3), '\n')
cat('NMI fused:', round(nmi, 3), '\n')
