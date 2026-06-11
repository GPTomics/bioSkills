# Reference: SNFtool 2.3+ | Verify API if version differs
# Patient stratification by Similarity Network Fusion. The clustering always returns the C
# requested, so the subtype count is a claim to defend (stability across a fixed grid, a
# fused-vs-best-single-omic check, and out-of-skill survival validation), not a discovery.

library(SNFtool)

data(Data1)   # shipped demo views: two omics, 200 samples, two true groups
data(Data2)

K <- 20        # neighbors for the local kernel bandwidth; SNFtool guidance 10-30
sigma <- 0.5   # affinityMatrix width (third arg is sigma, not alpha); guidance 0.3-0.8
t_iter <- 20   # cross-diffusion iterations; converges by ~10-20

views <- lapply(list(view1=Data1, view2=Data2), standardNormalization)   # per-feature z-score before distance
dists <- lapply(views, function(x) dist2(x, x)^(1/2))                     # dist2 returns SQUARED distance; take the root
affinities <- lapply(dists, function(d) affinityMatrix(d, K, sigma))

fused <- SNF(affinities, K, t_iter)

est <- estimateNumberOfClustersGivenGraph(fused, NUMC=2:6)   # FOUR estimates (eigengap + rotation cost); plausibility, not truth
print(est)

C <- 2                                                       # the central claim - defend it, do not assume it
clusters <- spectralClustering(fused, K=C, type=3)           # here K is the CLUSTER COUNT; type 3 = Ng-Jordan-Weiss default

concordance <- concordanceNetworkNMI(c(affinities, list(fused)), C)   # did fusion beat the best single omic?
print(round(concordance, 3))

feat_rank <- rankFeaturesByNMI(views, fused)                 # POST-HOC attribution, not a model that made the clusters
cat('top view1 feature rank:', head(order(feat_rank[[1]][[1]], decreasing=TRUE)), '\n')

# Stability: resample patients, re-cluster, and compare; a real subtype recurs across subsamples.
sub <- sample(nrow(fused), 0.8 * nrow(fused))
clusters_sub <- spectralClustering(fused[sub, sub], K=C, type=3)
cat('subsample-vs-full NMI:', round(calNMI(clusters[sub], clusters_sub), 3), '\n')
