# Reference: limma 3.58+, Mfuzz 2.62+, mgcv 1.9+, clusterProfiler 4.10+, MetaCycle 1.2+ | Verify API if version differs
## Time-course pipeline: temporal DE (limma splines) -> Mfuzz soft clustering ->
## mgcv GAM trajectory -> per-cluster clusterProfiler enrichment (temporal-gene background).
## Self-contained demo: synthetic data, outputs to a temp dir removed on exit (no stray files).
## The MetaCycle rhythm branch is OPTIONAL and gated OFF by default (see CIRCADIAN_DESIGN).

library(limma)
library(splines)
library(Mfuzz)
library(mgcv)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cluster)

# --- Configuration ---
FDR_THRESHOLD <- 0.05     # standard temporal-DE threshold; 0.1 for exploratory clustering only
N_CLUSTERS <- 4           # a CHOICE not a result; matches the 4 synthetic archetypes, sweep on real data
GAM_K <- 5                # GAM basis-dimension CEILING; keep < number of unique timepoints
CIRCADIAN_DESIGN <- FALSE # set TRUE ONLY under a real circadian design (see the gate below)
set.seed(42)

workdir <- tempfile('timecourse_')
dir.create(workdir)

# --- Step 0: synthetic data (8 timepoints x 3 reps; 4 archetypes + flat genes; real symbols for enrichment) ---
timepoints <- c(0, 3, 6, 9, 12, 18, 24, 36)
n_reps <- 3
times <- rep(timepoints, each = n_reps)
tnorm <- scale(timepoints)[, 1]
archetypes <- list(early_transient = exp(-((tnorm + 1)^2)),
                   late_sustained = 1 / (1 + exp(-3 * tnorm)),
                   monotone_down = -tnorm,
                   biphasic = sin(1.5 * tnorm))
all_symbols <- keys(org.Hs.eg.db, keytype = 'SYMBOL')
gene_ids <- sample(all_symbols, 240)
expr <- matrix(0, nrow = 240, ncol = length(times), dimnames = list(gene_ids, paste0('s', seq_along(times))))
gi <- 1
for (shape in archetypes) {
    for (k in 1:40) {
        expr[gi, ] <- rep(2 * shape, each = n_reps) + rnorm(length(times), sd = 0.3)
        gi <- gi + 1
    }
}
for (k in 1:80) expr[gi + k - 1, ] <- runif(1, 3, 7) + rnorm(length(times), sd = 0.3)
meta <- data.frame(sample = colnames(expr), time = times)
message(sprintf('Synthetic input: %d genes x %d samples, %d timepoints', nrow(expr), ncol(expr), length(timepoints)))

# --- Step 1: temporal DE (limma splines) ---
design <- model.matrix(~ ns(meta$time, df = 3))   # df=3 cubic; raise to 4-5 for >10 timepoints
fit <- lmFit(expr, design)
fit <- eBayes(fit)
temporal_results <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = 'F')
# topTable returns adj.P.Val (BH-corrected); use it directly (not $FDR, which does not exist)
sig_genes <- rownames(temporal_results[temporal_results$adj.P.Val < FDR_THRESHOLD, ])
message(sprintf('Significant temporal genes (FDR <%s): %d', FDR_THRESHOLD, length(sig_genes)))
if (length(sig_genes) < 100) message('WARNING: Few temporal genes. On real data, check replicates or relax FDR.')

# --- Step 2: filter to temporal genes; average replicates to timepoint means for clustering ---
expr_sig <- expr[sig_genes, ]
means <- t(apply(expr_sig, 1, function(y) tapply(y, meta$time, mean)))

# --- Step 3: Mfuzz soft clustering on z-scored profiles ---
eset <- ExpressionSet(assayData = means)
eset <- standardise(eset)   # per-gene mean 0, sd 1: distance is shape-based, not magnitude-based
# mestimate() (Schwaemmle & Jensen 2010) returns the smallest m that stops clustering of RANDOM data;
# it is dominated by the number of timepoints, so inspect it rather than hardcoding m=2.
m <- mestimate(eset)
message(sprintf('Estimated fuzzifier m = %.2f', m))
cl <- mfuzz(eset, c = N_CLUSTERS, m = m)
core_genes <- acore(eset, cl, min.acore = 0.5)   # membership >0.5 = core (confident) genes

cluster_sizes <- table(cl$cluster)
print(cluster_sizes)
if (any(cluster_sizes == 0)) message('WARNING: Empty clusters. Reduce N_CLUSTERS.')
sil <- silhouette(cl$cluster, dist(exprs(eset)))
message(sprintf('Mean silhouette: %.3f', mean(sil[, 3])))
write.csv(data.frame(gene = names(cl$cluster), cluster = cl$cluster),
          file.path(workdir, 'clusters.csv'), row.names = FALSE)

# --- Step 4a: OPTIONAL rhythm detection - GATED (skipped unless the design licenses it) ---
n_cycles <- (max(timepoints) - min(timepoints)) / 24
samples_per_cycle <- length(timepoints) / max(n_cycles, 1e-9)
gate_design <- n_cycles >= 2 && samples_per_cycle >= 6   # the COMPUTABLE part of the gate
if (CIRCADIAN_DESIGN && gate_design) {
    library(MetaCycle)
    meta_in <- means
    colnames(meta_in) <- sort(unique(meta$time))
    write.csv(meta_in, file.path(workdir, 'for_metacycle.csv'))
    # ARS/JTK need EVEN integer sampling and drop out silently on uneven data, leaving LS only
    meta2d(file.path(workdir, 'for_metacycle.csv'), filestyle = 'csv', minper = 20, maxper = 28,
           timepoints = sort(unique(meta$time)), outdir = file.path(workdir, 'metacycle'))
    message('MetaCycle rhythm detection complete (filter on meta2d_BH.Q AND meta2d_rAMP).')
} else if (!gate_design) {
    message(sprintf('Rhythm detection SKIPPED: inadequate design (%.1f cycles, %.1f samples/cycle; need >=2 and >=6-8/cycle).',
                    n_cycles, samples_per_cycle))
} else {
    message('Rhythm detection SKIPPED: design meets the cycle/sampling floor but CIRCADIAN_DESIGN is not set (randomized-order / circadian precondition unconfirmed).')
}

# --- Step 4b: GAM trajectory per cluster (standardized cluster-mean profiles -> Gaussian is fine) ---
for (cl_id in 1:N_CLUSTERS) {
    cl_names <- names(cl$cluster[cl$cluster == cl_id])
    if (length(cl_names) == 0) next
    mean_profile <- colMeans(means[cl_names, , drop = FALSE])
    df_gam <- data.frame(time = sort(unique(meta$time)), expr = mean_profile)
    # k is a flexibility CEILING (max basis dimension), NOT the number of knots/bends. REML (not GCV)
    # picks the wiggliness penalty; realized complexity is edf. Keep k < number of unique timepoints.
    gam_fit <- gam(expr ~ s(time, k = GAM_K), data = df_gam, method = 'REML')
    message(sprintf('Cluster %d: GAM R^2 = %.3f, EDF = %.2f (edf~1 => linear; edf~k-1 => highly non-linear)',
                    cl_id, summary(gam_fit)$r.sq, summary(gam_fit)$edf))
}

# --- Step 5: per-cluster enrichment (background = temporal genes, NOT the genome) ---
all_temporal_entrez <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
clusters_with_terms <- 0
for (i in seq_along(core_genes)) {
    genes <- core_genes[[i]]$NAME
    if (length(genes) < 5) next
    entrez <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene = entrez$ENTREZID, universe = all_temporal_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    n_terms <- nrow(as.data.frame(ego))
    if (n_terms > 0) {
        ego <- simplify(ego, cutoff = 0.7, by = 'p.adjust')   # collapse redundant parent-child GO terms
        clusters_with_terms <- clusters_with_terms + 1
    }
    message(sprintf('Cluster %d: %d significant GO BP terms', i, n_terms))
}
message(sprintf('Clusters with significant GO terms: %d / %d', clusters_with_terms, length(core_genes)))
if (clusters_with_terms < 3) message('WARNING: Few clusters enriched. On real data, check gene ID mapping or thresholds.')

unlink(workdir, recursive = TRUE)   # remove all pipeline outputs; leave no stray files
message(sprintf('Pipeline complete: %d temporal genes, %d clusters, %d enriched',
                length(sig_genes), N_CLUSTERS, clusters_with_terms))
