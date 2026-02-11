## Time-course analysis pipeline: temporal DE, Mfuzz clustering, optional MetaCycle,
## GAM trajectory fitting, and per-cluster pathway enrichment.

library(limma)
library(splines)
library(Mfuzz)
library(mgcv)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cluster)

# --- Configuration ---
COUNTS_FILE <- 'counts_normalized.csv'
METADATA_FILE <- 'metadata.csv'
OUTPUT_PREFIX <- 'timecourse_results'
# FDR <0.05: standard temporal DE threshold; relax to 0.1 for exploratory analysis
FDR_THRESHOLD <- 0.05
# 4-20 clusters typical; start with 8, refine via gap statistic or silhouette
N_CLUSTERS <- 8
# Membership >0.5: core gene assignment; lower to 0.3 for exploratory
MIN_MEMBERSHIP <- 0.5
# GAM basis dimension; 5 default, reduce to 3 for <6 time points
GAM_K <- 5
# Set TRUE if experiment covers 24h+ cycles with 2-4h resolution
CIRCADIAN_DESIGN <- FALSE

# --- Step 1: Load data ---
expr <- as.matrix(read.csv(COUNTS_FILE, row.names = 1))
meta <- read.csv(METADATA_FILE)
message(sprintf('Loaded: %d genes x %d samples, %d time points',
                nrow(expr), ncol(expr), length(unique(meta$time))))

# --- Step 2: Temporal DE (limma splines) ---
# df=3: cubic spline sufficient for most time courses; increase to 4-5 for >10 time points
design <- model.matrix(~ ns(meta$time, df = 3))
fit <- lmFit(expr, design)
fit <- eBayes(fit)
temporal_results <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = 'F')
# topTable already returns adj.P.Val (BH-corrected); use it directly
sig_genes <- rownames(temporal_results[temporal_results$adj.P.Val < FDR_THRESHOLD, ])
message(sprintf('Significant temporal genes (FDR <%s): %d', FDR_THRESHOLD, length(sig_genes)))

# QC gate: sufficient temporal genes
if (length(sig_genes) < 100) {
    message('WARNING: Few temporal genes detected. Check replicate count or relax FDR.')
}

expr_sig <- expr[sig_genes, ]

# --- Step 3: Mfuzz soft clustering ---
eset <- ExpressionSet(assayData = as.matrix(expr_sig))
eset <- standardise(eset)

# mestimate(): optimal fuzzifier from data geometry; typical range 1.5-2.5
m <- mestimate(eset)
message(sprintf('Estimated fuzzifier m = %.2f', m))

cl <- mfuzz(eset, c = N_CLUSTERS, m = m)

# QC gate: no empty clusters
cluster_sizes <- table(cl$cluster)
message('Cluster sizes:')
print(cluster_sizes)
if (any(cluster_sizes == 0)) {
    message('WARNING: Empty clusters found. Reduce N_CLUSTERS.')
}

core_genes <- acore(eset, cl, min.acore = MIN_MEMBERSHIP)
for (i in seq_along(core_genes)) {
    message(sprintf('Cluster %d: %d core genes (membership >%.1f)', i, nrow(core_genes[[i]]), MIN_MEMBERSHIP))
}

# Silhouette validation
sil <- silhouette(cl$cluster, dist(exprs(eset)))
message(sprintf('Mean silhouette score: %.3f', mean(sil[, 3])))

# Save cluster assignments
cluster_df <- data.frame(gene = names(cl$cluster), cluster = cl$cluster)
write.csv(cluster_df, paste0(OUTPUT_PREFIX, '_clusters.csv'), row.names = FALSE)

# --- Step 4a: Optional rhythm detection (MetaCycle) ---
if (CIRCADIAN_DESIGN) {
    library(MetaCycle)
    expr_for_meta <- expr_sig
    colnames(expr_for_meta) <- meta$time_hours
    write.csv(expr_for_meta, paste0(OUTPUT_PREFIX, '_for_metacycle.csv'))

    # Period range 20-28h: standard circadian window
    # Adjust to 4-12h for ultradian or >28h for infradian rhythms
    meta2d(paste0(OUTPUT_PREFIX, '_for_metacycle.csv'), filestyle = 'csv',
           minper = 20, maxper = 28,
           timepoints = sort(unique(meta$time_hours)),
           outdir = paste0(OUTPUT_PREFIX, '_metacycle'))
    message('MetaCycle rhythm detection complete.')
}

# --- Step 4b: GAM trajectory fitting ---
cluster_trajectories <- list()
for (cl_id in 1:N_CLUSTERS) {
    cl_gene_names <- names(cl$cluster[cl$cluster == cl_id])
    if (length(cl_gene_names) == 0) next

    mean_profile <- colMeans(expr_sig[cl_gene_names, , drop = FALSE])
    df_gam <- data.frame(time = meta$time, expr = mean_profile)

    gam_fit <- gam(expr ~ s(time, k = GAM_K), data = df_gam)
    r_sq <- summary(gam_fit)$r.sq
    edf <- summary(gam_fit)$edf

    cluster_trajectories[[cl_id]] <- list(fit = gam_fit, r_squared = r_sq, edf = edf)
    message(sprintf('Cluster %d: GAM R^2 = %.3f, EDF = %.2f', cl_id, r_sq, edf))
}

# --- Step 5: Per-cluster pathway enrichment ---
enrichment_results <- list()
clusters_with_terms <- 0

for (i in seq_along(core_genes)) {
    genes <- core_genes[[i]]$NAME
    if (length(genes) < 5) {
        message(sprintf('Cluster %d: too few genes (%d), skipping enrichment', i, length(genes)))
        next
    }

    entrez <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

    # pvalueCutoff 0.05, qvalueCutoff 0.05: standard enrichment thresholds
    ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = 'BP', pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                    readable = TRUE)

    enrichment_results[[i]] <- ego
    n_terms <- nrow(as.data.frame(ego))
    if (n_terms > 0) clusters_with_terms <- clusters_with_terms + 1
    message(sprintf('Cluster %d: %d significant GO BP terms', i, n_terms))
}

# QC gate: enrichment coverage
message(sprintf('Clusters with significant GO terms: %d / %d', clusters_with_terms, N_CLUSTERS))
if (clusters_with_terms < 3) {
    message('WARNING: Few clusters enriched. Check gene ID mapping or relax thresholds.')
}

# --- Save enrichment results ---
for (i in seq_along(enrichment_results)) {
    if (!is.null(enrichment_results[[i]])) {
        df_enr <- as.data.frame(enrichment_results[[i]])
        if (nrow(df_enr) > 0) {
            write.csv(df_enr, sprintf('%s_cluster%d_GO.csv', OUTPUT_PREFIX, i), row.names = FALSE)
        }
    }
}

message(sprintf('\nPipeline complete: %d temporal genes, %d clusters, %d enriched',
                length(sig_genes), N_CLUSTERS, clusters_with_terms))
