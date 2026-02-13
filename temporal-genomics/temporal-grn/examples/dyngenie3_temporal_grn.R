# Reference: numpy 1.26+, pandas 2.2+, statsmodels 0.14+ | Verify API if version differs
library(dynGENIE3)

set.seed(42)

# --- Simulate time-series expression with regulatory relationships ---
# 6 timepoints per replicate: 0, 4, 8, 12, 24, 48h
# 3 replicates: multiple series improve ODE derivative estimation
timepoints <- c(0, 4, 8, 12, 24, 48)
n_tfs <- 5
n_targets <- 20
n_genes <- n_tfs + n_targets
gene_names <- c(paste0('TF', seq_len(n_tfs)), paste0('target', seq_len(n_targets)))

simulate_series <- function(seed_offset) {
    set.seed(42 + seed_offset)
    mat <- matrix(0, nrow = n_genes, ncol = length(timepoints))
    rownames(mat) <- gene_names
    colnames(mat) <- paste0('t', timepoints)

    for (i in seq_len(n_tfs)) {
        mat[i, 1] <- runif(1, 6, 10)
        for (t in 2:length(timepoints)) {
            dt <- timepoints[t] - timepoints[t - 1]
            # Decay rate 0.05/h: moderate transcriptional decay
            mat[i, t] <- mat[i, t - 1] * exp(-0.05 * dt) + rnorm(1, 0, 0.3)
            mat[i, t] <- mat[i, t] + 2 * sin(2 * pi * timepoints[t] / 48)
        }
    }

    for (j in seq_len(n_targets)) {
        target_idx <- n_tfs + j
        mat[target_idx, 1] <- runif(1, 6, 10)
        # First 10 targets regulated by a TF with lag-proportional influence
        causal_tf <- ((j - 1) %% n_tfs) + 1
        for (t in 2:length(timepoints)) {
            if (j <= 10) {
                # 0.3 coefficient: moderate regulatory influence from TF at previous timepoint
                mat[target_idx, t] <- 0.6 * mat[target_idx, t - 1] +
                    0.3 * mat[causal_tf, t - 1] + rnorm(1, 0, 0.2)
            } else {
                mat[target_idx, t] <- 0.7 * mat[target_idx, t - 1] + rnorm(1, 0, 0.3)
            }
        }
    }
    mat
}

expr_list <- list(simulate_series(0), simulate_series(1), simulate_series(2))
time_list <- list(timepoints, timepoints, timepoints)

cat(sprintf('Input: %d genes x %d timepoints x %d replicates\n',
            n_genes, length(timepoints), length(expr_list)))

# --- Run dynGENIE3 ---
# regulators: restrict to TF indices; improves accuracy and reduces noise
tf_indices <- seq_len(n_tfs)

res <- dynGENIE3(
    TS.data = expr_list,
    time.points = time_list,
    regulators = tf_indices
)

# --- Extract top regulatory links ---
# threshold=100: top 100 edges; ~20 edges per TF for a focused network
# Increase to 500-1000 for genome-wide discovery with many TFs
link_list <- get.link.list(res$weight.matrix, report.max = 100)
cat(sprintf('\nTop 100 regulatory links:\n'))
print(head(link_list, 20))

# --- Evaluate against ground truth ---
true_edges <- data.frame(
    regulatoryGene = paste0('TF', ((seq_len(10) - 1) %% n_tfs) + 1),
    targetGene = paste0('target', seq_len(10)),
    stringsAsFactors = FALSE
)

predicted_set <- paste(link_list$regulatoryGene, link_list$targetGene, sep = '->')
true_set <- paste(true_edges$regulatoryGene, true_edges$targetGene, sep = '->')

tp <- sum(predicted_set %in% true_set)
fp <- sum(!(predicted_set %in% true_set))
fn <- sum(!(true_set %in% predicted_set))
precision <- tp / (tp + fp)
recall <- tp / (tp + fn)
cat(sprintf('\nPrecision: %.2f, Recall: %.2f\n', precision, recall))

# --- Weight matrix analysis ---
weight_mat <- res$weight.matrix
cat(sprintf('\nWeight matrix dimensions: %d regulators x %d targets\n',
            nrow(weight_mat), ncol(weight_mat)))

# Top regulators by total outgoing weight
tf_total_weight <- rowSums(weight_mat[tf_indices, ])
names(tf_total_weight) <- gene_names[tf_indices]
cat('\nTF influence ranking (total outgoing weight):\n')
print(sort(tf_total_weight, decreasing = TRUE))

# --- Visualization ---
pdf('dyngenie3_grn_results.pdf', width = 12, height = 5)
par(mfrow = c(1, 2))

top_n <- min(20, nrow(link_list))
barplot(
    rev(link_list$weight[seq_len(top_n)]),
    names.arg = rev(paste(link_list$regulatoryGene[seq_len(top_n)], '->',
                          link_list$targetGene[seq_len(top_n)])),
    horiz = TRUE, las = 1, cex.names = 0.6,
    col = ifelse(rev(predicted_set[seq_len(top_n)] %in% true_set), 'forestgreen', 'gray70'),
    xlab = 'Importance score',
    main = sprintf('Top %d edges (green = true)', top_n)
)

image(weight_mat[tf_indices, seq_len(min(20, ncol(weight_mat)))],
      axes = FALSE, col = heat.colors(50),
      main = 'TF-target weight matrix')
axis(1, at = seq(0, 1, length.out = min(20, ncol(weight_mat))),
     labels = colnames(weight_mat)[seq_len(min(20, ncol(weight_mat)))],
     las = 2, cex.axis = 0.6)
axis(2, at = seq(0, 1, length.out = n_tfs),
     labels = gene_names[tf_indices], las = 1, cex.axis = 0.8)

dev.off()
cat('\nPlot saved to dyngenie3_grn_results.pdf\n')

# --- Export results ---
write.csv(link_list, 'dyngenie3_edge_list.csv', row.names = FALSE)
cat('Edge list saved to dyngenie3_edge_list.csv\n')
