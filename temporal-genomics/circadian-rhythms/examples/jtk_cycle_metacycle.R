# Reference: R stats (base), pandas 2.2+, statsmodels 0.14+ | Verify API if version differs
library(MetaCycle)
library(data.table)

set.seed(42)

# --- Simulate circadian expression data ---
# 48h sampled every 4h = 13 timepoints, 2 complete 24h cycles
timepoints <- seq(0, 48, by = 4)
n_genes <- 200
n_rhythmic <- 50

expression_mat <- matrix(nrow = n_genes, ncol = length(timepoints))
rownames(expression_mat) <- paste0('gene_', seq_len(n_genes))
colnames(expression_mat) <- paste0('ZT', timepoints)

for (i in seq_len(n_genes)) {
    mesor <- runif(1, 5, 12)
    if (i <= n_rhythmic) {
        # 24h period: standard circadian oscillation
        amplitude <- runif(1, 1.0, 3.0)
        phase <- runif(1, 0, 2 * pi)
        values <- mesor + amplitude * cos(2 * pi * timepoints / 24 - phase)
    } else {
        values <- rep(mesor, length(timepoints))
    }
    # SD = 0.5: typical noise for normalized log-expression
    expression_mat[i, ] <- values + rnorm(length(timepoints), 0, 0.5)
}

# --- Write input file for MetaCycle ---
input_df <- data.frame(GeneID = rownames(expression_mat), expression_mat, check.names = FALSE)
input_file <- tempfile(fileext = '.csv')
write.csv(input_df, input_file, row.names = FALSE)

output_dir <- tempdir()

# --- Run MetaCycle with JTK_CYCLE ---
# cycMethod='JTK': JTK_CYCLE is non-parametric, robust, and fast for evenly sampled data
# minper=20, maxper=28: search window centered on 24h; captures near-circadian periods
meta2d(
    infile = input_file,
    filestyle = 'csv',
    outdir = output_dir,
    timepoints = timepoints,
    cycMethod = c('JTK'),
    minper = 20,
    maxper = 28,
    outputFile = TRUE,
    outRawData = FALSE
)

# --- Load and interpret results ---
result_file <- file.path(output_dir, paste0('meta2d_', basename(input_file)))
results <- fread(result_file)

# meta2d_BH.Q: BH-corrected q-value across all genes
# q < 0.05: standard FDR threshold for rhythmicity
rhythmic <- results[meta2d_BH.Q < 0.05]
cat(sprintf('Rhythmic genes (q < 0.05): %d / %d\n', nrow(rhythmic), nrow(results)))

gene_ids <- as.integer(gsub('gene_', '', rhythmic$CycID))
true_positives <- sum(gene_ids <= n_rhythmic)
cat(sprintf('True positives: %d / %d rhythmic genes in first %d\n',
            true_positives, nrow(rhythmic), n_rhythmic))

# --- Summary statistics ---
cat(sprintf('\nPhase range: %.1f - %.1f hours\n',
            min(rhythmic$meta2d_phase, na.rm = TRUE),
            max(rhythmic$meta2d_phase, na.rm = TRUE)))
cat(sprintf('Amplitude range: %.2f - %.2f\n',
            min(rhythmic$meta2d_AMP, na.rm = TRUE),
            max(rhythmic$meta2d_AMP, na.rm = TRUE)))

# --- Plot top rhythmic gene ---
top_gene <- rhythmic[which.min(meta2d_BH.Q), CycID]
top_idx <- which(rownames(expression_mat) == top_gene)

pdf('jtk_cycle_results.pdf', width = 10, height = 5)
par(mfrow = c(1, 2))

plot(timepoints, expression_mat[top_idx, ],
     pch = 19, col = 'steelblue', cex = 1.2,
     xlab = 'Time (hours)', ylab = 'Expression',
     main = sprintf('%s (q = %.2e)', top_gene,
                    rhythmic[CycID == top_gene, meta2d_BH.Q]))
t_fine <- seq(0, 48, length.out = 200)
top_row <- rhythmic[CycID == top_gene]
fitted <- top_row$meta2d_Base + top_row$meta2d_AMP *
    cos(2 * pi * t_fine / top_row$meta2d_period - 2 * pi * top_row$meta2d_phase / top_row$meta2d_period)
lines(t_fine, fitted, col = 'red', lwd = 2)

# Phase distribution of significant rhythmic genes
# Binned in 1h intervals across 24h cycle
hist(rhythmic$meta2d_phase, breaks = seq(0, 24, by = 1),
     col = 'coral', border = 'black',
     xlab = 'Phase (hours from ZT0)',
     ylab = 'Number of genes',
     main = 'Phase distribution')

dev.off()
cat('\nPlot saved to jtk_cycle_results.pdf\n')
