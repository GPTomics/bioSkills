# Reference: MetaCycle 1.2+, data.table 1.14+ | Verify API if version differs
# meta2d writes an OUTPUT DIRECTORY of result CSVs; this script uses tempdir() and removes them so no strays remain.
library(MetaCycle)
library(data.table)

set.seed(42)

# --- Simulate circadian expression data ---
# 48h sampled every 4h = 13 timepoints, 2 complete 24h cycles (>=2 cycles, >=6/cycle design minima).
timepoints <- seq(0, 48, by = 4)
n_genes <- 200
n_rhythmic <- 50

expression_mat <- matrix(nrow = n_genes, ncol = length(timepoints))
rownames(expression_mat) <- paste0('gene_', seq_len(n_genes))
colnames(expression_mat) <- paste0('ZT', timepoints)

for (i in seq_len(n_genes)) {
    mesor <- runif(1, 5, 12)
    if (i <= n_rhythmic) {
        amplitude <- runif(1, 1.0, 3.0)                  # moderate circadian amplitude
        phase <- runif(1, 0, 2 * pi)
        values <- mesor + amplitude * cos(2 * pi * timepoints / 24 - phase)
    } else {
        values <- rep(mesor, length(timepoints))
    }
    expression_mat[i, ] <- values + rnorm(length(timepoints), 0, 0.5)  # SD 0.5: normalized log-expression noise
}

# --- Write input and run MetaCycle into a temp output directory ---
input_df <- data.frame(GeneID = rownames(expression_mat), expression_mat, check.names = FALSE)
input_file <- tempfile(fileext = '.csv')
write.csv(input_df, input_file, row.names = FALSE)
output_dir <- file.path(tempdir(), 'metaout')
dir.create(output_dir, showWarnings = FALSE)

# cycMethod='JTK': non-parametric, robust, fast for evenly sampled integer-hour data.
# minper=maxper=24: this is a KNOWN-PERIOD test on entrained (LD) data, not a period search.
# For free-running (DD) data where tau != 24h, widen to minper=22, maxper=26.
meta2d(infile = input_file, filestyle = 'csv', outdir = output_dir,
       timepoints = timepoints, cycMethod = 'JTK', minper = 24, maxper = 24,
       outputFile = TRUE, outRawData = FALSE)

results <- fread(file.path(output_dir, paste0('meta2d_', basename(input_file))))

# meta2d_BH.Q ranks rhythmicity; pair it with relative amplitude (meta2d_rAMP) as an effect-size filter,
# since significance alone over-detects (Laloum & Robinson-Rechavi 2020).
rhythmic <- results[meta2d_BH.Q < 0.05 & meta2d_rAMP > 0.1]
gene_ids <- as.integer(gsub('gene_', '', rhythmic$CycID))
cat(sprintf('Rhythmic (BH.Q<0.05 & rAMP>0.1): %d / %d; true positives: %d / %d\n',
            nrow(rhythmic), nrow(results), sum(gene_ids <= n_rhythmic), n_rhythmic))

# --- Plot top gene + phase distribution into the temp dir, then clean up ---
top_gene <- rhythmic[which.min(meta2d_BH.Q), CycID]
top_idx <- which(rownames(expression_mat) == top_gene)
top_row <- rhythmic[CycID == top_gene]
plot_file <- file.path(output_dir, 'jtk_cycle_results.pdf')

pdf(plot_file, width = 10, height = 5)
par(mfrow = c(1, 2))
plot(timepoints, expression_mat[top_idx, ], pch = 19, col = 'steelblue', cex = 1.2,
     xlab = 'Time (hours)', ylab = 'Expression',
     main = sprintf('%s (q = %.2e)', top_gene, top_row$meta2d_BH.Q))
t_fine <- seq(0, 48, length.out = 200)
# meta2d_phase is peak time in HOURS from ZT0: reconstruct as Base + AMP*cos(2*pi*(t - phase)/period).
fitted <- top_row$meta2d_Base + top_row$meta2d_AMP *
    cos(2 * pi * (t_fine - top_row$meta2d_phase) / top_row$meta2d_period)
lines(t_fine, fitted, col = 'red', lwd = 2)
hist(rhythmic$meta2d_phase, breaks = seq(0, 24, by = 1), col = 'coral', border = 'black',
     xlab = 'Peak time (hours from ZT0)', ylab = 'Number of genes', main = 'Phase distribution')
dev.off()

cat(sprintf('Plot written to %s\n', plot_file))
unlink(output_dir, recursive = TRUE)
unlink(input_file)
