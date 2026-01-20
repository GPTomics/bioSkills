library(xcms)
library(MSnbase)
library(limma)
library(ggplot2)

# === CONFIGURATION ===
data_dir <- 'data/'
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. LOAD DATA ===
cat('Loading data...\n')
mzml_files <- list.files(data_dir, pattern = '\\.mzML$', full.names = TRUE)
sample_data <- read.csv('sample_metadata.csv')
raw_data <- readMSData(mzml_files, mode = 'onDisk')
pData(raw_data) <- sample_data
cat('Loaded', length(mzml_files), 'samples\n')

# === 2. PEAK DETECTION ===
cat('Detecting peaks...\n')
cwp <- CentWaveParam(peakwidth = c(5, 30), ppm = 25, snthresh = 10,
                     prefilter = c(3, 1000), noise = 1000)
xdata <- findChromPeaks(raw_data, param = cwp)
cat('Detected', nrow(chromPeaks(xdata)), 'peaks\n')

# === 3. RETENTION TIME ALIGNMENT ===
cat('Aligning...\n')
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

# === 4. FEATURE GROUPING ===
cat('Grouping features...\n')
pdp <- PeakDensityParam(sampleGroups = pData(xdata)$condition,
                        minFraction = 0.5, bw = 5, binSize = 0.025)
xdata <- groupChromPeaks(xdata, param = pdp)
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
cat('Total features:', nrow(featureDefinitions(xdata)), '\n')

# === 5. EXTRACT & NORMALIZE ===
feature_matrix <- featureValues(xdata, value = 'into', method = 'maxint')
feature_matrix[feature_matrix == 0] <- NA
log_matrix <- log2(feature_matrix)

valid_features <- rowSums(!is.na(log_matrix)) > ncol(log_matrix) * 0.5
filtered_matrix <- log_matrix[valid_features, ]
cat('Features after filtering:', nrow(filtered_matrix), '\n')

sample_medians <- apply(filtered_matrix, 2, median, na.rm = TRUE)
normalized <- sweep(filtered_matrix, 2, sample_medians - median(sample_medians))

# === 6. QC PLOTS ===
cat('Generating QC plots...\n')
pca <- prcomp(t(normalized), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Condition = pData(xdata)$condition)
ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) + theme_bw() + labs(title = 'Metabolomics PCA')
ggsave(file.path(output_dir, 'qc_pca.png'), width = 8, height = 6)

# === 7. DIFFERENTIAL ANALYSIS ===
cat('Running differential analysis...\n')
design <- model.matrix(~ 0 + condition, data = pData(xdata))
colnames(design) <- levels(factor(pData(xdata)$condition))

imputed <- normalized
imputed[is.na(imputed)] <- min(imputed, na.rm = TRUE) - 1

fit <- lmFit(imputed, design)
contrast <- makeContrasts(Treatment - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))

results <- topTable(fit2, number = Inf, adjust.method = 'BH')
results$significant <- abs(results$logFC) > 1 & results$adj.P.Val < 0.05
cat('Significant features:', sum(results$significant), '\n')

# === 8. VOLCANO PLOT ===
ggplot(results, aes(logFC, -log10(adj.P.Val), color = significant)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    scale_color_manual(values = c('gray', 'red')) +
    theme_bw() + labs(title = 'Differential Metabolites')
ggsave(file.path(output_dir, 'volcano.png'), width = 8, height = 6)

# === 9. SAVE RESULTS ===
write.csv(results, file.path(output_dir, 'differential_metabolites.csv'), row.names = TRUE)
write.csv(normalized, file.path(output_dir, 'normalized_matrix.csv'))
cat('Results saved to', output_dir, '\n')
