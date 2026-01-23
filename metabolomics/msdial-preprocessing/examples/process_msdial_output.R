library(tidyverse)

# === CONFIGURATION ===
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. LOAD MS-DIAL OUTPUT ===
cat('Loading MS-DIAL alignment results...\n')

# Simulate MS-DIAL output structure
# In practice: msdial_data <- read.csv('msdial_alignment_result.csv', check.names = FALSE)
# For realistic data, download from Metabolomics Workbench (www.metabolomicsworkbench.org)
# or use mtbls2 package: BiocManager::install('mtbls2'); library(mtbls2)

set.seed(42)
n_features <- 500
n_samples <- 12

# Create simulated MS-DIAL output
msdial_data <- data.frame(
    `Alignment ID` = 1:n_features,
    `Average Rt(min)` = runif(n_features, 1, 15),
    `Average Mz` = runif(n_features, 100, 800),
    `Metabolite name` = c(
        paste0('Metabolite_', 1:200),
        paste0('Lipid_', 1:150),
        rep('Unknown', 150)
    ),
    `Adduct type` = sample(c('[M+H]+', '[M+Na]+', '[M+NH4]+'), n_features, replace = TRUE),
    `Fill %` = sample(30:100, n_features, replace = TRUE),
    `MS/MS assigned` = sample(c(TRUE, FALSE), n_features, replace = TRUE, prob = c(0.3, 0.7)),
    `Annotation tag` = c(rep('Metabolite', 200), rep('Lipid', 150), rep('Unknown', 150)),
    `Formula` = NA,
    `INCHIKEY` = NA,
    check.names = FALSE
)

# Add sample intensity columns
sample_names <- paste0('Sample', 1:n_samples, ' Area')
conditions <- rep(c('Control', 'Treatment'), each = 6)

for (i in 1:n_samples) {
    base_intensity <- rlnorm(n_features, 12, 1.5)
    # Add treatment effect to first 50 features
    if (conditions[i] == 'Treatment') {
        base_intensity[1:50] <- base_intensity[1:50] * 2
    }
    msdial_data[[sample_names[i]]] <- round(base_intensity)
}

cat('Loaded', n_features, 'features from', n_samples, 'samples\n')

# === 2. EXTRACT COMPONENTS ===
cat('Extracting feature info and intensities...\n')

# Identify sample columns
sample_cols <- grep('Area$', colnames(msdial_data), value = TRUE)
meta_cols <- setdiff(colnames(msdial_data), sample_cols)

# Feature metadata
feature_info <- msdial_data[, meta_cols]

# Intensity matrix
intensity_matrix <- as.matrix(msdial_data[, sample_cols])
rownames(intensity_matrix) <- msdial_data$`Alignment ID`
colnames(intensity_matrix) <- gsub(' Area$', '', colnames(intensity_matrix))

# === 3. FILTER FEATURES ===
cat('Filtering features...\n')

# Filter by fill percentage
fill_threshold <- 50
good_fill <- feature_info$`Fill %` >= fill_threshold

# Filter by intensity (remove very low features)
min_intensity <- 1000
good_intensity <- apply(intensity_matrix, 1, max) >= min_intensity

# Combine filters
keep_features <- good_fill & good_intensity

filtered_matrix <- intensity_matrix[keep_features, ]
filtered_info <- feature_info[keep_features, ]

cat('After filtering:', sum(keep_features), '/', n_features, 'features\n')

# === 4. SUMMARIZE ANNOTATIONS ===
cat('\nAnnotation summary:\n')
print(table(filtered_info$`Annotation tag`))

# === 5. LOG TRANSFORM AND NORMALIZE ===
cat('\nNormalizing data...\n')

# Log2 transform
log_matrix <- log2(filtered_matrix + 1)

# Median normalization
sample_medians <- apply(log_matrix, 2, median)
normalized <- sweep(log_matrix, 2, sample_medians - median(sample_medians))

# === 6. QC PLOTS ===
cat('Generating QC plots...\n')

# PCA
pca <- prcomp(t(normalized), scale. = TRUE)
pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Sample = colnames(normalized),
    Condition = conditions
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 4) +
    theme_bw() +
    labs(title = 'PCA of MS-DIAL Preprocessed Data')
ggsave(file.path(output_dir, 'pca_plot.png'), width = 8, height = 6)

# Feature distribution by annotation
ggplot(filtered_info, aes(x = `Annotation tag`, fill = `Annotation tag`)) +
    geom_bar() +
    theme_bw() +
    labs(title = 'Features by Annotation Type', x = 'Annotation', y = 'Count')
ggsave(file.path(output_dir, 'annotation_distribution.png'), width = 6, height = 5)

# Retention time vs m/z
ggplot(filtered_info, aes(x = `Average Rt(min)`, y = `Average Mz`, color = `Annotation tag`)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    labs(title = 'Feature Distribution', x = 'Retention Time (min)', y = 'm/z')
ggsave(file.path(output_dir, 'rt_mz_plot.png'), width = 10, height = 6)

# === 7. EXPORT FOR DOWNSTREAM ANALYSIS ===
cat('Exporting results...\n')

# Feature table with metadata
export_table <- cbind(filtered_info, normalized)
write.csv(export_table, file.path(output_dir, 'msdial_processed.csv'), row.names = FALSE)

# Matrix only (for statistical analysis)
write.csv(normalized, file.path(output_dir, 'intensity_matrix.csv'))

# MetaboAnalyst format (samples as rows)
metaboanalyst_format <- data.frame(
    Sample = colnames(normalized),
    Group = conditions,
    t(normalized)
)
write.csv(metaboanalyst_format, file.path(output_dir, 'for_metaboanalyst.csv'), row.names = FALSE)

cat('Results saved to', output_dir, '\n')
