library(tidyverse)

# Load data
data <- read.csv('feature_table.csv', row.names = 1)
sample_info <- read.csv('sample_info.csv')

# Identify QC samples
qc_idx <- sample_info$sample_type == 'QC'
qc_names <- sample_info$sample_name[qc_idx]

cat('Samples:', nrow(data), '\n')
cat('QC samples:', sum(qc_idx), '\n')
cat('Features:', ncol(data), '\n')

# 1. Filter high-missing features
missing_pct <- colMeans(is.na(data))
data_filt <- data[, missing_pct < 0.2]
cat('Features after filtering:', ncol(data_filt), '\n')

# 2. Log transform
data_log <- log2(data_filt + 1)

# 3. PQN normalization
reference <- apply(data_log, 2, median, na.rm = TRUE)
quotients <- data_log / reference
factors <- apply(quotients, 1, median, na.rm = TRUE)
data_norm <- data_log / factors

# 4. Calculate QC RSD
qc_data <- data_norm[qc_names, ]
rsd <- apply(qc_data, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)

cat('Median QC RSD:', round(median(rsd, na.rm = TRUE), 1), '%\n')
cat('Features with RSD <30%:', sum(rsd < 30, na.rm = TRUE), '\n')

# 5. Pareto scale
data_centered <- scale(data_norm, center = TRUE, scale = FALSE)
data_scaled <- data_centered / sqrt(apply(data_norm, 2, sd, na.rm = TRUE))

# Save
write.csv(data_scaled, 'normalized_data.csv')
cat('Saved normalized data\n')
