library(flowCore)
library(ggplot2)

# === 1. LOAD DATA ===
cat('Loading FCS files...\n')

# Simulate FCS data (in practice: fs <- read.flowSet(...))
set.seed(42)
n_cells <- 50000

# Generate singlets (linear FSC-A vs FSC-H relationship)
fsc_h_singlets <- rnorm(n_cells * 0.95, 100000, 20000)
fsc_a_singlets <- fsc_h_singlets * 1.1 + rnorm(n_cells * 0.95, 0, 5000)

# Generate doublets (elevated FSC-A)
n_doublets <- n_cells * 0.05
fsc_h_doublets <- rnorm(n_doublets, 100000, 20000)
fsc_a_doublets <- fsc_h_doublets * 1.5 + rnorm(n_doublets, 20000, 5000)

# Combine
fsc_h <- c(fsc_h_singlets, fsc_h_doublets)
fsc_a <- c(fsc_a_singlets, fsc_a_doublets)
true_doublet <- c(rep(FALSE, length(fsc_h_singlets)), rep(TRUE, length(fsc_h_doublets)))

data <- data.frame(FSC_A = pmax(fsc_a, 0), FSC_H = pmax(fsc_h, 0), true_doublet = true_doublet)
cat('Total events:', nrow(data), '\n')
cat('True doublets:', sum(true_doublet), '(', round(mean(true_doublet) * 100, 1), '%)\n')

# === 2. RATIO-BASED DETECTION ===
cat('\nMethod 1: Ratio-based detection...\n')

data$ratio <- data$FSC_A / (data$FSC_H + 1)
ratio_threshold <- quantile(data$ratio, 0.95)
data$doublet_ratio <- data$ratio > ratio_threshold

sensitivity_ratio <- sum(data$doublet_ratio & data$true_doublet) / sum(data$true_doublet)
specificity_ratio <- sum(!data$doublet_ratio & !data$true_doublet) / sum(!data$true_doublet)

cat('  Threshold:', round(ratio_threshold, 3), '\n')
cat('  Sensitivity:', round(sensitivity_ratio * 100, 1), '%\n')
cat('  Specificity:', round(specificity_ratio * 100, 1), '%\n')

# === 3. RESIDUAL-BASED DETECTION ===
cat('\nMethod 2: Residual-based detection...\n')

fit <- lm(FSC_A ~ FSC_H, data = data)
data$residual <- abs(data$FSC_A - predict(fit))
residual_threshold <- quantile(data$residual, 0.95)
data$doublet_residual <- data$residual > residual_threshold

sensitivity_resid <- sum(data$doublet_residual & data$true_doublet) / sum(data$true_doublet)
specificity_resid <- sum(!data$doublet_residual & !data$true_doublet) / sum(!data$true_doublet)

cat('  Threshold:', round(residual_threshold, 0), '\n')
cat('  Sensitivity:', round(sensitivity_resid * 100, 1), '%\n')
cat('  Specificity:', round(specificity_resid * 100, 1), '%\n')

# === 4. VISUALIZATION ===
cat('\nGenerating plots...\n')

# Scatter plot with true labels
p1 <- ggplot(data, aes(x = FSC_H, y = FSC_A, color = true_doublet)) +
    geom_point(alpha = 0.2, size = 0.5) +
    scale_color_manual(values = c('gray40', 'red'), labels = c('Singlet', 'Doublet')) +
    theme_bw() +
    labs(title = 'True Doublet Status', x = 'FSC-H', y = 'FSC-A', color = 'Status')

# Scatter plot with detected labels
p2 <- ggplot(data, aes(x = FSC_H, y = FSC_A, color = doublet_residual)) +
    geom_point(alpha = 0.2, size = 0.5) +
    scale_color_manual(values = c('gray40', 'blue'), labels = c('Singlet', 'Doublet')) +
    theme_bw() +
    labs(title = 'Detected Doublets (Residual Method)', x = 'FSC-H', y = 'FSC-A', color = 'Status')

# Ratio distribution
p3 <- ggplot(data, aes(x = ratio, fill = true_doublet)) +
    geom_histogram(bins = 100, alpha = 0.7, position = 'identity') +
    geom_vline(xintercept = ratio_threshold, linetype = 'dashed', color = 'red') +
    scale_fill_manual(values = c('gray40', 'red'), labels = c('Singlet', 'Doublet')) +
    theme_bw() +
    labs(title = 'FSC Ratio Distribution', x = 'FSC-A / FSC-H Ratio', y = 'Count', fill = 'True Status')

# Save plots
ggsave('doublet_true.png', p1, width = 8, height = 6)
ggsave('doublet_detected.png', p2, width = 8, height = 6)
ggsave('doublet_ratio_dist.png', p3, width = 8, height = 5)

# === 5. SUMMARY ===
cat('\n=== SUMMARY ===\n')
cat('Total events:', nrow(data), '\n')
cat('True doublets:', sum(data$true_doublet), '\n')
cat('Detected doublets (residual):', sum(data$doublet_residual), '\n')
cat('True positives:', sum(data$doublet_residual & data$true_doublet), '\n')
cat('False positives:', sum(data$doublet_residual & !data$true_doublet), '\n')

# Final singlet data
singlets <- data[!data$doublet_residual, ]
cat('\nFinal singlet count:', nrow(singlets), '\n')
