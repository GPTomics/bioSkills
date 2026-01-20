library(lipidr)
library(ggplot2)
library(dplyr)
library(tidyr)

# === CONFIGURATION ===
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. LOAD DATA ===
cat('Loading lipidomics data...\n')

# From LipidSearch export (or create from matrix)
# lipid_data <- read_lipidomes('lipidsearch_export.csv', data_type = 'LipidSearch')

# Example: Create from intensity matrix
set.seed(42)
n_samples <- 12
n_lipids <- 100
sample_names <- paste0('Sample', 1:n_samples)
conditions <- rep(c('Control', 'Treatment'), each = 6)

lipid_names <- c(
    paste0('PC(', 30:39, ':', sample(0:2, 10, replace = TRUE), ')'),
    paste0('PE(', 32:41, ':', sample(0:3, 10, replace = TRUE), ')'),
    paste0('TG(', 48:57, ':', sample(0:4, 10, replace = TRUE), ')'),
    paste0('SM(d18:1/', 14:23, ':0)'),
    paste0('Cer(d18:1/', 16:25, ':0)'),
    paste0('PS(', 34:43, ':', sample(0:2, 10, replace = TRUE), ')'),
    paste0('PI(', 34:43, ':', sample(0:3, 10, replace = TRUE), ')'),
    paste0('PG(', 32:41, ':', sample(0:2, 10, replace = TRUE), ')'),
    paste0('DG(', 32:41, ':', sample(0:2, 10, replace = TRUE), ')'),
    paste0('CE(', 16:25, ':', sample(0:1, 10, replace = TRUE), ')')
)

intensity_matrix <- matrix(rlnorm(n_samples * n_lipids, 10, 1), nrow = n_lipids)
colnames(intensity_matrix) <- sample_names
rownames(intensity_matrix) <- lipid_names

# Add treatment effect to some lipids
intensity_matrix[1:20, 7:12] <- intensity_matrix[1:20, 7:12] * 2

cat('Loaded', n_lipids, 'lipids from', n_samples, 'samples\n')

# === 2. PARSE LIPID ANNOTATIONS ===
cat('Parsing lipid annotations...\n')

parse_lipid <- function(name) {
    pattern <- '(\\w+)\\(.*?(\\d+):(\\d+)'
    match <- regmatches(name, regexec(pattern, name))[[1]]
    if (length(match) >= 4) {
        return(data.frame(
            lipid = name,
            class = match[2],
            carbons = as.numeric(match[3]),
            unsaturation = as.numeric(match[4])
        ))
    }
    return(data.frame(lipid = name, class = NA, carbons = NA, unsaturation = NA))
}

lipid_info <- do.call(rbind, lapply(lipid_names, parse_lipid))
cat('Lipid classes:', paste(unique(lipid_info$class), collapse = ', '), '\n')

# === 3. NORMALIZATION ===
cat('Normalizing data...\n')
log_data <- log2(intensity_matrix + 1)
sample_medians <- apply(log_data, 2, median)
normalized <- sweep(log_data, 2, sample_medians - median(sample_medians))

# === 4. DIFFERENTIAL ANALYSIS ===
cat('Running differential analysis...\n')

results <- data.frame(lipid = rownames(normalized))
results$class <- lipid_info$class
results$carbons <- lipid_info$carbons
results$unsaturation <- lipid_info$unsaturation

ctrl_cols <- which(conditions == 'Control')
treat_cols <- which(conditions == 'Treatment')

results$mean_ctrl <- rowMeans(normalized[, ctrl_cols])
results$mean_treat <- rowMeans(normalized[, treat_cols])
results$logFC <- results$mean_treat - results$mean_ctrl

pvals <- apply(normalized, 1, function(x) {
    t.test(x[treat_cols], x[ctrl_cols])$p.value
})
results$pvalue <- pvals
results$padj <- p.adjust(pvals, method = 'BH')
results$significant <- abs(results$logFC) > 1 & results$padj < 0.05

cat('Significant lipids:', sum(results$significant), '\n')

# === 5. VISUALIZATION ===
cat('Generating plots...\n')

# Volcano plot
ggplot(results, aes(x = logFC, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    scale_color_manual(values = c('gray', 'red')) +
    theme_bw() +
    labs(title = 'Lipidomics Differential Analysis', x = 'Log2 Fold Change', y = '-Log10(adj. p-value)')
ggsave(file.path(output_dir, 'volcano_plot.png'), width = 8, height = 6)

# Class composition
class_summary <- results %>%
    group_by(class) %>%
    summarise(n_lipids = n(), n_sig = sum(significant), mean_fc = mean(logFC))

ggplot(class_summary, aes(x = reorder(class, -n_lipids), y = n_lipids, fill = class)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Lipids per Class', x = 'Lipid Class', y = 'Count')
ggsave(file.path(output_dir, 'class_distribution.png'), width = 8, height = 6)

# Saturation analysis
results$sat_class <- ifelse(results$unsaturation == 0, 'Saturated',
                            ifelse(results$unsaturation == 1, 'MUFA', 'PUFA'))

sat_summary <- results %>%
    group_by(sat_class) %>%
    summarise(mean_fc = mean(logFC), se = sd(logFC) / sqrt(n()))

ggplot(sat_summary, aes(x = sat_class, y = mean_fc, fill = sat_class)) +
    geom_bar(stat = 'identity') +
    geom_errorbar(aes(ymin = mean_fc - se, ymax = mean_fc + se), width = 0.2) +
    theme_bw() +
    labs(title = 'Mean Fold Change by Saturation', x = 'Saturation', y = 'Mean Log2 FC')
ggsave(file.path(output_dir, 'saturation_analysis.png'), width = 6, height = 5)

# === 6. SAVE RESULTS ===
write.csv(results, file.path(output_dir, 'lipidomics_results.csv'), row.names = FALSE)
cat('Results saved to', output_dir, '\n')
