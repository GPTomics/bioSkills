# Complete proteomics workflow: MaxQuant to differential proteins
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# === CONFIGURATION ===
input_file <- 'proteinGroups.txt'
output_prefix <- 'proteomics_results'
fdr_threshold <- 0.05
lfc_threshold <- 1

# Sample groups (modify for your experiment)
sample_groups <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')

# === 1. DATA IMPORT ===
cat('=== Data Import ===\n')
proteins <- read.delim(input_file, stringsAsFactors = FALSE)
cat('Loaded', nrow(proteins), 'protein groups\n')

proteins <- proteins[proteins$Potential.contaminant != '+' &
                      proteins$Reverse != '+' &
                      proteins$Only.identified.by.site != '+', ]
cat('After filtering:', nrow(proteins), 'proteins\n')

lfq_cols <- grep('^LFQ\\.intensity\\.', colnames(proteins), value = TRUE)
intensities <- proteins[, lfq_cols]
rownames(intensities) <- proteins$Majority.protein.IDs
colnames(intensities) <- gsub('LFQ\\.intensity\\.', '', colnames(intensities))
cat('Samples:', ncol(intensities), '\n')

# === 2. TRANSFORM & NORMALIZE ===
cat('\n=== Normalization ===\n')
intensities[intensities == 0] <- NA
log2_int <- log2(intensities)

sample_medians <- apply(log2_int, 2, median, na.rm = TRUE)
cat('Sample medians before:', round(sample_medians, 2), '\n')
normalized <- sweep(log2_int, 2, sample_medians - median(sample_medians))
cat('Sample medians after:', round(apply(normalized, 2, median, na.rm = TRUE), 2), '\n')

# === 3. FILTER & IMPUTE ===
cat('\n=== Filtering & Imputation ===\n')
missing_pct <- rowSums(is.na(normalized)) / ncol(normalized)
filtered <- normalized[missing_pct < 0.5, ]
cat('Proteins after filtering (< 50% missing):', nrow(filtered), '\n')

impute_minprob <- function(x) {
    nas <- is.na(x)
    if (all(nas) || sum(!nas) < 2) return(x)
    x[nas] <- rnorm(sum(nas), mean(x, na.rm = TRUE) - 1.8 * sd(x, na.rm = TRUE), 0.3 * sd(x, na.rm = TRUE))
    x
}
set.seed(42)
imputed <- as.data.frame(t(apply(filtered, 1, impute_minprob)))

# === 4. QC ===
cat('\n=== Quality Control ===\n')
pca <- prcomp(t(imputed), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = rownames(pca$x), Group = sample_groups)
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4) + theme_minimal() +
    labs(x = paste0('PC1 (', var_exp[1], '%)'), y = paste0('PC2 (', var_exp[2], '%)'), title = 'PCA of Protein Abundances')
ggsave(paste0(output_prefix, '_pca.pdf'), p_pca, width = 7, height = 6)

# === 5. DIFFERENTIAL ANALYSIS ===
cat('\n=== Differential Analysis ===\n')
sample_info <- data.frame(sample = colnames(imputed), condition = factor(sample_groups, levels = c('Control', 'Treatment')))
design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(as.matrix(imputed), design)
contrast <- makeContrasts(Treatment - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))

results <- topTable(fit2, number = Inf, adjust.method = 'BH')
results$protein <- rownames(results)
results$significant <- abs(results$logFC) > lfc_threshold & results$adj.P.Val < fdr_threshold

cat('Total proteins tested:', nrow(results), '\n')
cat('Significant:', sum(results$significant), '\n')
cat('  Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('  Down-regulated:', sum(results$significant & results$logFC < 0), '\n')

# === 6. VISUALIZATION ===
p_volcano <- ggplot(results, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = 'dashed') +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = 'dashed') +
    scale_color_manual(values = c('grey60', 'firebrick')) +
    theme_minimal() + labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-value')
ggsave(paste0(output_prefix, '_volcano.pdf'), p_volcano, width = 8, height = 6)

# Heatmap of significant proteins
if (sum(results$significant) > 1) {
    sig_proteins <- rownames(results)[results$significant]
    mat <- as.matrix(imputed[sig_proteins, ])
    mat_scaled <- t(scale(t(mat)))
    annotation_col <- data.frame(Group = sample_groups, row.names = colnames(mat))
    pheatmap(mat_scaled, annotation_col = annotation_col, show_rownames = nrow(mat_scaled) < 50,
             filename = paste0(output_prefix, '_heatmap.pdf'), width = 8, height = 10)
}

# === 7. EXPORT ===
write.csv(results, paste0(output_prefix, '.csv'), row.names = FALSE)
cat('\n=== Output Files ===\n')
cat(paste0(output_prefix, '.csv\n'))
cat(paste0(output_prefix, '_pca.pdf\n'))
cat(paste0(output_prefix, '_volcano.pdf\n'))
cat(paste0(output_prefix, '_heatmap.pdf\n'))
