# Reference: MSnbase 2.28+, ggplot2 3.5+, limma 3.58+ | Verify API if version differs
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

# === 2. LOG2, THEN INSPECT RAW DISTRIBUTIONS (before any normalization) ===
cat('\n=== Raw distribution QC ===\n')
intensities[intensities == 0] <- NA
log2_int <- log2(intensities)

# Median-centering rescales every sample onto a common median, so it mathematically erases the 3x-low
# load that marks a failed injection. A failed sample must be found and dropped HERE, not afterwards.
names(sample_groups) <- colnames(log2_int)   # key by column so a drop cannot misalign the design below
id_counts <- colSums(!is.na(log2_int))
print(data.frame(id_count = id_counts, raw_median_log2 = round(apply(log2_int, 2, median, na.rm = TRUE), 2)))
pdf(paste0(output_prefix, '_raw_boxplot.pdf'), width = 8, height = 5)
boxplot(log2_int, las = 2, main = 'RAW log2 LFQ (pre-normalization)', ylab = 'log2 intensity')
dev.off()

# <50% of the cohort median ID count is a failed injection / low load, not biology.
failed <- names(id_counts)[id_counts < 0.5 * median(id_counts)]
if (length(failed) > 0) {
    cat('Dropping failed samples:', paste(failed, collapse = ', '), '\n')
    log2_int <- log2_int[, !colnames(log2_int) %in% failed, drop = FALSE]
    sample_groups <- sample_groups[colnames(log2_int)]   # keep annotation aligned to surviving columns
}

# === 3. NORMALIZE (only after the raw inspection above) ===
cat('\n=== Normalization ===\n')
sample_medians <- apply(log2_int, 2, median, na.rm = TRUE)
cat('Sample medians before:', round(sample_medians, 2), '\n')
normalized <- sweep(log2_int, 2, sample_medians - median(sample_medians))
cat('Sample medians after:', round(apply(normalized, 2, median, na.rm = TRUE), 2), '\n')

# === 4. FILTER & IMPUTE ===
# NOTE: this is the SIMPLIFIED FALLBACK path. The defensible route (see SKILL.md) filters on
# PER-GROUP completeness (valid in >=50-70% of replicates in >=1 condition, not a blanket rule),
# then MODELS the MNAR dropout with proDA rather than downshift-imputing -- downshift on an on/off
# protein manufactures false positives (the volcano-wing artifact). Downshift is used here only to
# keep the example dependency-light; prefer proDA for real data.
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

# === 5. QC ===
cat('\n=== Quality Control ===\n')
pca <- prcomp(t(imputed), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = rownames(pca$x), Group = sample_groups)
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4) + theme_minimal() +
    labs(x = paste0('PC1 (', var_exp[1], '%)'), y = paste0('PC2 (', var_exp[2], '%)'), title = 'PCA of Protein Abundances')
ggsave(paste0(output_prefix, '_pca.pdf'), p_pca, width = 7, height = 6)

# === 6. DIFFERENTIAL ANALYSIS ===
cat('\n=== Differential Analysis ===\n')
sample_info <- data.frame(sample = colnames(imputed), condition = factor(sample_groups, levels = c('Control', 'Treatment')))
design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(as.matrix(imputed), design)
contrast <- makeContrasts(Treatment - Control, levels = design)
# treat() folds the minimum fold-change INTO the test (a proper hypothesis against |lfc| > threshold),
# instead of a post-hoc logFC AND adj.P double filter -- the double filter is a collider/selection
# effect whose realized FDR can exceed the nominal rate (SKILL.md). Significance is then adj.P alone.
fit2 <- treat(contrasts.fit(fit, contrast), lfc = lfc_threshold, trend = TRUE, robust = TRUE)

results <- topTreat(fit2, coef = 1, number = Inf, adjust.method = 'BH')
results$protein <- rownames(results)
results$significant <- results$adj.P.Val < fdr_threshold

cat('Total proteins tested:', nrow(results), '\n')
cat('Significant:', sum(results$significant), '\n')
cat('  Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('  Down-regulated:', sum(results$significant & results$logFC < 0), '\n')

# === 7. VISUALIZATION ===
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

# === 8. EXPORT ===
write.csv(results, paste0(output_prefix, '.csv'), row.names = FALSE)
cat('\n=== Output Files ===\n')
cat(paste0(output_prefix, '.csv\n'))
cat(paste0(output_prefix, '_raw_boxplot.pdf\n'))
cat(paste0(output_prefix, '_pca.pdf\n'))
cat(paste0(output_prefix, '_volcano.pdf\n'))
cat(paste0(output_prefix, '_heatmap.pdf\n'))
