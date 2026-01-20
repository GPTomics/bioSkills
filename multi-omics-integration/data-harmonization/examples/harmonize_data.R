library(DESeq2)
library(sva)

# Simulate multi-omics data
set.seed(42)
n_samples <- 60
samples <- paste0('Sample', 1:n_samples)

# RNA-seq counts
rna_counts <- matrix(rnbinom(n_samples * 1000, size = 10, mu = 100), nrow = 1000, dimnames = list(paste0('Gene', 1:1000), samples))

# Proteomics intensity
protein_int <- matrix(2^rnorm(n_samples * 500, mean = 20, sd = 3), nrow = 500, dimnames = list(paste0('Protein', 1:500), samples))
protein_int[sample(length(protein_int), 0.1 * length(protein_int))] <- NA  # Add missing

# Sample info
sample_info <- data.frame(
    SampleID = samples,
    Condition = rep(c('Control', 'Treatment'), each = n_samples / 2),
    Batch = rep(c('B1', 'B2'), n_samples / 2)
)

# === HARMONIZATION ===

# 1. Normalize RNA-seq
dds <- DESeqDataSetFromMatrix(rna_counts, colData = sample_info, design = ~ 1)
vst_rna <- assay(vst(dds))

# 2. Normalize proteomics
log2_protein <- log2(protein_int)
log2_protein[is.infinite(log2_protein)] <- NA
medians <- apply(log2_protein, 2, median, na.rm = TRUE)
norm_protein <- sweep(log2_protein, 2, medians - median(medians))

# 3. Filter missing values in proteomics
keep_prot <- rowMeans(is.na(norm_protein)) < 0.3
norm_protein <- norm_protein[keep_prot, ]

# 4. Impute remaining missing
for (i in 1:ncol(norm_protein)) {
    nas <- is.na(norm_protein[, i])
    if (any(nas)) {
        q01 <- quantile(norm_protein[, i], 0.01, na.rm = TRUE)
        norm_protein[nas, i] <- rnorm(sum(nas), q01, abs(q01) * 0.1)
    }
}

# 5. Scale
scaled_rna <- t(scale(t(vst_rna)))
scaled_protein <- t(scale(t(norm_protein)))

# Save
cat('Harmonization complete\n')
cat('RNA:', nrow(scaled_rna), 'features x', ncol(scaled_rna), 'samples\n')
cat('Protein:', nrow(scaled_protein), 'features x', ncol(scaled_protein), 'samples\n')

saveRDS(list(RNA = scaled_rna, Protein = scaled_protein, sample_info = sample_info), 'harmonized_data.rds')
