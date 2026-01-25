# Batch Design Examples
# Demonstrates experimental design to minimize batch effects

library(sva)

# =============================================================================
# Creating a Balanced Design
# =============================================================================

# Example: 24 samples, 2 conditions, 3 batches
set.seed(42)

# BAD DESIGN: Condition confounded with batch
bad_design <- data.frame(
  sample = paste0('S', 1:24),
  condition = c(rep('Control', 12), rep('Treated', 12)),
  batch = c(rep('Batch1', 12), rep('Batch2', 12))  # All controls in batch 1!
)

cat('=== BAD DESIGN (Confounded) ===\n')
print(table(bad_design$condition, bad_design$batch))
cat('Problem: Cannot distinguish batch effect from treatment effect!\n\n')

# GOOD DESIGN: Balanced across batches
samples <- data.frame(
  sample = paste0('S', 1:24),
  condition = rep(c('Control', 'Treated'), each = 12)
)

# Randomly assign to 3 batches, balanced by condition
samples$batch <- NA
for (cond in c('Control', 'Treated')) {
  idx <- which(samples$condition == cond)
  samples$batch[idx] <- sample(rep(paste0('Batch', 1:3), length.out = length(idx)))
}

cat('=== GOOD DESIGN (Balanced) ===\n')
print(table(samples$condition, samples$batch))
cat('Each batch has both conditions - batch effect is correctable\n\n')

# =============================================================================
# Detecting Batch Effects with SVA
# =============================================================================

# Simulate count data with batch effect
n_genes <- 1000
n_samples <- 24

# Base expression
counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes, ncol = n_samples)
rownames(counts) <- paste0('Gene', 1:n_genes)
colnames(counts) <- samples$sample

# Add batch effect (multiply batch 2 by 1.5, batch 3 by 0.7)
batch_effects <- c(Batch1 = 1, Batch2 = 1.5, Batch3 = 0.7)
for (i in 1:n_samples) {
  counts[, i] <- counts[, i] * batch_effects[samples$batch[i]]
}

# Add treatment effect to first 50 genes
treatment_genes <- 1:50
counts[treatment_genes, samples$condition == 'Treated'] <-
  counts[treatment_genes, samples$condition == 'Treated'] * 2

# Normalize (simple log2)
counts_norm <- log2(counts + 1)

# Estimate number of surrogate variables
mod <- model.matrix(~condition, data = samples)
mod0 <- model.matrix(~1, data = samples)

n_sv <- num.sv(counts_norm, mod, method = 'leek')
cat('Estimated number of hidden factors (SVs):', n_sv, '\n')

# =============================================================================
# PCA to Visualize Batch Effects
# =============================================================================

pca <- prcomp(t(counts_norm), scale. = TRUE)

cat('\n=== PCA Results ===\n')
cat('PC1 variance explained:', round(summary(pca)$importance[2,1] * 100, 1), '%\n')
cat('PC2 variance explained:', round(summary(pca)$importance[2,2] * 100, 1), '%\n')

# Check if PCs correlate with batch
pc1_batch_cor <- cor(pca$x[,1], as.numeric(factor(samples$batch)))
pc1_cond_cor <- cor(pca$x[,1], as.numeric(factor(samples$condition)))

cat('\nPC1 correlation with batch:', round(pc1_batch_cor, 3), '\n')
cat('PC1 correlation with condition:', round(pc1_cond_cor, 3), '\n')

if (abs(pc1_batch_cor) > 0.5) {
  cat('WARNING: Strong batch effect detected in PC1!\n')
}

# =============================================================================
# Batch Correction with ComBat
# =============================================================================

library(sva)

# ComBat correction (preserving condition effect)
combat_counts <- ComBat(dat = counts_norm,
                        batch = samples$batch,
                        mod = mod)

# Check PCA after correction
pca_corrected <- prcomp(t(combat_counts), scale. = TRUE)
pc1_batch_cor_after <- cor(pca_corrected$x[,1], as.numeric(factor(samples$batch)))
pc1_cond_cor_after <- cor(pca_corrected$x[,1], as.numeric(factor(samples$condition)))

cat('\n=== After ComBat Correction ===\n')
cat('PC1 correlation with batch:', round(pc1_batch_cor_after, 3), '\n')
cat('PC1 correlation with condition:', round(pc1_cond_cor_after, 3), '\n')

# =============================================================================
# Design Checklist
# =============================================================================

cat('\n=== Experimental Design Checklist ===\n')
cat('Before experiment:\n')
cat('  [ ] Balance conditions across batches\n')
cat('  [ ] Balance known covariates (sex, age) across batches\n')
cat('  [ ] Randomize sample order within batches\n')
cat('  [ ] Plan for technical replicates across batches\n')
cat('  [ ] Document batch assignments\n\n')

cat('During experiment:\n')
cat('  [ ] Record date, operator, reagent lots\n')
cat('  [ ] Note any protocol deviations\n')
cat('  [ ] Include positive/negative controls\n\n')

cat('During analysis:\n')
cat('  [ ] Check PCA for batch clustering\n')
cat('  [ ] Test batch-condition correlation\n')
cat('  [ ] Apply appropriate correction if needed\n')
cat('  [ ] Include batch in statistical model\n')
