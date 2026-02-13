# Reference: vegan 2.6+ | Verify API if version differs
library(LEA)
library(qvalue)

# --- Step 1: Convert VCF to LEA formats ---
vcf2lfmm('variants.vcf', 'genotypes.lfmm')
vcf2geno('variants.vcf', 'genotypes.geno')

# --- Step 2: Estimate K with sNMF ---
# K=1:10: test 1 to 10 ancestral populations
# repetitions=5: run 5 replicates per K for stability
# entropy=TRUE: compute cross-entropy for model selection
snmf_result <- snmf('genotypes.geno', K = 1:10, repetitions = 5,
                     entropy = TRUE, project = 'new')

# Cross-entropy plot: lower = better fit; look for elbow
pdf('snmf_cross_entropy.pdf', width = 7, height = 5)
plot(snmf_result, col = 'blue', pch = 19, cex = 1.2,
     xlab = 'K (number of ancestral populations)',
     ylab = 'Cross-entropy criterion')
dev.off()

# Select best K (example: K=3)
best_K <- 3
cat('Selected K:', best_K, '\n')

# --- Step 3: Prepare environmental data ---
# env_vars: matrix with rows = individuals, columns = environmental variables
# Must match order of individuals in genotype file
env_vars <- as.matrix(read.table('environment.env', header = TRUE))
cat('Environmental variables:', ncol(env_vars), '\n')
cat('Individuals:', nrow(env_vars), '\n')

# --- Step 4: Run LFMM2 ---
genotypes <- read.lfmm('genotypes.lfmm')

# K: number of latent factors (= estimated K from sNMF)
# Latent factors capture population structure without explicit assignment
lfmm_result <- lfmm2(input = genotypes, env = env_vars, K = best_K)

# --- Step 5: Test associations and calibrate ---
pvalues <- lfmm2.test(lfmm_result, input = genotypes, env = env_vars,
                       full = TRUE, genomic.control = TRUE)

# Genomic inflation factor (GIF) per environmental variable
# lambda ~ 1.0: well-calibrated; >1.5: increase K; <0.5: decrease K
for (i in 1:ncol(env_vars)) {
    gif <- median(qchisq(1 - pvalues$pvalues[, i], df = 1)) / qchisq(0.5, df = 1)
    cat(sprintf('Variable %d GIF (lambda): %.3f\n', i, gif))
}

# --- Step 6: Multiple testing correction ---
results_list <- list()
for (i in 1:ncol(env_vars)) {
    pvals_i <- pvalues$pvalues[, i]
    # Replace any NA p-values with 1 (non-significant)
    pvals_i[is.na(pvals_i)] <- 1

    # q-value < 0.05: Storey FDR, more powerful than Benjamini-Hochberg
    qvals_i <- qvalue(pvals_i)$qvalues

    candidates_i <- which(qvals_i < 0.05)
    cat(sprintf('Variable %d: %d candidate loci (q < 0.05)\n', i, length(candidates_i)))
    results_list[[i]] <- data.frame(locus = candidates_i, pvalue = pvals_i[candidates_i],
                                     qvalue = qvals_i[candidates_i], env_var = i)
}

all_candidates <- do.call(rbind, results_list)
cat('\nTotal candidate associations:', nrow(all_candidates), '\n')
cat('Unique candidate loci:', length(unique(all_candidates$locus)), '\n')

# --- Step 7: Manhattan plot ---
pdf('lfmm_manhattan.pdf', width = 12, height = 5)
par(mfrow = c(1, ncol(env_vars)))
for (i in 1:ncol(env_vars)) {
    pvals_i <- pvalues$pvalues[, i]
    pvals_i[is.na(pvals_i)] <- 1
    plot(-log10(pvals_i), pch = 19, cex = 0.3, col = 'grey40',
         xlab = 'Locus index', ylab = '-log10(p-value)',
         main = paste('Environmental variable', i))
    # Bonferroni threshold for reference (conservative)
    abline(h = -log10(0.05 / length(pvals_i)), col = 'red', lty = 2)

    qvals_i <- qvalue(pvals_i)$qvalues
    sig <- which(qvals_i < 0.05)
    if (length(sig) > 0) {
        points(sig, -log10(pvals_i[sig]), pch = 19, cex = 0.5, col = 'red')
    }
}
dev.off()

# --- Step 8: Export candidate loci ---
write.csv(all_candidates, 'lfmm_candidates.csv', row.names = FALSE)
cat('Results written to lfmm_candidates.csv\n')
