# Basic edgeR workflow for differential expression analysis

library(edgeR)

# Simulate example data
set.seed(42)
n_genes <- 1000
n_samples <- 6

counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('gene', 1:n_genes),
                                paste0('sample', 1:n_samples)))

group <- factor(rep(c('control', 'treated'), each = 3))

# Create DGEList
y <- DGEList(counts = counts, group = group)

# Filter low-expression genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
cat('Genes after filtering:', nrow(y), '\n')

# Normalize
y <- calcNormFactors(y)
cat('Normalization factors:', round(y$samples$norm.factors, 3), '\n')

# Create design matrix
design <- model.matrix(~ group)

# Estimate dispersion
y <- estimateDisp(y, design)
cat('Common dispersion:', round(y$common.dispersion, 4), '\n')

# Fit quasi-likelihood model
fit <- glmQLFit(y, design)

# Test for differential expression
qlf <- glmQLFTest(fit, coef = 2)

# View top genes
cat('\nTop 10 differentially expressed genes:\n')
print(topTags(qlf, n = 10))

# Summary of DE genes
cat('\nSummary of DE genes at FDR < 0.05:\n')
print(summary(decideTests(qlf)))
