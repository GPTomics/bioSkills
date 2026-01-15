# edgeR with batch effect correction

library(edgeR)

# Simulate data with batch effects
set.seed(42)
n_genes <- 1000
n_samples <- 12

counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('gene', 1:n_genes),
                                paste0('sample', 1:n_samples)))

sample_info <- data.frame(
    condition = factor(rep(c('control', 'treated'), 6)),
    batch = factor(rep(c('batch1', 'batch2', 'batch3'), each = 4)),
    row.names = colnames(counts)
)

# Create DGEList
y <- DGEList(counts = counts)

# Filter and normalize
keep <- filterByExpr(y, group = sample_info$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Design with batch correction
design <- model.matrix(~ batch + condition, data = sample_info)
cat('Design matrix columns:\n')
print(colnames(design))

# Estimate dispersion and fit
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Test condition effect (last coefficient, controlling for batch)
qlf <- glmQLFTest(fit, coef = 'conditiontreated')

cat('\nTop DE genes (condition effect, controlling for batch):\n')
print(topTags(qlf, n = 10))

cat('\nSummary:\n')
print(summary(decideTests(qlf)))
