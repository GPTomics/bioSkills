# edgeR with multiple contrasts

library(edgeR)

# Simulate data with three conditions
set.seed(42)
n_genes <- 1000
n_samples <- 9

counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('gene', 1:n_genes),
                                paste0('sample', 1:n_samples)))

group <- factor(rep(c('control', 'drugA', 'drugB'), each = 3))

# Create DGEList and process
y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Design matrix without intercept (for easy contrasts)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Estimate dispersion and fit
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Define contrasts
contrasts <- makeContrasts(
    DrugA_vs_Control = drugA - control,
    DrugB_vs_Control = drugB - control,
    DrugA_vs_DrugB = drugA - drugB,
    levels = design
)

# Test each contrast
cat('Drug A vs Control:\n')
qlf_A <- glmQLFTest(fit, contrast = contrasts[, 'DrugA_vs_Control'])
print(topTags(qlf_A, n = 5))

cat('\nDrug B vs Control:\n')
qlf_B <- glmQLFTest(fit, contrast = contrasts[, 'DrugB_vs_Control'])
print(topTags(qlf_B, n = 5))

cat('\nDrug A vs Drug B:\n')
qlf_AB <- glmQLFTest(fit, contrast = contrasts[, 'DrugA_vs_DrugB'])
print(topTags(qlf_AB, n = 5))
