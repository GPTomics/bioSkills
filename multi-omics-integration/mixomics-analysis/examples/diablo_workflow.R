# Reference: mixOmics 6.26+ | Verify API if version differs
# DIABLO supervised multi-block signature on the shipped breast.TCGA demo. The discipline:
# choose the design from the goal, tune ncomp then keepX with BER inside CV folds, and treat
# the selected features as candidates whose accuracy needs an external/held-out estimate.

library(mixOmics)

data(breast.TCGA)
X_blocks <- list(mRNA=breast.TCGA$data.train$mrna,
                 miRNA=breast.TCGA$data.train$mirna,
                 protein=breast.TCGA$data.train$protein)   # samples x features, matched rows
Y <- breast.TCGA$data.train$subtype
stopifnot(identical(rownames(X_blocks$mRNA), rownames(X_blocks$protein)))   # matched rownames or the result is garbage

design <- matrix(0.5, nrow=length(X_blocks), ncol=length(X_blocks),
                 dimnames=list(names(X_blocks), names(X_blocks)))   # 0.5-1 favors cross-block correlation; <0.5 favors prediction
diag(design) <- 0

ncomp_perf <- perf(block.plsda(X_blocks, Y, ncomp=4, design=design),
                   validation='Mfold', folds=10, nrepeat=10)        # tune ncomp on a NON-sparse model; nrepeat>=10
ncomp <- ncomp_perf$choice.ncomp$WeightedVote['Overall.BER', 'max.dist']   # read the elbow; here it lands at 2

tune <- tune.block.splsda(X_blocks, Y, ncomp=ncomp, design=design,
                          test.keepX=list(mRNA=c(8, 16), miRNA=c(8, 16), protein=c(8, 16)),
                          validation='Mfold', folds=10, nrepeat=10,
                          measure='BER', BPPARAM=BiocParallel::MulticoreParam(workers=2))   # BER for imbalanced subtypes; selection happens inside folds (cpus= is defunct)

diablo <- block.splsda(X_blocks, Y, ncomp=ncomp, keepX=tune$choice.keepX, design=design)
sel_mrna <- selectVar(diablo, block='mRNA', comp=1)$mRNA$name        # candidate features, not validated biomarkers
print(head(sel_mrna))

# Honest accuracy comes from the held-out block, never touched during tuning:
pred <- predict(diablo, newdata=list(mRNA=breast.TCGA$data.test$mrna, miRNA=breast.TCGA$data.test$mirna))
