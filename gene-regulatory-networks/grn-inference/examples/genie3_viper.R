# GENIE3 edge inference + VIPER TF-activity (the "activity, not edges" workflow)
# Reference: GENIE3 1.24+, viper 1.36+ (Bioconductor) | Verify API if version differs
#
# Demonstrates: tree-ensemble network inference, regulon assembly with Mode of Regulation,
# msVIPER master-regulator analysis, and a per-sample VIPER activity matrix.

library(GENIE3)
library(viper)
library(Biobase)

# GENIE3 convention: genes in ROWS, samples in COLUMNS (transpose of the WGCNA layout).
expr <- as.matrix(read.csv('normalized_counts.csv', row.names = 1))
regulators <- intersect(readLines('tf_list.txt'), rownames(expr))
cat('Inferring network over', length(regulators), 'TFs x', nrow(expr), 'genes\n')

# Tree ensembles are stochastic; seed for reproducibility. Restricting predictors to TFs
# is what orients the edges -- this is an assumption, not an inferred causal direction.
set.seed(42)
weight_matrix <- GENIE3(expr, regulators = regulators, treeMethod = 'RF',
                        K = 'sqrt', nTrees = 1000, nCores = 4)
link_list <- getLinkList(weight_matrix)        # ranked TF-target edges (not thresholded)
write.csv(link_list, 'genie3_links.csv', row.names = FALSE)

# Assemble a regulon: keep the top targets PER TF (not a global quantile) so regulons clear
# VIPER's minsize = 25, then assign a Mode-of-Regulation sign from the TF-target correlation.
n_top <- 50
regulon <- lapply(split(link_list, as.character(link_list$regulatoryGene)), function(df) {
    df <- df[order(df$weight, decreasing = TRUE), ]
    df <- head(df, n_top)
    tf <- as.character(df$regulatoryGene[1])
    targets <- as.character(df$targetGene)
    tfmode <- sign(sapply(targets, function(g) cor(expr[g, ], expr[tf, ])))
    list(tfmode = setNames(tfmode, targets),
         likelihood = setNames(df$weight / max(df$weight), targets))
})
regulon <- regulon[sapply(regulon, function(r) length(r$tfmode) >= 25)]   # VIPER minsize
class(regulon) <- 'regulon'

# Two-group signature -> msVIPER master regulators. A non-DE TF can still top the ranking
# because activity is read from the coordinated shift of its targets, not its own mRNA.
pheno <- read.csv('sample_info.csv', row.names = 1)
eset <- ExpressionSet(assayData = expr,
                      phenoData = AnnotatedDataFrame(pheno))
signature <- rowTtest(eset, pheno = 'condition', group1 = 'tumor', group2 = 'normal')
sig_z <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
nullmodel <- ttestNull(eset, pheno = 'condition', group1 = 'tumor', group2 = 'normal', per = 1000)

mra <- msviper(sig_z, regulon, nullmodel)
cat('\nTop master regulators (msVIPER NES):\n')
print(head(summary(mra), 15))

# Per-sample TF-activity matrix for stratification/clustering.
activity <- viper(eset, regulon, method = 'scale')
write.csv(exprs(activity), 'viper_activity.csv')
cat('\nWrote per-sample activity for', nrow(activity), 'regulators\n')
