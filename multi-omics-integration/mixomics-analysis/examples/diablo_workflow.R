library(mixOmics)

# Load breast cancer example data
data(breast.TCGA)
X_mrna <- breast.TCGA$data.train$mrna
X_mirna <- breast.TCGA$data.train$mirna
X_protein <- breast.TCGA$data.train$protein
Y <- breast.TCGA$data.train$subtype

# Prepare blocks
X_blocks <- list(mRNA = X_mrna, miRNA = X_mirna, Protein = X_protein)

# Design matrix (moderate correlation between blocks)
design <- matrix(0.1, ncol = 3, nrow = 3, dimnames = list(names(X_blocks), names(X_blocks)))
diag(design) <- 0

# Run DIABLO
diablo <- block.splsda(X_blocks, Y, ncomp = 2,
                        keepX = list(mRNA = c(10, 10), miRNA = c(10, 10), Protein = c(10, 10)),
                        design = design)

# Performance
cat('DIABLO model trained\n')
cat('Components:', diablo$ncomp, '\n')

# Visualize consensus
pdf('diablo_consensus.pdf', width = 8, height = 6)
plotIndiv(diablo, comp = c(1, 2), blocks = 'consensus', group = Y, legend = TRUE, title = 'DIABLO Consensus')
dev.off()

# Circos plot
pdf('diablo_circos.pdf', width = 10, height = 10)
circosPlot(diablo, cutoff = 0.7, line = TRUE)
dev.off()

# Extract selected features
selected <- lapply(names(X_blocks), function(b) selectVar(diablo, block = b, comp = 1)[[b]]$name)
names(selected) <- names(X_blocks)
cat('Selected features per block:\n')
print(sapply(selected, length))
