# Reference: FigR 0.1.0+, FNN 1.1+, Signac 1.12+, Seurat 5.0+, SummarizedExperiment 1.30+ | Verify API if version differs
# FigR: DORC-based TF-gene regulatory inference from paired scRNA+scATAC

library(FigR)
library(FNN)
library(Seurat)
library(Signac)
library(SummarizedExperiment)

seurat_obj <- readRDS('multiome_seurat.rds')
cat('Loaded', ncol(seurat_obj), 'cells\n')

# FigR expects ATAC as a RangedSummarizedExperiment (peak counts + peak GRanges),
# not a Signac ChromatinAssay -- construct it from the assay counts and peak ranges.
atac_counts <- GetAssayData(seurat_obj, assay = 'ATAC', slot = 'counts')
atac_se <- SummarizedExperiment(assays = list(counts = atac_counts),
                                rowRanges = granges(seurat_obj[['ATAC']]))
rna_mat <- GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')

# Cell kNN graph for smoothing (built from an LSI/integrated reduction; drop LSI comp 1).
lsi <- Embeddings(seurat_obj, 'lsi')[, 2:30]
cellkNN <- get.knn(lsi, k = 30)$nn.index
rownames(cellkNN) <- colnames(seurat_obj)

# Step 1: peak-gene correlations, then filter to significant links (pvalZ <= 0.05).
cisCor <- runGenePeakcorr(ATAC.se = atac_se, RNAmat = rna_mat, genome = 'hg38', nCores = 8)
cisCor.filt <- cisCor[cisCor$pvalZ <= 0.05, ]

# Step 2: call DORCs (genes with >= cutoff significant peaks), then score and SMOOTH.
# scATAC sparsity makes raw DORC/RNA scores meaningless -- smoothing over the kNN is required.
dorcGenes <- dorcJPlot(cisCor.filt, cutoff = 10, returnGeneList = TRUE)
dorcMat <- getDORCScores(atac_se, cisCor.filt, geneList = dorcGenes, nCores = 8)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN, mat = dorcMat, nCores = 8)
rnaMat.s  <- smoothScoresNN(NNmat = cellkNN, mat = rna_mat, nCores = 8)

# Step 3: TF-DORC regulation scores (signed: activator/repressor) on the smoothed matrices.
figR <- runFigRGRN(ATAC.se = atac_se, dorcTab = cisCor.filt, dorcMat = dorcMat.s,
                   rnaMat = rnaMat.s, genome = 'hg38', nCores = 8)

top_links <- figR[order(abs(figR$Score), decreasing = TRUE), ]
cat('Top regulatory TF-gene links:\n')
print(head(top_links, 20))

write.csv(figR, 'figr_results.csv', row.names = FALSE)
cat('Saved figr_results.csv\n')
