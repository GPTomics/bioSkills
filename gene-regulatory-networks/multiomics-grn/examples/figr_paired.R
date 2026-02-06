# FigR: TF-gene regulatory inference from paired scRNA+scATAC

library(FigR)
library(Seurat)
library(Signac)

seurat_obj <- readRDS('multiome_seurat.rds')
cat('Loaded', ncol(seurat_obj), 'cells\n')

# Step 1: peak-gene correlations
# p.cut 0.05: significance threshold for peak-gene links
peak_gene_cors <- runGenePeakcorr(
    ATAC.se = seurat_obj@assays$ATAC,
    RNAmat = seurat_obj@assays$RNA@data,
    genome = 'hg38',
    nCores = 8,
    p.cut = 0.05
)
cat('Peak-gene correlations:', nrow(peak_gene_cors), '\n')

# Step 2: DORC (domains of regulatory chromatin) scores
dorc_scores <- getDORCScores(seurat_obj@assays$ATAC, peak_gene_cors)

# Step 3: TF-gene regulation scores
fig_results <- runFigR(
    ATAC.se = seurat_obj@assays$ATAC,
    dorcTab = peak_gene_cors,
    genome = 'hg38',
    dorcMat = dorc_scores,
    rnaMat = seurat_obj@assays$RNA@data,
    nCores = 8
)

top_links <- fig_results[order(abs(fig_results$Score), decreasing = TRUE), ]
cat('Top regulatory TF-gene links:\n')
print(head(top_links, 20))

write.csv(fig_results, 'figr_results.csv', row.names = FALSE)
cat('Saved figr_results.csv\n')
