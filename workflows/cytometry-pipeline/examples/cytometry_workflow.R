library(CATALYST)
library(diffcyt)
library(SingleCellExperiment)
library(flowCore)
library(ggplot2)

# === CONFIGURATION ===
data_dir <- 'data/'
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. DEFINE PANEL ===
panel <- data.frame(
    fcs_colname = c('FSC-A', 'SSC-A', 'CD45', 'CD3', 'CD4', 'CD8', 'CD19', 'CD14', 'Ki67'),
    antigen = c('FSC', 'SSC', 'CD45', 'CD3', 'CD4', 'CD8', 'CD19', 'CD14', 'Ki67'),
    marker_class = c('none', 'none', 'type', 'type', 'type', 'type', 'type', 'type', 'state')
)

# === 2. LOAD METADATA ===
fcs_files <- list.files(data_dir, pattern = '\\.fcs$')
md <- data.frame(
    file_name = fcs_files,
    sample_id = gsub('\\.fcs$', '', fcs_files),
    condition = ifelse(grepl('ctrl|control', fcs_files, ignore.case = TRUE), 'Control', 'Treatment')
)
cat('Found', nrow(md), 'FCS files\n')

# === 3. LOAD AND PREPARE DATA ===
cat('Loading data...\n')
fcs_paths <- file.path(data_dir, md$file_name)
fs <- read.flowSet(fcs_paths)
sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 150, FACS = TRUE)
cat('Loaded', ncol(sce), 'cells\n')

# === 4. CLUSTERING ===
cat('Clustering...\n')
sce <- cluster(sce, features = 'type', xdim = 10, ydim = 10, maxK = 20, seed = 42)

# === 5. UMAP ===
cat('Running UMAP...\n')
sce <- runDR(sce, dr = 'UMAP', features = 'type')

# === 6. PLOTS ===
cat('Generating plots...\n')
plotDR(sce, dr = 'UMAP', color_by = 'meta20')
ggsave(file.path(output_dir, 'umap_clusters.png'), width = 8, height = 6)

plotExprHeatmap(sce, features = 'type', k = 'meta20', by = 'cluster_id', scale = 'last')
ggsave(file.path(output_dir, 'heatmap.png'), width = 12, height = 8)

# === 7. DIFFERENTIAL ANALYSIS ===
cat('Running differential analysis...\n')
design <- createDesignMatrix(ei(sce), cols_design = 'condition')
contrast <- createContrast(c(0, 1))

res_DA <- testDA_edgeR(sce, design, contrast, cluster_id = 'meta20')
da_results <- as.data.frame(rowData(res_DA))
da_results <- da_results[order(da_results$p_adj), ]

cat('\nSignificant clusters (FDR < 0.05):\n')
print(da_results[da_results$p_adj < 0.05, c('cluster_id', 'logFC', 'p_adj')])

# === 8. VISUALIZATION ===
plotDiffHeatmap(sce, res_DA, all = TRUE, fdr = 0.05)
ggsave(file.path(output_dir, 'da_heatmap.png'), width = 10, height = 8)

plotAbundances(sce, k = 'meta20', by = 'cluster_id', group_by = 'condition')
ggsave(file.path(output_dir, 'abundances.png'), width = 12, height = 8)

# === 9. SAVE RESULTS ===
write.csv(da_results, file.path(output_dir, 'da_results.csv'), row.names = FALSE)
saveRDS(sce, file.path(output_dir, 'cytometry_analysis.rds'))
cat('Analysis complete! Results saved to', output_dir, '\n')
