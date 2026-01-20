library(CATALYST)
library(SingleCellExperiment)

# Panel definition
panel <- data.frame(
    fcs_colname = c('CD45', 'CD3', 'CD4', 'CD8', 'CD20', 'CD14'),
    antigen = c('CD45', 'CD3', 'CD4', 'CD8', 'CD20', 'CD14'),
    marker_class = rep('type', 6)
)

# Sample metadata
md <- data.frame(
    file_name = list.files('data', pattern = '\\.fcs$'),
    sample_id = paste0('S', 1:4),
    condition = c('Control', 'Control', 'Treatment', 'Treatment')
)

# Load and prepare data
fs <- read.flowSet(file.path('data', md$file_name))
sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 5)

cat('Loaded', ncol(sce), 'cells\n')

# Clustering
sce <- cluster(sce, features = 'type', xdim = 10, ydim = 10, maxK = 20, seed = 42)
cat('Clustering complete\n')

# UMAP
sce <- runDR(sce, dr = 'UMAP', features = 'type')

# Summary
cat('\nCluster sizes (meta20):\n')
print(table(cluster_ids(sce, 'meta20')))

# Save
saveRDS(sce, 'sce_clustered.rds')
cat('\nSaved to sce_clustered.rds\n')
