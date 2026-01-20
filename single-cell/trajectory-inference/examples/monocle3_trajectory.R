library(Seurat)
library(monocle3)
library(dplyr)

seurat_obj <- readRDS('seurat_clustered.rds')

cds <- as.cell_data_set(seurat_obj)
cds <- cluster_cells(cds, reduction_method = 'UMAP')
cds <- learn_graph(cds)

get_earliest_principal_node <- function(cds, cluster_name) {
    cell_ids <- which(colData(cds)$seurat_clusters == cluster_name)
    closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[cell_ids, ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[['UMAP']])$name[
        as.numeric(names(which.max(table(closest_vertex))))
    ]
    root_pr_nodes
}

root_cluster <- 'HSC'
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, root_cluster))

pdf('trajectory_pseudotime.pdf', width = 10, height = 8)
plot_cells(cds, color_cells_by = 'pseudotime', label_branch_points = TRUE,
           label_leaves = TRUE, label_roots = TRUE)
dev.off()

graph_test_res <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
sig_genes <- graph_test_res %>% filter(q_value < 0.05) %>% arrange(desc(morans_I))
write.csv(sig_genes, 'trajectory_genes.csv')

top_genes <- head(sig_genes$gene_short_name, 10)
pdf('gene_dynamics.pdf', width = 12, height = 8)
plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% top_genes, ],
                         color_cells_by = 'seurat_clusters', min_expr = 0.5)
dev.off()

colData(cds)$pseudotime <- pseudotime(cds)
seurat_obj$pseudotime <- pseudotime(cds)[Cells(seurat_obj)]
saveRDS(seurat_obj, 'seurat_with_pseudotime.rds')
