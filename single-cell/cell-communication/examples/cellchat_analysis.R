library(CellChat)
library(Seurat)
library(patchwork)

seurat_obj <- readRDS('seurat_annotated.rds')

cellchat <- createCellChat(object = seurat_obj, group.by = 'cell_type')
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = 'Secreted Signaling')
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = 'triMean')
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

pdf('cellchat_networks.pdf', width = 12, height = 10)
par(mfrow = c(1, 2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = 'Number of interactions')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = 'Interaction strength')
dev.off()

pdf('cellchat_heatmap.pdf', width = 10, height = 8)
netVisual_heatmap(cellchat, color.heatmap = 'Reds')
dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')
pdf('signaling_roles.pdf', width = 10, height = 8)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

top_pathways <- cellchat@netP$pathways[1:5]
pdf('top_pathways.pdf', width = 15, height = 12)
par(mfrow = c(2, 3))
for (pathway in top_pathways) {
    netVisual_aggregate(cellchat, signaling = pathway, layout = 'circle')
}
dev.off()

interactions_df <- subsetCommunication(cellchat)
write.csv(interactions_df, 'cellchat_interactions.csv', row.names = FALSE)
saveRDS(cellchat, 'cellchat_object.rds')
