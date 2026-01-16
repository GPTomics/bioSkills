library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('sample_peaks.narrowPeak')

peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = 'org.Hs.eg.db')

peak_anno

plotAnnoPie(peak_anno)
plotAnnoBar(peak_anno)
plotDistToTSS(peak_anno, title = 'Distribution of peaks relative to TSS')

anno_df <- as.data.frame(peak_anno)
write.csv(anno_df, 'annotated_peaks.csv', row.names = FALSE)

promoter_peaks <- anno_df[grep('Promoter', anno_df$annotation), ]
promoter_genes <- unique(promoter_peaks$SYMBOL)
length(promoter_genes)

genes <- unique(anno_df$ENTREZID[!is.na(anno_df$ENTREZID)])
ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05)
head(as.data.frame(ego))
