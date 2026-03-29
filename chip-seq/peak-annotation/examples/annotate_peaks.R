# Reference: ChIPseeker 1.38+, GenomicFeatures 1.54+, rtracklayer 1.62+ | Verify API if version differs
library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)

# --- Custom GTF approach (use when a specific GTF is provided) ---
txdb <- makeTxDbFromGFF('genes.gtf.gz', format = 'gtf')
peaks <- readPeakFile('peaks.bed')

# overlap='all': host-gene convention -- peaks inside a gene body are assigned to
# that gene rather than the nearest TSS gene. Use overlap='TSS' (default) to match
# HOMER behavior where gene assignment is always by nearest TSS.
peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-2000, 2000), overlap = 'all')
anno_df <- as.data.frame(peak_anno)

# Map gene symbols from GTF (annoDb does not work with custom TxDb)
gtf <- import('genes.gtf.gz')
gene_map <- unique(data.frame(
    gene_id = sub('\\..*', '', gtf$gene_id),
    symbol = gtf$gene_name, stringsAsFactors = FALSE))
gene_map <- gene_map[!is.na(gene_map$symbol), ]
anno_df$geneId_base <- sub('\\..*', '', anno_df$geneId)
anno_df$SYMBOL <- gene_map$symbol[match(anno_df$geneId_base, gene_map$gene_id)]

# Collapse annotation categories: promoter, exon, intron, intergenic
collapse_annotation <- function(ann) {
    ifelse(grepl('Promoter', ann), 'promoter',
    ifelse(grepl("5' UTR|3' UTR|Exon", ann), 'exon',
    ifelse(grepl('Intron', ann), 'intron', 'intergenic')))
}
anno_df$feature <- collapse_annotation(anno_df$annotation)

# Export as TSV
output <- data.frame(chr = anno_df$seqnames, start = anno_df$start, end = anno_df$end,
    nearest_gene = anno_df$SYMBOL, distance_to_tss = anno_df$distanceToTSS,
    feature = anno_df$feature)
write.table(output, 'annotations.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

# Visualization
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno, title = 'Distribution of peaks relative to TSS')

# --- Standard TxDb approach (when no custom GTF is needed) ---
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = 'org.Hs.eg.db')
