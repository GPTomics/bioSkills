# Reference: DESeq2 1.42+, Salmon 1.10+, edgeR 4.0+, kallisto 0.50+, scanpy 1.10+ | Verify API if version differs
library(GenomicFeatures)

gtf_file <- 'Homo_sapiens.GRCh38.110.gtf.gz'

cat('Creating TxDb from GTF...\n')
txdb <- makeTxDbFromGFF(gtf_file)

cat('Extracting transcript-gene mapping...\n')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- select(txdb, k, 'GENEID', 'TXNAME')
tx2gene <- tx2gene[, c('TXNAME', 'GENEID')]

cat('Saving tx2gene.csv...\n')
write.csv(tx2gene, 'tx2gene.csv', row.names = FALSE)

cat('Done! Created tx2gene.csv with', nrow(tx2gene), 'transcripts\n')
