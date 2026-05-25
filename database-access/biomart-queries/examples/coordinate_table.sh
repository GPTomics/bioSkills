#!/bin/bash
# Reference: R biomaRt 2.58+ (Bioconductor) | Verify API if version differs
# R biomaRt example: version-pinned coordinate + GO term tables.

cat > biomart_query.R <<'EOF'
# Reference: Bioconductor biomaRt 2.58+ | Verify API if version differs
suppressMessages(library(biomaRt))

# Pin to Ensembl release 110 for reproducibility
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 110)

# Coordinate table for chr17 protein-coding genes
coord_df <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name',
                   'chromosome_name', 'start_position', 'end_position',
                   'strand', 'biotype'),
    filters = c('chromosome_name', 'biotype'),
    values = list('17', 'protein_coding'),
    mart = ensembl
)
cat(sprintf('chr17 protein-coding genes: %d\n', nrow(coord_df)))
write.table(coord_df, 'chr17_pc_genes.tsv', sep='\t', row.names=FALSE, quote=FALSE)

# GO annotations (long format: one row per gene-GO pair)
genes_of_interest <- c('TP53', 'BRCA1', 'BRCA2', 'MYC', 'EGFR')
go_df <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name',
                   'go_id', 'name_1006', 'namespace_1003'),
    filters = 'external_gene_name',
    values = genes_of_interest,
    mart = ensembl
)
cat(sprintf('GO annotation rows: %d (long format)\n', nrow(go_df)))
write.table(go_df, 'go_annotations.tsv', sep='\t', row.names=FALSE, quote=FALSE)

# Per-gene GO term count
go_counts <- aggregate(go_id ~ external_gene_name, data=go_df, FUN=length)
print(go_counts)
EOF

echo "=== Running biomart_query.R ==="
Rscript biomart_query.R

echo
echo "=== Outputs ==="
head chr17_pc_genes.tsv
echo
head go_annotations.tsv
