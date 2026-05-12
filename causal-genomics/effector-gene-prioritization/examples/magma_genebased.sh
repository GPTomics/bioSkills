#!/usr/bin/env bash
# Reference: MAGMA 1.10+, PLINK 1.9+, 1000 Genomes EUR reference | Verify CLI flags if version differs
#
# MAGMA gene-based and gene-set analysis from GWAS summary statistics.
# Three stages:
#   1. SNP-to-gene annotation (window-based; FUMA default 35kb upstream + 10kb downstream)
#   2. Gene-based association (multi-SNP joint test per gene; LD-corrected via reference panel)
#   3. Gene-set enrichment (competitive test against MSigDB C2)

set -euo pipefail

GWAS_PVAL='gwas.pval.tsv'
GWAS_NCOL='N'
GENE_LOC='NCBI37.3.gene.loc'
SNP_LOC='gwas.snploc'
REF_BFILE='g1000_eur'
GENESET_GMT='msigdb_v7.5_C2.gmt'
OUT_PREFIX='magma_run'

# Window: 35 kb upstream + 10 kb downstream (FUMA recommendation; balances regulatory capture vs gene-dense dilution)
WINDOW='35,10'

# Step 1: SNP-to-gene annotation
magma --annotate window="$WINDOW" \
    --snp-loc "$SNP_LOC" \
    --gene-loc "$GENE_LOC" \
    --out "${OUT_PREFIX}_annot"

# Step 2: Gene-based association
magma --bfile "$REF_BFILE" \
    --pval "$GWAS_PVAL" "ncol=${GWAS_NCOL}" \
    --gene-annot "${OUT_PREFIX}_annot.genes.annot" \
    --out "${OUT_PREFIX}_gene"

# Step 3: Gene-set enrichment (competitive test; preferred over self-contained)
magma --gene-results "${OUT_PREFIX}_gene.genes.raw" \
    --set-annot "$GENESET_GMT" \
    --out "${OUT_PREFIX}_geneset"

# Outputs:
#   ${OUT_PREFIX}_gene.genes.out: per-gene Z, p, n_SNPs
#   ${OUT_PREFIX}_geneset.gsa.out: per-set p, beta, beta_se

# Bonferroni threshold for ~20k autosomal protein-coding genes
N_GENES=$(wc -l < "${OUT_PREFIX}_gene.genes.out")
echo "Bonferroni threshold at alpha=0.05 across ${N_GENES} genes: $(echo "0.05 / ${N_GENES}" | bc -l)"

# Top 50 genes by p
# MAGMA .genes.out columns: GENE CHR START STOP NSNPS NPARAM N ZSTAT P
# P is the 9th whitespace-separated field; sort ascending numerically.
{ head -1 "${OUT_PREFIX}_gene.genes.out"
  tail -n +2 "${OUT_PREFIX}_gene.genes.out" | sort -k9,9g
} | head -51 > "${OUT_PREFIX}_top50.tsv"

echo "MAGMA gene-based and gene-set complete."
echo "Pair with PoPS (FinucaneLab/pops) using ${OUT_PREFIX}_gene.genes.raw as --magma_prefix input."
