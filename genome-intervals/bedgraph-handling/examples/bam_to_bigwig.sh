#!/bin/bash
# Reference: deeptools 3.5+, bedtools 2.31+, ucsc-bedgraphtobigwig 445+ | Verify API if version differs
# BAM -> normalized bigWig. Route A (preferred): deepTools bamCoverage in one step.
# Route B: hand-built bedGraph -> bedGraphToBigWig under the strict sort/overlap/chrom.sizes contract.

set -euo pipefail

BAM=$1
CHROM_SIZES=$2          # derive from the EXACT aligned-to FASTA: samtools faidx ref.fa && cut -f1,2 ref.fa.fai
EFFGENOME=${3:-2913022398}   # GRCh38 non-N length (faCount); use the read-length unique-kmer value if MAPQ/uniqueness-filtered

if [ -z "${BAM:-}" ] || [ -z "${CHROM_SIZES:-}" ]; then
    echo "Usage: $0 <bam_file> <chrom.sizes> [effective_genome_size]"
    exit 1
fi

NAME=$(basename "$BAM" .bam)

BIN_SIZE=25             # bp; match to feature width (sharp TF/ATAC 10-25; broad marks 50-200). Compared tracks must share this.

echo "Route A: bamCoverage -> RPGC-normalized bigWig (one step)"
bamCoverage -b "$BAM" -o "${NAME}.bw" \
    --binSize $BIN_SIZE --normalizeUsing RPGC --effectiveGenomeSize "$EFFGENOME" \
    --extendReads --ignoreForNormalization chrX chrM -p 8

echo "Route B (only when a hand-edited text bedGraph is needed): genomecov -> bedGraphToBigWig"
bedtools genomecov -ibam "$BAM" -bga -split > "${NAME}.bedgraph"
# C locale is mandatory: a locale-aware sort triggers "is not case-sensitive sorted"
LC_COLLATE=C sort -k1,1 -k2,2n "${NAME}.bedgraph" > "${NAME}.sorted.bedgraph"
# Inspect the text before converting - the last human-readable checkpoint before opaque binary
head "${NAME}.sorted.bedgraph"
bedGraphToBigWig "${NAME}.sorted.bedgraph" "$CHROM_SIZES" "${NAME}.fromtext.bw"

rm -f "${NAME}.bedgraph" "${NAME}.sorted.bedgraph"

echo "Done: ${NAME}.bw (RPGC) and ${NAME}.fromtext.bw (raw coverage from text route)"
