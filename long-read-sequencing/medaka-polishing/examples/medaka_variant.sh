#!/bin/bash
# Reference: medaka 2.2+, bcftools 1.19+ | Verify API if version differs
# Haploid consensus/variant calling for microbial, mitochondrial, or viral ONT samples.
# In medaka v2, medaka_variant is the renamed haploid wrapper (was medaka_haploid_variant).
# For DIPLOID/germline ONT calling, use Clair3 instead (medaka diploid was deprecated).
set -euo pipefail

READS=${1:?Usage: $0 <reads.fq.gz> <reference.fa> [output_dir] [threads]}
REFERENCE=${2:?reference required}
OUTPUT_DIR=${3:-medaka_variants}
THREADS=${4:-8}

# Model auto-detected from the basecaller annotation; supply -m only if it cannot resolve.
medaka_variant -i "$READS" -r "$REFERENCE" -o "$OUTPUT_DIR" -t "$THREADS"

VCF="${OUTPUT_DIR}/medaka.annotated.vcf"
[ -f "$VCF" ] || { echo 'Error: variant calling failed'; exit 1; }

echo 'Variant summary:'
bcftools stats "$VCF" | grep '^SN'

# QUAL>20 is a permissive starting filter for haploid consensus variants, not a hard rule.
bcftools filter -i 'QUAL>20' "$VCF" > "${OUTPUT_DIR}/medaka.filtered.vcf"
echo "Filtered (QUAL>20): ${OUTPUT_DIR}/medaka.filtered.vcf"
