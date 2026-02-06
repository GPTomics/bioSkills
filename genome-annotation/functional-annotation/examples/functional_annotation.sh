#!/bin/bash
# Functional annotation with eggNOG-mapper and InterProScan
set -euo pipefail

PROTEINS=$1
EGGNOG_DB=${2:-/path/to/eggnog_db}
OUTDIR=${3:-functional_out}
THREADS=${4:-16}

mkdir -p $OUTDIR

PROT_COUNT=$(grep -c '^>' $PROTEINS)
echo "=== Functional Annotation Pipeline ==="
echo "Input proteins: $PROT_COUNT"

# eggNOG-mapper
echo ""
echo "Running eggNOG-mapper..."
emapper.py \
    -i $PROTEINS \
    --output annot \
    --output_dir ${OUTDIR}/eggnog \
    --data_dir $EGGNOG_DB \
    --cpu $THREADS \
    -m diamond \
    --go_evidence non-electronic \
    --override

EGGNOG_HITS=$(grep -v '^#' ${OUTDIR}/eggnog/annot.emapper.annotations | grep -c -v '^$' || echo 0)
echo "eggNOG hits: $EGGNOG_HITS / $PROT_COUNT"

# InterProScan
echo ""
echo "Running InterProScan..."
interproscan.sh \
    -i $PROTEINS \
    -o ${OUTDIR}/interpro/interpro_results.tsv \
    -f tsv,gff3 \
    -cpu $THREADS \
    -goterms \
    -pa

INTERPRO_HITS=$(cut -f1 ${OUTDIR}/interpro/interpro_results.tsv | sort -u | wc -l)
echo "InterProScan hits: $INTERPRO_HITS / $PROT_COUNT"

# Summary
echo ""
echo "=========================================="
echo "Annotation Summary"
echo "=========================================="
echo "Total proteins: $PROT_COUNT"
echo "eggNOG annotated: $EGGNOG_HITS ($(echo "scale=1; $EGGNOG_HITS * 100 / $PROT_COUNT" | bc)%)"
echo "InterProScan annotated: $INTERPRO_HITS ($(echo "scale=1; $INTERPRO_HITS * 100 / $PROT_COUNT" | bc)%)"
echo ""
echo "Run merge_annotations.py to combine results."
echo "Results in: $OUTDIR"
