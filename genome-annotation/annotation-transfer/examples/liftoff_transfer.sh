#!/bin/bash
# Reference: BioPython 1.83+, pandas 2.2+ | Verify API if version differs
# Same-species annotation transfer with Liftoff
set -euo pipefail

TARGET=$1
REFERENCE=$2
REF_ANNOTATION=$3
OUTDIR=${4:-liftoff_out}
THREADS=${5:-16}

mkdir -p $OUTDIR

echo "=== Liftoff Annotation Transfer ==="
echo "Target assembly: $TARGET"
echo "Reference genome: $REFERENCE"
echo "Reference annotation: $REF_ANNOTATION"

# Count reference features
REF_GENES=$(grep -c $'\tgene\t' $REF_ANNOTATION || echo 0)
echo "Reference genes: $REF_GENES"

# Run Liftoff with strict parameters for same-species transfer
# -sc 0.95: 95% coverage threshold (fraction of ref feature that must align)
# -s 0.90: 90% sequence identity threshold
# -exclude_partial: Exclude partially mapped features
liftoff \
    -g $REF_ANNOTATION \
    -o ${OUTDIR}/lifted_annotation.gff3 \
    -u ${OUTDIR}/unmapped_features.txt \
    -dir ${OUTDIR}/intermediates \
    -sc 0.95 \
    -s 0.90 \
    -exclude_partial \
    -p $THREADS \
    $TARGET \
    $REFERENCE

# Report transfer statistics
LIFTED_GENES=$(grep -c $'\tgene\t' ${OUTDIR}/lifted_annotation.gff3 || echo 0)
UNMAPPED=$(wc -l < ${OUTDIR}/unmapped_features.txt)

echo ""
echo "=========================================="
echo "Transfer Summary"
echo "=========================================="
echo "Reference genes: $REF_GENES"
echo "Transferred genes: $LIFTED_GENES"
echo "Unmapped features: $UNMAPPED"

if [ $REF_GENES -gt 0 ]; then
    RATE=$(echo "scale=1; $LIFTED_GENES * 100 / $REF_GENES" | bc)
    echo "Transfer rate: ${RATE}%"
fi

echo ""
echo "Transferred annotation: ${OUTDIR}/lifted_annotation.gff3"
echo "Unmapped features: ${OUTDIR}/unmapped_features.txt"
echo ""
echo "Run compare_annotations.py to validate transferred gene models."
