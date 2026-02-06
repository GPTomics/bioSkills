#!/bin/bash
# Repeat annotation with RepeatModeler and RepeatMasker
set -euo pipefail

ASSEMBLY=$1
OUTDIR=${2:-repeat_out}
THREADS=${3:-16}
DB_NAME=${4:-my_genome}

mkdir -p $OUTDIR

echo "=== Repeat Annotation Pipeline ==="
echo "Assembly: $ASSEMBLY"

# Build RepeatModeler database
echo ""
echo "Building RepeatModeler database..."
BuildDatabase -name ${OUTDIR}/${DB_NAME} -engine ncbi $ASSEMBLY

# Run RepeatModeler (de novo repeat library construction)
# -LTRStruct enables LTR element structural identification
echo ""
echo "Running RepeatModeler (this may take hours to days)..."
cd $OUTDIR
RepeatModeler -database ${DB_NAME} -pa $THREADS -LTRStruct
cd -

REPEAT_LIB=${OUTDIR}/${DB_NAME}-families.fa
if [ ! -f "$REPEAT_LIB" ]; then
    echo "ERROR: RepeatModeler did not produce a repeat library."
    exit 1
fi

LIB_COUNT=$(grep -c '^>' $REPEAT_LIB)
echo "De novo repeat library: $LIB_COUNT consensus sequences"

# Run RepeatMasker with de novo library
# -xsmall: Softmask (lowercase) repeats instead of hardmasking (N's)
# -gff: Produce GFF output for downstream use
echo ""
echo "Running RepeatMasker..."
RepeatMasker \
    -lib $REPEAT_LIB \
    -pa $THREADS \
    -xsmall \
    -gff \
    -dir ${OUTDIR}/repeatmasker \
    $ASSEMBLY

# Report
echo ""
echo "=========================================="
echo "Repeat Masking Summary"
echo "=========================================="
cat ${OUTDIR}/repeatmasker/$(basename $ASSEMBLY).tbl

# Verify softmasking
SOFTMASKED=${OUTDIR}/repeatmasker/$(basename $ASSEMBLY).masked
LOWER=$(grep -v '^>' $SOFTMASKED | tr -cd 'a-z' | wc -c)
TOTAL=$(grep -v '^>' $SOFTMASKED | tr -cd 'a-zA-Z' | wc -c)
MASK_PCT=$(echo "scale=2; $LOWER * 100 / $TOTAL" | bc)
echo ""
echo "Softmasking: ${MASK_PCT}% of genome in lowercase"
echo "Softmasked assembly: $SOFTMASKED"
echo ""
echo "This softmasked assembly is ready for gene prediction with BRAKER3."
