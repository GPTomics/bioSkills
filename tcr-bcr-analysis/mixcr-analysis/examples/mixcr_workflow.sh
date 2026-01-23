#!/bin/bash
# MiXCR workflow for TCR/BCR repertoire analysis

R1=$1
R2=$2
OUTPUT_PREFIX=${3:-"sample"}
SPECIES=${4:-"human"}
CHAIN=${5:-"TRB"}  # TRA, TRB, TRG, TRD, IGH, IGK, IGL

echo "MiXCR analysis: $OUTPUT_PREFIX"
echo "Species: $SPECIES, Chain: $CHAIN"

# Method 1: One-command workflow (recommended for most cases)
# Amplicon protocol for targeted TCR/BCR sequencing
mixcr analyze amplicon \
    -s $SPECIES \
    --starting-material rna \
    --5-end v-primers \
    --3-end c-primers \
    --adapters adapters-cdr3 \
    $R1 $R2 \
    ${OUTPUT_PREFIX}

# Method 2: Step-by-step (for custom analysis)
# Uncomment to use instead of analyze command

# # Step 1: Align reads to V(D)J reference
# mixcr align \
#     -s $SPECIES \
#     -p rna-seq \
#     -OallowPartialAlignments=true \
#     $R1 $R2 \
#     ${OUTPUT_PREFIX}.vdjca

# # Step 2: Rescue partial alignments (improves recovery)
# mixcr assemblePartial \
#     ${OUTPUT_PREFIX}.vdjca \
#     ${OUTPUT_PREFIX}_rescued.vdjca

# # Step 3: Extend alignments
# mixcr extend \
#     ${OUTPUT_PREFIX}_rescued.vdjca \
#     ${OUTPUT_PREFIX}_extended.vdjca

# # Step 4: Assemble clonotypes
# mixcr assemble \
#     ${OUTPUT_PREFIX}_extended.vdjca \
#     ${OUTPUT_PREFIX}.clns

# Export clonotypes with all relevant information
mixcr exportClones \
    -c $CHAIN \
    -cloneId \
    -count \
    -fraction \
    -nFeature CDR3 \
    -aaFeature CDR3 \
    -vGene \
    -dGene \
    -jGene \
    ${OUTPUT_PREFIX}.clns \
    ${OUTPUT_PREFIX}_clones.txt

# Export alignment report
mixcr exportReports ${OUTPUT_PREFIX}.clns > ${OUTPUT_PREFIX}_report.txt

# Summary statistics
echo ""
echo "Analysis complete!"
echo "Alignment statistics:"
grep -E "Successfully aligned|CDR3" ${OUTPUT_PREFIX}_report.txt 2>/dev/null || echo "See ${OUTPUT_PREFIX}_report.txt"

CLONE_COUNT=$(wc -l < ${OUTPUT_PREFIX}_clones.txt)
echo "Clonotypes identified: $((CLONE_COUNT - 1))"  # Subtract header

echo ""
echo "Output files:"
ls -lh ${OUTPUT_PREFIX}*
