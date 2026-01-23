#!/bin/bash
# ORF detection from Ribo-seq using RiboCode

BAM=$1
GTF=$2
GENOME=$3
OUTPUT_DIR=${4:-"ribocode_output"}

mkdir -p $OUTPUT_DIR

echo "ORF detection with RiboCode"
echo "Input: $BAM"

# One-step RiboCode analysis
# -l: Read lengths to use (28-30 nt are typical ribosome footprints)
# Adjust based on your library's read length distribution
RiboCode_onestep \
    -g $GTF \
    -r $BAM \
    -f $GENOME \
    -l 27,28,29,30 \
    -o $OUTPUT_DIR \
    2>&1 | tee ${OUTPUT_DIR}/ribocode.log

# Alternative: Step-by-step approach
# # Step 1: Prepare annotation
# prepare_transcripts \
#     -g $GTF \
#     -f $GENOME \
#     -o ${OUTPUT_DIR}/ribocode_annot

# # Step 2: Create config file
# echo -e "sample\t${BAM}\tyes" > ${OUTPUT_DIR}/config.txt

# # Step 3: Run RiboCode
# RiboCode \
#     -a ${OUTPUT_DIR}/ribocode_annot \
#     -c ${OUTPUT_DIR}/config.txt \
#     -l 27,28,29,30 \
#     -o ${OUTPUT_DIR}/ribocode

# Summarize results
RESULT_FILE=$(ls ${OUTPUT_DIR}/*_ORF_result.txt 2>/dev/null | head -1)

if [ -f "$RESULT_FILE" ]; then
    echo ""
    echo "ORF detection complete!"
    echo ""

    # Count ORFs by type
    echo "ORF counts by type:"
    # ORF types:
    # - annotated: Matches known CDS
    # - uORF: Upstream of main CDS (regulatory)
    # - dORF: Downstream of main CDS
    # - novel: In unannotated regions (potential micropeptides)
    awk -F'\t' 'NR>1 {types[$7]++} END {for (t in types) print "  "t": "types[t]}' $RESULT_FILE

    TOTAL=$(wc -l < $RESULT_FILE)
    echo "Total ORFs: $((TOTAL - 1))"

    echo ""
    echo "Output files:"
    ls -lh ${OUTPUT_DIR}/*ORF*
else
    echo "No results found. Check ribocode.log for errors."
fi
