#!/bin/bash
# Reference: stjude/ROSE (Python 3), samtools 1.19+ | Verify API if version differs
# Super-enhancer calling with the stjude/ROSE Python-3 fork (genomeDict includes HG38;
# run from the ROSE repo dir so ROSE_main.py finds its annotation/ files)

set -euo pipefail

PEAKS=${1:-H3K27ac_peaks.bed}
BAM=${2:-H3K27ac.bam}
GENOME=${3:-hg38}
OUTPUT_DIR=${4:-rose_output}

echo "=== ROSE Super-Enhancer Calling ==="
echo "Peaks: $PEAKS"
echo "BAM: $BAM"
echo "Genome: $GENOME"

# Convert BED to GFF if needed
if [[ "$PEAKS" == *.bed ]]; then
    echo "Converting BED to GFF..."
    awk 'BEGIN{OFS="\t"} {print $1,"MACS","peak",($2+1),$3,".",$6,".","peak_"NR}' $PEAKS > peaks.gff
    PEAKS_GFF=peaks.gff
else
    PEAKS_GFF=$PEAKS
fi

# Run stjude/ROSE (Python-3 fork; ROSE_main.py supports -g HG38)
python ROSE_main.py \
    -g $GENOME \
    -i $PEAKS_GFF \
    -r $BAM \
    -o $OUTPUT_DIR \
    -s 12500 \
    -t 2500

echo ""
echo "=== Results ==="
echo "All enhancers: ${OUTPUT_DIR}/*_AllEnhancers.table.txt"
echo "Super-enhancers: ${OUTPUT_DIR}/*_SuperEnhancers.table.txt"

# Count super-enhancers
n_se=$(tail -n +7 ${OUTPUT_DIR}/*_SuperEnhancers.table.txt 2>/dev/null | wc -l || echo "0")
echo "Number of super-enhancers: $n_se"
