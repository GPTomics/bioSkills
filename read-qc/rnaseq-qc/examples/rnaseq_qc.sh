#!/bin/bash
# Reference: RSeQC 5.0+, samtools 1.19+ | Verify API if version differs
# RNA-seq specific post-alignment QC. Strandedness is inferred FIRST because it gates
# correct quantification (dUTP = reverse = featureCounts -s 2 = salmon ISR). Duplication
# is reported as a diagnostic only -- do NOT dedup non-UMI bulk RNA-seq. Needs a BED12 gene model.

BAM=$1
GENES_BED=$2

if [ -z "$BAM" ] || [ -z "$GENES_BED" ]; then
    echo "Usage: $0 <aligned.bam> <genes.bed12>"
    exit 1
fi

NAME=$(basename $BAM .bam)

echo "=== RNA-seq QC: $NAME ==="

echo -e "\n--- Strandedness (set featureCounts/salmon/Picard to match) ---"
infer_experiment.py -i $BAM -r $GENES_BED 2>/dev/null

echo -e "\n--- Read Distribution ---"
read_distribution.py -i $BAM -r $GENES_BED 2>/dev/null

echo -e "\n--- Gene Body Coverage ---"
geneBody_coverage.py -i $BAM -r $GENES_BED -o ${NAME}_coverage 2>/dev/null
echo "Plot: ${NAME}_coverage.geneBodyCoverage.curves.pdf"

echo -e "\n--- TIN Scores ---"
# tin.py writes <bam_basename>.summary.txt (mean/median/stdev TIN) and <bam_basename>.tin.xls
# (per-transcript; TIN is column 5). It does NOT print the table to stdout, so read the summary.
tin.py -i $BAM -r $GENES_BED 2>/dev/null
echo "medTIN summary:"
cat ${NAME}.summary.txt 2>/dev/null || cat *.summary.txt 2>/dev/null
