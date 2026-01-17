#!/bin/bash
# RNA-seq specific QC

BAM=$1
GENES_BED=$2

if [ -z "$BAM" ] || [ -z "$GENES_BED" ]; then
    echo "Usage: $0 <aligned.bam> <genes.bed>"
    exit 1
fi

NAME=$(basename $BAM .bam)

echo "=== RNA-seq QC: $NAME ==="

echo -e "\n--- Strandedness ---"
infer_experiment.py -i $BAM -r $GENES_BED 2>/dev/null

echo -e "\n--- Read Distribution ---"
read_distribution.py -i $BAM -r $GENES_BED 2>/dev/null

echo -e "\n--- Gene Body Coverage ---"
geneBody_coverage.py -i $BAM -r $GENES_BED -o ${NAME}_coverage 2>/dev/null
echo "Plot: ${NAME}_coverage.geneBodyCoverage.curves.pdf"

echo -e "\n--- TIN Scores ---"
tin.py -i $BAM -r $GENES_BED 2>/dev/null | tee ${NAME}_tin.txt
mean_tin=$(awk 'NR>1 {sum+=$3; count++} END {print sum/count}' ${NAME}_tin.txt)
echo "Mean TIN: $mean_tin"
