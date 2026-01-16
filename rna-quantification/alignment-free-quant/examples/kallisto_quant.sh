#!/bin/bash
# kallisto quantification workflow

TRANSCRIPTOME="transcripts.fa"
INDEX="kallisto_index"
THREADS=8
BOOTSTRAPS=100

# Build index (one-time)
if [ ! -f "$INDEX" ]; then
    echo "Building kallisto index..."
    kallisto index -i $INDEX $TRANSCRIPTOME
fi

# Quantify all samples
for r1 in *_R1.fastq.gz; do
    sample=$(basename $r1 _R1.fastq.gz)
    r2="${sample}_R2.fastq.gz"

    if [ ! -d "${sample}_quant" ]; then
        echo "Quantifying $sample..."
        kallisto quant -i $INDEX \
            -o ${sample}_quant \
            -t $THREADS \
            -b $BOOTSTRAPS \
            $r1 $r2
    fi
done

echo "Done! Results in *_quant/ directories"
