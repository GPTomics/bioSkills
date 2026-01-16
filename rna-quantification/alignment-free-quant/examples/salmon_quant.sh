#!/bin/bash
# Salmon quantification workflow

TRANSCRIPTOME="Homo_sapiens.GRCh38.cdna.all.fa"
GENOME="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INDEX="salmon_index"
THREADS=8

# Build decoy-aware index (one-time)
if [ ! -d "$INDEX" ]; then
    echo "Building Salmon index with decoys..."
    grep "^>" $GENOME | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
    cat $TRANSCRIPTOME $GENOME > gentrome.fa
    salmon index -t gentrome.fa -d decoys.txt -i $INDEX -p $THREADS
    rm gentrome.fa decoys.txt
fi

# Quantify all samples
for r1 in *_R1.fastq.gz; do
    sample=$(basename $r1 _R1.fastq.gz)
    r2="${sample}_R2.fastq.gz"

    if [ ! -d "${sample}_quant" ]; then
        echo "Quantifying $sample..."
        salmon quant -i $INDEX -l A \
            -1 $r1 -2 $r2 \
            -o ${sample}_quant \
            -p $THREADS \
            --gcBias --seqBias \
            --validateMappings
    fi
done

echo "Done! Results in *_quant/ directories"
