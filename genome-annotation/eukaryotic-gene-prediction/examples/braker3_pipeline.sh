#!/bin/bash
# Reference: BUSCO 5.5+, HISAT2 2.2.1+, pandas 2.2+, samtools 1.19+ | Verify API if version differs
# Eukaryotic gene prediction with BRAKER3 using RNA-seq and protein evidence
set -euo pipefail

GENOME=$1
RNASEQ_R1=$2
RNASEQ_R2=$3
PROTEINS=${4:-Viridiplantae.fa}
OUTDIR=${5:-braker3_out}
SPECIES=${6:-my_species}
THREADS=${7:-16}

echo "=== BRAKER3 Gene Prediction Pipeline ==="
echo "Genome: $GENOME"
echo "RNA-seq: $RNASEQ_R1, $RNASEQ_R2"
echo "Proteins: $PROTEINS"

# Verify softmasking
# At least 5% lowercase expected for typical eukaryotic genomes
LOWER=$(grep -v '^>' $GENOME | tr -cd 'a-z' | wc -c)
TOTAL=$(grep -v '^>' $GENOME | tr -cd 'a-zA-Z' | wc -c)
MASK_PCT=$(echo "scale=2; $LOWER * 100 / $TOTAL" | bc)
echo "Repeat masking: ${MASK_PCT}%"

if [ $(echo "$MASK_PCT < 1" | bc) -eq 1 ]; then
    echo "WARNING: Very low masking (<1%). Run RepeatMasker first."
    echo "Unmasked genomes produce many false positive gene predictions."
    exit 1
fi

# Align RNA-seq with HISAT2
echo ""
echo "Aligning RNA-seq with HISAT2..."
HISAT2_IDX=${OUTDIR}/hisat2_index/genome
mkdir -p ${OUTDIR}/hisat2_index
hisat2-build -p $THREADS $GENOME $HISAT2_IDX
hisat2 -x $HISAT2_IDX -1 $RNASEQ_R1 -2 $RNASEQ_R2 \
    --dta -p $THREADS | samtools sort -@ 4 -o ${OUTDIR}/rnaseq.bam
samtools index ${OUTDIR}/rnaseq.bam

ALIGN_RATE=$(samtools flagstat ${OUTDIR}/rnaseq.bam | grep 'mapped (' | head -1 | awk '{print $5}' | tr -d '(')
echo "RNA-seq alignment rate: $ALIGN_RATE"

# Run BRAKER3
echo ""
echo "Running BRAKER3..."
braker.pl \
    --genome=$GENOME \
    --bam=${OUTDIR}/rnaseq.bam \
    --prot_seq=$PROTEINS \
    --softmasking \
    --threads=$THREADS \
    --species=$SPECIES \
    --workingdir=$OUTDIR \
    --gff3

# Summary statistics
echo ""
echo "=========================================="
echo "Gene Prediction Summary"
echo "=========================================="
GENE_COUNT=$(grep -c $'\tgene\t' ${OUTDIR}/braker.gff3 || echo 0)
MRNA_COUNT=$(grep -c $'\tmRNA\t' ${OUTDIR}/braker.gff3 || echo 0)
EXON_COUNT=$(grep -c $'\texon\t' ${OUTDIR}/braker.gff3 || echo 0)

echo "Genes: $GENE_COUNT"
echo "mRNAs: $MRNA_COUNT"
echo "Exons: $EXON_COUNT"
if [ $GENE_COUNT -gt 0 ]; then
    echo "Avg transcripts/gene: $(echo "scale=2; $MRNA_COUNT / $GENE_COUNT" | bc)"
    echo "Avg exons/transcript: $(echo "scale=1; $EXON_COUNT / $MRNA_COUNT" | bc)"
fi

# BUSCO evaluation on predicted proteins
echo ""
echo "Running BUSCO on predicted proteins..."
busco -i ${OUTDIR}/braker.aa -m proteins -l embryophyta_odb10 -o ${OUTDIR}/busco_eval -c $THREADS --offline 2>/dev/null || echo "BUSCO skipped (lineage dataset not available)"

echo ""
echo "Results in: $OUTDIR"
