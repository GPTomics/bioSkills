#!/bin/bash
# Complete MeRIP-seq pipeline

SAMPLE_SHEET=$1    # CSV: sample,read1,read2,type (IP or Input),condition
STAR_INDEX=$2
GTF=$3
OUTPUT_DIR=${4:-"merip_results"}
THREADS=${5:-8}

mkdir -p ${OUTPUT_DIR}/{aligned,peaks}

echo "=== Step 1: Alignment ==="
# Align all IP and Input samples
while IFS=',' read -r sample r1 r2 type condition; do
    echo "Aligning $sample ($type, $condition)..."

    STAR --genomeDir $STAR_INDEX \
        --readFilesIn $r1 $r2 \
        --readFilesCommand zcat \
        --runThreadN $THREADS \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUTPUT_DIR}/aligned/${sample}_

    samtools index ${OUTPUT_DIR}/aligned/${sample}_Aligned.sortedByCoord.out.bam

done < $SAMPLE_SHEET

echo "=== Step 2: QC - IP Enrichment ==="
# Check IP vs Input coverage patterns
# IP should show peaks, Input should be uniform
for bam in ${OUTPUT_DIR}/aligned/*IP*.bam; do
    samtools flagstat $bam > ${bam%.bam}_flagstat.txt
done

echo "=== Step 3: Peak Calling (exomePeak2) ==="
# Run in R - creates R script for peak calling
cat > ${OUTPUT_DIR}/peaks/call_peaks.R << 'RSCRIPT'
library(exomePeak2)

# Get BAM files
ip_bams <- list.files('aligned', pattern = 'IP.*\\.bam$', full.names = TRUE)
input_bams <- list.files('aligned', pattern = 'Input.*\\.bam$', full.names = TRUE)

# Peak calling
result <- exomePeak2(
    bam_ip = ip_bams,
    bam_input = input_bams,
    gff = Sys.getenv('GTF'),
    genome = 'hg38',
    paired_end = TRUE,
    p_cutoff = 0.05,
    log2FC_cutoff = 1
)

# Export
exportResults(result, format = 'BED', file = 'peaks/m6a_peaks.bed')
exportResults(result, format = 'CSV', file = 'peaks/m6a_peaks.csv')

print(result)
RSCRIPT

cd $OUTPUT_DIR
GTF=$GTF Rscript peaks/call_peaks.R
cd -

echo "=== Step 4: Visualization ==="
cat > ${OUTPUT_DIR}/peaks/metagene_plot.R << 'RSCRIPT'
library(Guitar)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

peaks <- import('peaks/m6a_peaks.bed')

GuitarPlot(
    peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    saveToPDFprefix = 'peaks/m6a_metagene'
)
RSCRIPT

cd $OUTPUT_DIR
Rscript peaks/metagene_plot.R
cd -

echo "=== Pipeline Complete ==="
echo "Results in: $OUTPUT_DIR"
echo ""
echo "Key outputs:"
echo "  - Aligned BAMs: ${OUTPUT_DIR}/aligned/"
echo "  - m6A peaks: ${OUTPUT_DIR}/peaks/m6a_peaks.bed"
echo "  - Metagene plot: ${OUTPUT_DIR}/peaks/m6a_metagene.pdf"
