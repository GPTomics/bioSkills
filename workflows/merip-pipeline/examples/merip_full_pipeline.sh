#!/bin/bash
# Reference: STAR 2.7.11+, samtools 1.19+, exomePeak2 1.14+, Guitar 2.18+ | Verify API if version differs
# Complete MeRIP-seq pipeline (single-condition peak call; for differential populate bam_treated_ip/input)
set -e

SAMPLE_SHEET=$1    # CSV: sample,read1,read2,type (IP or Input),condition
STAR_INDEX=$2
GTF=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")    # absolute: the R step runs after cd into OUTPUT_DIR
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

# Peak calling. gff_dir (NOT gff) is the annotation argument. exomePeak2 auto-writes
# Mod.bed / Mod.csv / Mod.rds under save_dir/; also export a plain BED for the Guitar step.
result <- exomePeak2(
    bam_ip = ip_bams,
    bam_input = input_bams,
    gff_dir = Sys.getenv('GTF'),
    genome = 'hg38',
    paired_end = TRUE,
    save_dir = 'exomePeak2_output'
)

library(rtracklayer)
export(rowRanges(result), 'peaks/m6a_peaks.bed')   # peak GRangesList -> BED12
print(result)
RSCRIPT

cd "$OUTPUT_DIR"
GTF=$GTF Rscript peaks/call_peaks.R
cd -

echo "=== Step 4: Visualization ==="
cat > ${OUTPUT_DIR}/peaks/metagene_plot.R << 'RSCRIPT'
library(Guitar)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Guitar 2.x API: stBedFiles (site sets), txTxdb (annotation), miscOutFilePrefix (output prefix)
GuitarPlot(
    stBedFiles = list('peaks/m6a_peaks.bed'),
    txTxdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    miscOutFilePrefix = 'peaks/m6a_metagene'
)
RSCRIPT

cd "$OUTPUT_DIR"
Rscript peaks/metagene_plot.R
cd -

echo "=== Pipeline Complete ==="
echo "Results in: $OUTPUT_DIR"
echo ""
echo "Key outputs:"
echo "  - Aligned BAMs: ${OUTPUT_DIR}/aligned/"
echo "  - m6A peaks: ${OUTPUT_DIR}/peaks/m6a_peaks.bed"
echo "  - Metagene plot: ${OUTPUT_DIR}/peaks/m6a_metagene.pdf"
