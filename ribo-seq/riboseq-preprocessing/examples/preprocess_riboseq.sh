#!/bin/bash
# Reference: cutadapt 4.4+, umi_tools 1.1+, STAR 2.7.11+, bowtie2 2.5.3+, SortMeRNA 4.3+, samtools 1.19+ | Verify API if version differs
# Ribo-seq preprocessing: UMI-extract -> trim -> contaminant-remove -> align -> dedup -> QC.
# Canonical order from nf-core/riboseq and McGlincy & Ingolia 2017 (Methods 126:112-129).

set -euo pipefail

INPUT=$1
OUTPUT_PREFIX=${2:-riboseq}
ADAPTER=${3:-CTGTAGGCACCATCAAT}   # Example only; the real linker is protocol/kit-specific
HAS_UMI=${4:-no}                  # "yes" if the library carries UMIs
CONTAM_INDEX=${5:-contaminant_index}   # bowtie2 index of rRNA+tRNA+snoRNA+snRNA
STAR_INDEX=${6:-STAR_index}

mkdir -p logs

# Step 1: UMI extraction (BEFORE trimming, so the UMI survives in the read name).
# Pattern is library-specific; NNNNN is a placeholder for a 5-nt inline UMI.
READS=$INPUT
if [ "$HAS_UMI" = "yes" ]; then
    echo "Extracting UMIs..."
    umi_tools extract --bc-pattern=NNNNN \
        --stdin "$INPUT" --stdout "${OUTPUT_PREFIX}.umi.fastq.gz" \
        --log logs/${OUTPUT_PREFIX}_umi_extract.log
    READS="${OUTPUT_PREFIX}.umi.fastq.gz"
fi

# Step 2: Trim the 3' linker. --discard-untrimmed enriches for real footprints
# (read-through is universal, so no-adapter reads are almost never footprints).
# -m 15 is a PERMISSIVE floor: inspect the length distribution before narrowing.
echo "Trimming adapter..."
cutadapt -a "$ADAPTER" --discard-untrimmed -m 15 -M 40 -j 0 \
    -o "${OUTPUT_PREFIX}.trimmed.fastq.gz" "$READS" \
    > logs/${OUTPUT_PREFIX}_cutadapt.log 2>&1

# Step 3: Contaminant removal BEFORE alignment (rRNA is commonly 50-90% of reads).
echo "Removing contaminants..."
bowtie2 -x "$CONTAM_INDEX" \
    -U "${OUTPUT_PREFIX}.trimmed.fastq.gz" \
    --un-gz "${OUTPUT_PREFIX}.noncontam.fastq.gz" \
    -S /dev/null -p 8 \
    2> logs/${OUTPUT_PREFIX}_contaminant.log

# Step 4: Alignment with Ribo-seq-appropriate STAR flags.
# --alignEndsType EndToEnd (no soft-clipping) preserves the ends used for P-site offsets.
# --seedSearchStartLmax 15 fixes seeding for ~30 nt reads (default 50 is for long reads).
# No --alignIntronMax 1: that would forbid splicing on a genome alignment.
echo "Aligning with STAR..."
STAR --runMode alignReads \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "${OUTPUT_PREFIX}.noncontam.fastq.gz" \
    --readFilesCommand zcat \
    --alignEndsType EndToEnd \
    --seedSearchStartLmax 15 \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 10 --outSAMmultNmax 1 --outMultimapperOrder Random \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUTPUT_PREFIX}_" --runThreadN 8 \
    > logs/${OUTPUT_PREFIX}_star.log 2>&1

BAM="${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam"
samtools index "$BAM"

# Step 5: Deduplicate ONLY with UMIs. Without UMIs this step is skipped on purpose:
# many distinct ribosomes share the same 5' position and length (real co-occupancy).
if [ "$HAS_UMI" = "yes" ]; then
    echo "Deduplicating on UMIs..."
    umi_tools dedup --stdin "$BAM" --stdout "${OUTPUT_PREFIX}.dedup.bam" \
        --method directional --log logs/${OUTPUT_PREFIX}_umi_dedup.log
    BAM="${OUTPUT_PREFIX}.dedup.bam"
    samtools index "$BAM"
fi

# Step 6: QC. The read-length distribution is the key plot: expect a sharp ~28-30 nt peak
# (mammals), sometimes a ~20-22 nt shoulder (open-A-site footprints, Lareau 2014 eLife).
echo ""
echo "Read-length distribution:"
samtools view "$BAM" | awk '{print length($10)}' | sort -n | uniq -c
echo ""
echo "Alignment stats:"
samtools flagstat "$BAM"
