#!/bin/bash
# Complete ATAC-seq workflow
set -e

THREADS=8
INDEX="bt2_index/genome"
GENOME="genome.fa"
GENOME_SIZE="hs"
SAMPLES="sample1 sample2 sample3"
OUTDIR="atac_results"

# Nextera adapters used in ATAC-seq
NEXTERA="CTGTCTCTTATACACATCT"

mkdir -p ${OUTDIR}/{trimmed,aligned,peaks,qc,bigwig}

echo "=== ATAC-seq Pipeline ==="
echo "Samples: ${SAMPLES}"

# === Step 1: Quality Control ===
echo "=== Step 1: Quality Control ==="
for sample in $SAMPLES; do
    echo "QC: ${sample}"
    fastp \
        -i ${sample}_R1.fastq.gz \
        -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --adapter_sequence ${NEXTERA} \
        --adapter_sequence_r2 ${NEXTERA} \
        --qualified_quality_phred 20 \
        --length_required 25 \
        --html ${OUTDIR}/qc/${sample}_fastp.html \
        -w ${THREADS}
done

# === Step 2: Alignment ===
echo "=== Step 2: Alignment ==="
for sample in $SAMPLES; do
    echo "Aligning: ${sample}"
    bowtie2 -p ${THREADS} -x ${INDEX} \
        -1 ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -2 ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --very-sensitive \
        --no-mixed --no-discordant \
        -X 2000 \
        2> ${OUTDIR}/qc/${sample}_bowtie2.log | \
    samtools view -@ ${THREADS} -bS -q 30 -f 2 - | \
    samtools sort -@ ${THREADS} -o ${OUTDIR}/aligned/${sample}.sorted.bam

    echo "Alignment stats:"
    samtools flagstat ${OUTDIR}/aligned/${sample}.sorted.bam | head -5
done

# === Step 3: BAM Processing ===
echo "=== Step 3: BAM Processing ==="
for sample in $SAMPLES; do
    echo "Processing: ${sample}"

    # Remove mitochondrial reads
    samtools view -h ${OUTDIR}/aligned/${sample}.sorted.bam | \
        grep -v chrM | \
        samtools view -b - > ${OUTDIR}/aligned/${sample}.noMT.bam

    # Mark and remove duplicates
    samtools fixmate -@ ${THREADS} -m ${OUTDIR}/aligned/${sample}.noMT.bam - | \
    samtools sort -@ ${THREADS} - | \
    samtools markdup -@ ${THREADS} -r - ${OUTDIR}/aligned/${sample}.dedup.bam

    # Shift reads for Tn5
    alignmentSieve \
        -b ${OUTDIR}/aligned/${sample}.dedup.bam \
        -o ${OUTDIR}/aligned/${sample}.shifted.bam \
        --ATACshift \
        -p ${THREADS}

    samtools index ${OUTDIR}/aligned/${sample}.shifted.bam

    # Cleanup
    rm ${OUTDIR}/aligned/${sample}.sorted.bam \
       ${OUTDIR}/aligned/${sample}.noMT.bam \
       ${OUTDIR}/aligned/${sample}.dedup.bam

    # Stats
    mt_reads=$(samtools view ${OUTDIR}/aligned/${sample}.sorted.bam 2>/dev/null | grep -c chrM || echo 0)
    total_reads=$(samtools view -c ${OUTDIR}/aligned/${sample}.shifted.bam)
    echo "Final reads (no MT, deduped, shifted): ${total_reads}"
done

# === Step 4: Peak Calling ===
echo "=== Step 4: Peak Calling ==="

# Call peaks per sample
for sample in $SAMPLES; do
    macs3 callpeak \
        -t ${OUTDIR}/aligned/${sample}.shifted.bam \
        -f BAMPE \
        -g ${GENOME_SIZE} \
        -n ${sample} \
        --outdir ${OUTDIR}/peaks \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --keep-dup all \
        -q 0.01

    echo "${sample} peaks: $(wc -l < ${OUTDIR}/peaks/${sample}_peaks.narrowPeak)"
done

# Consensus peaks
macs3 callpeak \
    -t ${OUTDIR}/aligned/*.shifted.bam \
    -f BAMPE \
    -g ${GENOME_SIZE} \
    -n consensus \
    --outdir ${OUTDIR}/peaks \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    -q 0.01

echo "Consensus peaks: $(wc -l < ${OUTDIR}/peaks/consensus_peaks.narrowPeak)"

# === Step 5: QC Metrics ===
echo "=== Step 5: QC Metrics ==="

for sample in $SAMPLES; do
    # FRiP
    total=$(samtools view -c ${OUTDIR}/aligned/${sample}.shifted.bam)
    in_peaks=$(bedtools intersect \
        -a ${OUTDIR}/aligned/${sample}.shifted.bam \
        -b ${OUTDIR}/peaks/${sample}_peaks.narrowPeak -u | samtools view -c)
    frip=$(echo "scale=4; $in_peaks / $total" | bc)
    echo "${sample} FRiP: ${frip}"

    # Fragment size distribution
    samtools view ${OUTDIR}/aligned/${sample}.shifted.bam | \
        awk '{if($9>0) print $9}' | \
        sort -n | uniq -c | \
        awk '{print $2"\t"$1}' > ${OUTDIR}/qc/${sample}_fragment_sizes.txt
done

# Generate bigWig for visualization
for sample in $SAMPLES; do
    bamCoverage \
        -b ${OUTDIR}/aligned/${sample}.shifted.bam \
        -o ${OUTDIR}/bigwig/${sample}.bw \
        --normalizeUsing RPKM \
        --ignoreDuplicates \
        -p ${THREADS}
done

echo "=== Pipeline Complete ==="
echo "Results in: ${OUTDIR}/"
echo "  - Alignments: aligned/"
echo "  - Peaks: peaks/"
echo "  - BigWig: bigwig/"
echo "  - QC: qc/"
