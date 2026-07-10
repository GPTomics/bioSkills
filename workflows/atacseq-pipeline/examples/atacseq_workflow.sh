#!/bin/bash
# Reference: Bowtie2 2.5.3+, MACS3 3.0+, bedtools 2.31+, deepTools 3.5+, fastp 0.23+, samtools 1.19+ | Verify API if version differs
# Complete ATAC-seq workflow
set -e

THREADS=8
INDEX="bt2_index/genome"
GENOME_SIZE="hs"
SAMPLES="sample1 sample2 sample3"
OUTDIR="atac_results"
BLACKLIST="ENCODE_blacklist.bed"   # ENCODE blacklist for THIS build; ATAC has no input, so this
                                   # + the shift-extend model ARE the background -- subtract before calling

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

    # Remove mitochondrial reads (BEFORE dedup/peaks; mito can be 20-50% of an ATAC library).
    # Tab-anchored so 'chrM'/'MT' matches the RNAME field, not a substring of a read name.
    samtools view -h ${OUTDIR}/aligned/${sample}.sorted.bam | \
        grep -v -e $'\tchrM\t' -e $'\tMT\t' | \
        samtools view -b - > ${OUTDIR}/aligned/${sample}.noMT.bam

    # Library-complexity QC on the PRE-dedup BAM. `markdup -r` below physically removes duplicates, after
    # which NRF/PBC1 are identically 1.0 and meaningless. Computed on noMT.bam, not sorted.bam: mito reads
    # are over-amplified, so including them measures chrM chemistry rather than nuclear-library complexity.
    # NRF = distinct/total, PBC1 = read-once/distinct (ENCODE); >0.8 acceptable, >0.9 preferred.
    bedtools bamtobed -i ${OUTDIR}/aligned/${sample}.noMT.bam 2>/dev/null \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
      | awk -v s=${sample} '{tot+=$1; dist++; if($1==1) one++} END{if(tot) printf "  %s NRF=%.3f PBC1=%.3f\n", s, dist/tot, one/dist}'

    # Dedup: collate FIRST -- fixmate -m needs name-grouped input, but noMT.bam is coordinate-sorted,
    # so running fixmate directly would mis-pair mates and markdup would mis-flag duplicates.
    samtools collate -@ ${THREADS} -O -u ${OUTDIR}/aligned/${sample}.noMT.bam | \
    samtools fixmate -@ ${THREADS} -m -u - - | \
    samtools sort -@ ${THREADS} -u - | \
    samtools markdup -@ ${THREADS} -r - ${OUTDIR}/aligned/${sample}.dedup.bam
    samtools index ${OUTDIR}/aligned/${sample}.dedup.bam   # alignmentSieve needs an indexed input BAM

    # Tn5 +4/-5 shift ONCE
    alignmentSieve \
        -b ${OUTDIR}/aligned/${sample}.dedup.bam \
        -o ${OUTDIR}/aligned/${sample}.shift.bam \
        --ATACshift \
        -p ${THREADS}
    samtools index ${OUTDIR}/aligned/${sample}.shift.bam

    # Subtract the ENCODE blacklist (committed step; ATAC has no input control)
    if [ -f "$BLACKLIST" ]; then
        bedtools intersect -v -a ${OUTDIR}/aligned/${sample}.shift.bam -b "$BLACKLIST" \
            > ${OUTDIR}/aligned/${sample}.shifted.bam
    else
        echo "  WARNING: $BLACKLIST not found -- peaks will contain blacklist artifacts; supply the ENCODE blacklist BED"
        mv ${OUTDIR}/aligned/${sample}.shift.bam ${OUTDIR}/aligned/${sample}.shifted.bam
    fi
    samtools index ${OUTDIR}/aligned/${sample}.shifted.bam

    # Cleanup intermediates
    rm -f ${OUTDIR}/aligned/${sample}.sorted.bam \
       ${OUTDIR}/aligned/${sample}.noMT.bam \
       ${OUTDIR}/aligned/${sample}.dedup.bam ${OUTDIR}/aligned/${sample}.dedup.bam.bai \
       ${OUTDIR}/aligned/${sample}.shift.bam ${OUTDIR}/aligned/${sample}.shift.bam.bai

    total_reads=$(samtools view -c ${OUTDIR}/aligned/${sample}.shifted.bam)
    echo "Final reads (no MT, deduped, shifted): ${total_reads}"
done

# === Step 4: Peak Calling ===
echo "=== Step 4: Peak Calling ==="

# Call peaks per sample. Cut-site mode: -f BAM (single-end-ified shifted reads) so --shift/--extsize
# apply. Do NOT use -f BAMPE with --shift/--extsize -- BAMPE ignores them silently.
for sample in $SAMPLES; do
    macs3 callpeak \
        -t ${OUTDIR}/aligned/${sample}.shifted.bam \
        -f BAM \
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

# Pooled peaks (a fixed-width Corces consensus for differential is built separately; see
# atac-seq/consensus-peakset). Same cut-site mode.
macs3 callpeak \
    -t ${OUTDIR}/aligned/*.shifted.bam \
    -f BAM \
    -g ${GENOME_SIZE} \
    -n consensus \
    --outdir ${OUTDIR}/peaks \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    --keep-dup all \
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
