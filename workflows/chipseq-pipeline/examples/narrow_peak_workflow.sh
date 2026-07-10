#!/bin/bash
# Reference: Bowtie2 2.5.3+, MACS3 3.0+, bedtools 2.31+, fastp 0.23+, samtools 1.19+ | Verify API if version differs
# ChIP-seq workflow for narrow peaks (TFs, H3K4me3): align -> pre-dedup complexity QC ->
# dedup -> chrM + ENCODE-blacklist filter -> MACS3 (IP vs input) -> FRiP -> bigWig
set -e

THREADS=8
INDEX="bt2_index/genome"
GENOME_SIZE="hs"  # hs, mm, ce, dm -- must match the build/read length
OUTDIR="chipseq_results"
BLACKLIST="ENCODE_blacklist.bed"  # ENCODE blacklist for THIS build; subtracted before peak calling

IP_SAMPLES="IP_rep1 IP_rep2"
INPUT_SAMPLES="Input_rep1 Input_rep2"

mkdir -p ${OUTDIR}/{trimmed,aligned,peaks,qc,bigwig}

echo "=== ChIP-seq Pipeline (Narrow Peaks) ==="

# Step 1: QC
echo "=== Step 1: Quality Control ==="
for sample in $IP_SAMPLES $INPUT_SAMPLES; do
    fastp \
        -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --length_required 25 \
        --html ${OUTDIR}/qc/${sample}_fastp.html \
        -w ${THREADS}
done

# Step 2: Alignment, PRE-dedup complexity QC, dedup, chrM + blacklist filtering
echo "=== Step 2: Alignment + filtering ==="
for sample in $IP_SAMPLES $INPUT_SAMPLES; do
    echo "Aligning ${sample}..."
    A=${OUTDIR}/aligned
    bowtie2 -p ${THREADS} -x ${INDEX} \
        -1 ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -2 ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --no-mixed --no-discordant \
        --maxins 1000 \
        2> ${OUTDIR}/qc/${sample}_bowtie2.log | \
    samtools view -@ ${THREADS} -bS -q 30 -f 2 - | \
    samtools fixmate -@ ${THREADS} -m - - | \
    samtools sort -@ ${THREADS} -o ${A}/${sample}.sorted.bam

    # Library-complexity QC on the PRE-dedup BAM (NRF/PBC1 are meaningless after markdup -r)
    bedtools bamtobed -i ${A}/${sample}.sorted.bam 2>/dev/null \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
      | awk -v s=${sample} '{tot+=$1; dist++; if($1==1) one++} END{if(tot) printf "  %s NRF=%.3f PBC1=%.3f\n", s, dist/tot, one/dist}'

    samtools markdup -r -@ ${THREADS} ${A}/${sample}.sorted.bam ${A}/${sample}.dedup.bam
    samtools index ${A}/${sample}.dedup.bam

    # Drop chrM, then SUBTRACT the ENCODE blacklist (committed step, not optional cleanup)
    samtools idxstats ${A}/${sample}.dedup.bam | cut -f1 | grep -v -e '^chrM$' -e '^MT$' \
      | xargs samtools view -b ${A}/${sample}.dedup.bam > ${A}/${sample}.nochrM.bam
    if [ -f "$BLACKLIST" ]; then
        bedtools intersect -v -a ${A}/${sample}.nochrM.bam -b "$BLACKLIST" > ${A}/${sample}.bam
    else
        echo "  WARNING: $BLACKLIST not found -- peaks will contain blacklist artifacts; supply the ENCODE blacklist BED"
        mv ${A}/${sample}.nochrM.bam ${A}/${sample}.bam
    fi
    samtools index ${A}/${sample}.bam
    rm -f ${A}/${sample}.sorted.bam ${A}/${sample}.dedup.bam ${A}/${sample}.dedup.bam.bai ${A}/${sample}.nochrM.bam

    echo "Alignment stats for ${sample}:"
    samtools flagstat ${A}/${sample}.bam | head -5
done

# Step 3: Peak calling
echo "=== Step 3: Peak Calling ==="
ip_bams=""
for s in $IP_SAMPLES; do ip_bams="${ip_bams} ${OUTDIR}/aligned/${s}.bam"; done
input_bams=""
for s in $INPUT_SAMPLES; do input_bams="${input_bams} ${OUTDIR}/aligned/${s}.bam"; done

macs3 callpeak \
    -t ${ip_bams} \
    -c ${input_bams} \
    -f BAMPE \
    -g ${GENOME_SIZE} \
    -n experiment \
    --outdir ${OUTDIR}/peaks \
    -q 0.01 \
    --keep-dup all

echo "Peaks called: $(wc -l < ${OUTDIR}/peaks/experiment_peaks.narrowPeak)"

# Step 4: QC metrics
echo "=== Step 4: QC Metrics ==="
first_ip=$(echo $IP_SAMPLES | awk '{print $1}')
total=$(samtools view -c ${OUTDIR}/aligned/${first_ip}.bam)
in_peaks=$(bedtools intersect -a ${OUTDIR}/aligned/${first_ip}.bam \
    -b ${OUTDIR}/peaks/experiment_peaks.narrowPeak -u -wa | samtools view -c)
frip=$(echo "scale=4; $in_peaks / $total" | bc)
echo "FRiP for ${first_ip}: ${frip}"

# Step 5: Generate bigWig
echo "=== Step 5: Generate bigWig tracks ==="
for sample in $IP_SAMPLES; do
    bamCoverage \
        -b ${OUTDIR}/aligned/${sample}.bam \
        -o ${OUTDIR}/bigwig/${sample}.bw \
        --normalizeUsing RPKM \
        --ignoreDuplicates \
        -p ${THREADS}
done

echo "=== Pipeline Complete ==="
echo "Peaks: ${OUTDIR}/peaks/experiment_peaks.narrowPeak"
echo "BigWig: ${OUTDIR}/bigwig/"
