#!/bin/bash
# Reference: cutadapt 4.4+, umi_tools 1.1+, STAR 2.7.11+, bowtie2 2.5.3+, samtools 1.19+ | Verify API if version differs
# End-to-end Ribo-seq pipeline: preprocess -> periodicity QC -> P-site -> ORF -> TE.

set -euo pipefail

FASTQ=$1
RRNA_INDEX=$2        # bowtie2 contaminant index (rRNA + tRNA + snoRNA)
STAR_INDEX=$3
ANNOTATION=$4        # GTF
OUTDIR=${5:-riboseq_results}
ADAPTER=${6:-CTGTAGGCACCATCAAT}
HAS_UMI=${7:-no}

mkdir -p "${OUTDIR}"/{trimmed,aligned,psite}

echo "=== Step 1: UMI extract (if present) + trim ==="
READS=$FASTQ
if [ "$HAS_UMI" = "yes" ]; then
    umi_tools extract --bc-pattern=NNNNN --stdin "$FASTQ" \
        --stdout "${OUTDIR}/trimmed/umi.fastq.gz" --log "${OUTDIR}/trimmed/umi.log"
    READS="${OUTDIR}/trimmed/umi.fastq.gz"
fi
# Permissive floor + discard untrimmed (read-through is universal for footprints)
cutadapt -a "$ADAPTER" --discard-untrimmed -m 15 -M 40 \
    -o "${OUTDIR}/trimmed/trimmed.fastq.gz" "$READS" > "${OUTDIR}/trimmed/cutadapt.log"

echo "=== Step 2: rRNA removal (often >80% of reads) ==="
bowtie2 -x "$RRNA_INDEX" -U "${OUTDIR}/trimmed/trimmed.fastq.gz" \
    --un-gz "${OUTDIR}/trimmed/noncontam.fastq.gz" -S /dev/null \
    2> "${OUTDIR}/trimmed/rrna.log"

echo "=== Step 3: Alignment (EndToEnd; no soft-clipping) ==="
STAR --genomeDir "$STAR_INDEX" \
    --readFilesIn "${OUTDIR}/trimmed/noncontam.fastq.gz" --readFilesCommand zcat \
    --alignEndsType EndToEnd --seedSearchStartLmax 15 --outFilterMismatchNmax 2 \
    --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUTDIR}/aligned/" --runThreadN 8
BAM="${OUTDIR}/aligned/Aligned.sortedByCoord.out.bam"
samtools index "$BAM"

# Deduplicate only with UMIs (same position+length is mostly real co-occupancy)
if [ "$HAS_UMI" = "yes" ]; then
    umi_tools dedup --stdin "$BAM" --stdout "${OUTDIR}/aligned/dedup.bam" \
        --method directional --log "${OUTDIR}/aligned/dedup.log"
    BAM="${OUTDIR}/aligned/dedup.bam"
    samtools index "$BAM"
fi

echo "=== Step 4: P-site calibration (plastid CLI) ==="
metagene generate "${OUTDIR}/psite/cds_start" --landmark cds_start --annotation_files "$ANNOTATION"
psite "${OUTDIR}/psite/cds_start_rois.txt" "${OUTDIR}/psite/offsets" \
    --min_length 26 --max_length 34 --require_upstream --count_files "$BAM"

echo "=== Pipeline scaffold complete ==="
echo "Read-length distribution (frame-0 fraction gates downstream analysis):"
samtools view "$BAM" | awk '{print length($10)}' | sort -n | uniq -c
echo ""
echo "Transcriptome BAM for RiboCode/riboWaltz: ${OUTDIR}/aligned/Aligned.toTranscriptome.out.bam"
echo "Next: periodicity QC (riboWaltz) -> RiboCode ORFs -> riborex TE (see the ribo-seq skills)."
