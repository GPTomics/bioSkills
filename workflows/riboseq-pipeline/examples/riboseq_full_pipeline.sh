#!/bin/bash
# Reference: cutadapt 4.4+, umi_tools 1.1+, STAR 2.7.11+, bowtie2 2.5.3+, plastid 0.6+, samtools 1.19+ | Verify API if version differs
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
TX_BAM="${OUTDIR}/aligned/Aligned.toTranscriptome.out.bam"
if [ "$HAS_UMI" = "yes" ]; then
    umi_tools dedup --stdin "$BAM" --stdout "${OUTDIR}/aligned/dedup.bam" \
        --method directional --log "${OUTDIR}/aligned/dedup.log"
    BAM="${OUTDIR}/aligned/dedup.bam"
    samtools index "$BAM"
    # The transcriptome BAM feeds RiboCode/riboWaltz -- dedup it too or its ORF/periodicity inputs
    # stay PCR-inflated. STAR's transcriptome BAM is unsorted, so coordinate-sort+index before umi_tools.
    samtools sort -@ 4 -o "${OUTDIR}/aligned/tx.sorted.bam" "$TX_BAM"
    samtools index "${OUTDIR}/aligned/tx.sorted.bam"
    # Position-aware dedup (NOT --per-contig/--per-gene: those treat every read on a transcript as
    # the same position, collapsing the per-codon footprints periodicity analysis depends on).
    umi_tools dedup --stdin "${OUTDIR}/aligned/tx.sorted.bam" \
        --stdout "${OUTDIR}/aligned/tx.dedup.bam" --method directional
    TX_BAM="${OUTDIR}/aligned/tx.dedup.bam"
    samtools index "$TX_BAM"
fi

echo "=== Step 4: P-site calibration (plastid CLI) ==="
metagene generate "${OUTDIR}/psite/cds_start" --landmark cds_start --annotation_files "$ANNOTATION"
# 26-34 nt: the canonical ribosome-protected footprint length window (monosome RPFs); reads
# outside it are mostly non-ribosomal contaminants and lack clean 3-nt periodicity.
psite "${OUTDIR}/psite/cds_start_rois.txt" "${OUTDIR}/psite/offsets" \
    --min_length 26 --max_length 34 --require_upstream --count_files "$BAM"

echo "=== Pipeline scaffold complete ==="
echo "Read-length distribution (frame-0 fraction gates downstream analysis):"
samtools view "$BAM" | awk '{print length($10)}' | sort -n | uniq -c
echo ""
echo "Transcriptome BAM for RiboCode/riboWaltz: ${TX_BAM}"
echo "Next: periodicity QC (riboWaltz) -> RiboCode ORFs -> riborex TE (see the ribo-seq skills)."
