#!/usr/bin/env bash
# Reference: STAR 2.7.11+, HISAT2 2.2.1+, samtools 1.19+, fastp 0.23+ | Verify with `STAR --version`, `hisat2 --version`, `samtools --version`, `fastp --version` if installed releases differ.
# MeRIP-seq end-to-end alignment: trim -> STAR splice-aware -> sort -> index -> per-BAM QC.
# Standard non-UMI MeRIP: do NOT deduplicate (collapses real coverage at high-expression transcripts).

set -euo pipefail

# Sample inventory: pairs of (IP, matched_input) per biological replicate.
SAMPLES=(IP_rep1 IP_rep2 IP_rep3 Input_rep1 Input_rep2 Input_rep3)

REFERENCE_GENOME='refs/genome.fa'
REFERENCE_GTF='refs/annotation.gtf'
STAR_INDEX='refs/star_index'
THREADS=12
MIN_READ_LENGTH=25

mkdir -p trimmed aligned qc

# Build STAR index once if not present.
if [ ! -d "${STAR_INDEX}/SA" ]; then
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${REFERENCE_GENOME}" \
        --sjdbGTFfile "${REFERENCE_GTF}" \
        --sjdbOverhang 100 \
        --runThreadN "${THREADS}"
fi

for sample in "${SAMPLES[@]}"; do
    # Adapter trimming with fastp. Standard non-UMI MeRIP: no --umi flag.
    fastp \
        --in1 "raw/${sample}_R1.fastq.gz" \
        --in2 "raw/${sample}_R2.fastq.gz" \
        --out1 "trimmed/${sample}_R1.fq.gz" \
        --out2 "trimmed/${sample}_R2.fq.gz" \
        --html "qc/${sample}_fastp.html" \
        --json "qc/${sample}_fastp.json" \
        --length_required "${MIN_READ_LENGTH}" \
        --detect_adapter_for_pe \
        --thread "${THREADS}"

    # STAR splice-aware alignment. Multi-mapper retention up to 20 for MeRIP at multi-isoform loci.
    STAR \
        --runMode alignReads \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "trimmed/${sample}_R1.fq.gz" "trimmed/${sample}_R2.fq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 20 \
        --outSAMattributes NH HI AS nM NM MD \
        --outFileNamePrefix "aligned/${sample}_" \
        --runThreadN "${THREADS}"

    bam="aligned/${sample}_Aligned.sortedByCoord.out.bam"

    samtools index -@ 4 "${bam}"
    samtools flagstat "${bam}" > "qc/${sample}.flagstat"
    samtools idxstats "${bam}" > "qc/${sample}.idxstats"
done

echo 'Alignment complete. Skipped dedup (standard non-UMI MeRIP).'
echo 'Next: deepTools multiBamSummary + plotFingerprint + PreSeq lc_extrap; then exomePeak2 / MeTPeak / MACS3 peak calling.'
