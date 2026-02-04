#!/bin/bash
# Complete alternative splicing analysis pipeline
set -e

# Configuration
SAMPLES="sample1 sample2 sample3 sample4 sample5 sample6"
CONTROL_SAMPLES="sample1 sample2 sample3"
TREATMENT_SAMPLES="sample4 sample5 sample6"
GTF="annotation.gtf"
STAR_INDEX="star_index/"
BED="annotation.bed"
THREADS=8
READ_LENGTH=150

# Create output directories
mkdir -p qc aligned rmats_output sashimi_plots

echo "Step 1: Read QC and trimming"
for sample in $SAMPLES; do
    fastp \
        -i ${sample}_R1.fastq.gz \
        -I ${sample}_R2.fastq.gz \
        -o qc/${sample}_R1.fq.gz \
        -O qc/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --thread $THREADS \
        -h qc/${sample}_fastp.html
done

echo "Step 2: STAR 2-pass alignment - Pass 1"
for sample in $SAMPLES; do
    STAR \
        --runThreadN $THREADS \
        --genomeDir $STAR_INDEX \
        --readFilesIn qc/${sample}_R1.fq.gz qc/${sample}_R2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix aligned/${sample}_p1_ \
        --outSAMtype None \
        --outSJfilterOverhangMin 8 8 8 8 \
        --alignSJDBoverhangMin 1
done

# Combine splice junctions from all samples
cat aligned/*_p1_SJ.out.tab | \
    awk '$7 >= 3' | \
    cut -f1-6 | sort -u > aligned/combined_SJ.out.tab

echo "Step 2b: STAR 2-pass alignment - Pass 2"
for sample in $SAMPLES; do
    STAR \
        --runThreadN $THREADS \
        --genomeDir $STAR_INDEX \
        --readFilesIn qc/${sample}_R1.fq.gz qc/${sample}_R2.fq.gz \
        --readFilesCommand zcat \
        --sjdbFileChrStartEnd aligned/combined_SJ.out.tab \
        --outFileNamePrefix aligned/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSJfilterOverhangMin 8 8 8 8 \
        --alignSJDBoverhangMin 1

    samtools index aligned/${sample}_Aligned.sortedByCoord.out.bam
done

echo "Step 3: Junction QC checkpoint"
for sample in $SAMPLES; do
    junction_saturation.py \
        -i aligned/${sample}_Aligned.sortedByCoord.out.bam \
        -r $BED \
        -o qc/${sample}_junc_sat
done
echo "CHECK: Verify junction saturation curves plateau before proceeding"

echo "Step 4: Create sample lists"
rm -f control_bams.txt treatment_bams.txt
for sample in $CONTROL_SAMPLES; do
    echo "aligned/${sample}_Aligned.sortedByCoord.out.bam" >> control_bams.txt
done
for sample in $TREATMENT_SAMPLES; do
    echo "aligned/${sample}_Aligned.sortedByCoord.out.bam" >> treatment_bams.txt
done

echo "Step 5: rMATS differential splicing"
rmats.py \
    --b1 control_bams.txt \
    --b2 treatment_bams.txt \
    --gtf $GTF \
    -t paired \
    --readLength $READ_LENGTH \
    --nthread $THREADS \
    --od rmats_output \
    --tmp rmats_tmp

echo "Step 6: Filter significant events"
for event in SE A5SS A3SS MXE RI; do
    awk -F'\t' 'NR==1 || ($20 < 0.05 && ($23 > 0.1 || $23 < -0.1))' \
        rmats_output/${event}.MATS.JC.txt > rmats_output/${event}_significant.txt
    echo "$event significant events: $(wc -l < rmats_output/${event}_significant.txt)"
done

echo "Pipeline complete!"
echo "Results in: rmats_output/"
echo "QC reports in: qc/"
