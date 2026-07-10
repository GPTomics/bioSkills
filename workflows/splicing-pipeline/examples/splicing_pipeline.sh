#!/bin/bash
# Reference: STAR 2.7.11+, fastp 0.23+, rMATS-turbo 4.3+, RSeQC 5.0+, samtools 1.19+ | Verify API if version differs
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
# --outSJfilterOverhangMin 8 8 8 8 RELAXES STAR's default (30 12 12 12) so shorter-overhang novel
# junctions survive into SJ.out.tab -- sensitivity, not stringency; the >=3 unique-read filter below
# is what removes noise. --alignSJDBoverhangMin 1 is likewise permissive for ANNOTATED junctions
# (rMATS-style; STAR default is 3) so known short-overhang junctions are not lost.
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

# Combine splice junctions from ALL samples into ONE shared DB for pass 2 (cohort-consistent PSI).
# col 7 = uniquely-mapping reads spanning the junction; require >=3 to drop noisy singletons.
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
# rMATS --b1/--b2 expect ONE line of comma-separated BAM paths, NOT one path per line.
control_bams=$(for s in $CONTROL_SAMPLES; do printf 'aligned/%s_Aligned.sortedByCoord.out.bam,' "$s"; done | sed 's/,$//')
treatment_bams=$(for s in $TREATMENT_SAMPLES; do printf 'aligned/%s_Aligned.sortedByCoord.out.bam,' "$s"; done | sed 's/,$//')
echo "$control_bams" > control_bams.txt
echo "$treatment_bams" > treatment_bams.txt

echo "Step 5: rMATS differential splicing"
rmats.py \
    --b1 control_bams.txt \
    --b2 treatment_bams.txt \
    --gtf $GTF \
    -t paired \
    --readLength $READ_LENGTH \
    --variable-read-length \
    --nthread $THREADS \
    --od rmats_output \
    --tmp rmats_tmp

echo "Step 6: Filter significant events (|deltaPSI|>0.1, FDR<0.05, >=10 junction reads/replicate)"
# Column positions DIFFER by event type: MXE has 2 extra exon-coordinate columns, so FDR and
# IncLevelDifference shift +2 for MXE. Resolve every column by HEADER NAME, never by hardcoded position.
# The read floor is not optional: PSI is a ratio, so an event backed by 2 reads can show |dPSI|=0.9
# and pass FDR on noise alone. IJC_/SJC_SAMPLE_{1,2} are comma-separated per replicate; average total
# support per replicate and require >=10, the field-convention floor this pipeline commits to.
for event in SE A5SS A3SS MXE RI; do
    f=rmats_output/${event}.MATS.JC.txt
    awk -F'\t' '
        function sumlist(s,   n,a,i,t){n=split(s,a,","); t=0; for(i=1;i<=n;i++) t+=a[i]; return t}
        NR==1{for(i=1;i<=NF;i++) h[$i]=i; print; next}
        {
            nrep = split($h["IJC_SAMPLE_1"], r1, ",") + split($h["IJC_SAMPLE_2"], r2, ",")
            support = (sumlist($h["IJC_SAMPLE_1"]) + sumlist($h["SJC_SAMPLE_1"]) + \
                       sumlist($h["IJC_SAMPLE_2"]) + sumlist($h["SJC_SAMPLE_2"])) / nrep
            fdr = $h["FDR"] + 0; d = $h["IncLevelDifference"] + 0
            if (fdr < 0.05 && (d > 0.1 || d < -0.1) && support >= 10) print
        }' "$f" > rmats_output/${event}_significant.txt
    echo "$event significant events: $(( $(wc -l < rmats_output/${event}_significant.txt) - 1 ))"
done

echo "Pipeline complete!"
echo "Results in: rmats_output/"
echo "QC reports in: qc/"
