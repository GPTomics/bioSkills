#!/bin/bash
# Complete metagenomics workflow: Kraken2 + Bracken + HUMAnN
set -e

THREADS=8
KRAKEN_DB="/path/to/kraken2_standard"
HOST_INDEX="/path/to/human_bt2"
SAMPLES="sample1 sample2 sample3"
OUTDIR="metagenomics_results"

mkdir -p ${OUTDIR}/{trimmed,host_removed,kraken,bracken,humann,qc,viz}

echo "=== Metagenomics Pipeline ==="
echo "Samples: ${SAMPLES}"

# === Step 1: Quality Control ===
echo "=== Step 1: QC ==="
for sample in $SAMPLES; do
    echo "QC: ${sample}"
    fastp \
        -i ${sample}_R1.fastq.gz \
        -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --complexity_threshold 30 \
        --html ${OUTDIR}/qc/${sample}_fastp.html \
        -w ${THREADS}
done

# === Step 2: Host Removal ===
echo "=== Step 2: Host Removal ==="
for sample in $SAMPLES; do
    echo "Removing host: ${sample}"
    bowtie2 -p ${THREADS} -x ${HOST_INDEX} \
        -1 ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -2 ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --un-conc-gz ${OUTDIR}/host_removed/${sample}_R%.fq.gz \
        --very-sensitive \
        > /dev/null 2> ${OUTDIR}/qc/${sample}_host_removal.log

    # Stats
    total=$(zcat ${OUTDIR}/trimmed/${sample}_R1.fq.gz | wc -l)
    remaining=$(zcat ${OUTDIR}/host_removed/${sample}_R1.fq.gz | wc -l)
    host_pct=$(echo "scale=2; 100 - ($remaining / $total * 100)" | bc)
    echo "${sample}: ${host_pct}% host reads removed"
done

# === Step 3: Kraken2 Classification ===
echo "=== Step 3: Kraken2 Classification ==="
for sample in $SAMPLES; do
    echo "Classifying: ${sample}"
    kraken2 \
        --db ${KRAKEN_DB} \
        --threads ${THREADS} \
        --paired \
        --report ${OUTDIR}/kraken/${sample}.report \
        --report-minimizer-data \
        --output ${OUTDIR}/kraken/${sample}.output \
        ${OUTDIR}/host_removed/${sample}_R1.fq.gz \
        ${OUTDIR}/host_removed/${sample}_R2.fq.gz

    # Classification rate
    classified=$(grep -m1 "root" ${OUTDIR}/kraken/${sample}.report | awk '{print $1}')
    echo "${sample}: ${classified}% classified"
done

# === Step 4: Bracken Abundance ===
echo "=== Step 4: Bracken Abundance ==="
for sample in $SAMPLES; do
    echo "Estimating abundance: ${sample}"

    # Species level
    bracken -d ${KRAKEN_DB} \
        -i ${OUTDIR}/kraken/${sample}.report \
        -o ${OUTDIR}/bracken/${sample}.species.txt \
        -w ${OUTDIR}/bracken/${sample}.species.report \
        -r 150 -l S -t 10

    # Genus level
    bracken -d ${KRAKEN_DB} \
        -i ${OUTDIR}/kraken/${sample}.report \
        -o ${OUTDIR}/bracken/${sample}.genus.txt \
        -w ${OUTDIR}/bracken/${sample}.genus.report \
        -r 150 -l G -t 10
done

# Combine tables
echo "Combining abundance tables..."
combine_bracken_outputs.py \
    --files ${OUTDIR}/bracken/*.species.txt \
    -o ${OUTDIR}/bracken/combined_species.txt

# === Step 5: HUMAnN Functional Profiling (Optional) ===
echo "=== Step 5: HUMAnN Functional Profiling ==="
for sample in $SAMPLES; do
    echo "Functional profiling: ${sample}"

    # Concatenate paired reads
    cat ${OUTDIR}/host_removed/${sample}_R1.fq.gz \
        ${OUTDIR}/host_removed/${sample}_R2.fq.gz > \
        ${OUTDIR}/host_removed/${sample}_concat.fq.gz

    humann \
        --input ${OUTDIR}/host_removed/${sample}_concat.fq.gz \
        --output ${OUTDIR}/humann/${sample} \
        --threads ${THREADS}

    rm ${OUTDIR}/host_removed/${sample}_concat.fq.gz
done

# Join HUMAnN tables
humann_join_tables -i ${OUTDIR}/humann -o ${OUTDIR}/humann/merged_pathabundance.tsv \
    --file_name pathabundance
humann_join_tables -i ${OUTDIR}/humann -o ${OUTDIR}/humann/merged_genefamilies.tsv \
    --file_name genefamilies

# Normalize
humann_renorm_table -i ${OUTDIR}/humann/merged_pathabundance.tsv \
    -o ${OUTDIR}/humann/merged_pathabundance_cpm.tsv -u cpm

echo "=== Pipeline Complete ==="
echo "Results:"
echo "  Kraken2 reports: ${OUTDIR}/kraken/"
echo "  Bracken abundances: ${OUTDIR}/bracken/"
echo "  HUMAnN pathways: ${OUTDIR}/humann/"
