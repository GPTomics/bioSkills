#!/bin/bash
# CNVkit workflow for exome CNV detection
set -e

GENOME="reference.fa"
TARGETS="exome_targets.bed"
REFFLAT="refFlat.txt"
OUTDIR="cnv_results"

# Sample files
NORMAL_BAMS="normal1.bam normal2.bam normal3.bam"
TUMOR_BAMS="tumor1.bam tumor2.bam"

mkdir -p ${OUTDIR}/{coverage,cnv,plots,export}

echo "=== CNVkit Pipeline ==="

# === Step 1: Prepare Targets ===
echo "=== Step 1: Prepare Targets ==="

cnvkit.py target ${TARGETS} \
    --annotate ${REFFLAT} \
    --split \
    -o ${OUTDIR}/targets.bed

cnvkit.py access ${GENOME} -o ${OUTDIR}/access.bed

cnvkit.py antitarget ${OUTDIR}/targets.bed \
    --access ${OUTDIR}/access.bed \
    -o ${OUTDIR}/antitargets.bed

echo "Targets: $(wc -l < ${OUTDIR}/targets.bed)"
echo "Antitargets: $(wc -l < ${OUTDIR}/antitargets.bed)"

# === Step 2: Calculate Coverage (Normals) ===
echo "=== Step 2: Coverage (Normals) ==="

for bam in ${NORMAL_BAMS}; do
    sample=$(basename $bam .bam)
    echo "Coverage: ${sample}"

    cnvkit.py coverage $bam ${OUTDIR}/targets.bed \
        -o ${OUTDIR}/coverage/${sample}.targetcoverage.cnn

    cnvkit.py coverage $bam ${OUTDIR}/antitargets.bed \
        -o ${OUTDIR}/coverage/${sample}.antitargetcoverage.cnn
done

# === Step 3: Create Reference ===
echo "=== Step 3: Create Reference ==="

cnvkit.py reference \
    ${OUTDIR}/coverage/normal*.targetcoverage.cnn \
    ${OUTDIR}/coverage/normal*.antitargetcoverage.cnn \
    --fasta ${GENOME} \
    -o ${OUTDIR}/reference.cnn

# === Step 4: Process Tumor Samples ===
echo "=== Step 4: Process Tumors ==="

for bam in ${TUMOR_BAMS}; do
    sample=$(basename $bam .bam)
    echo "Processing: ${sample}"

    # Coverage
    cnvkit.py coverage $bam ${OUTDIR}/targets.bed \
        -o ${OUTDIR}/coverage/${sample}.targetcoverage.cnn

    cnvkit.py coverage $bam ${OUTDIR}/antitargets.bed \
        -o ${OUTDIR}/coverage/${sample}.antitargetcoverage.cnn

    # Fix
    cnvkit.py fix \
        ${OUTDIR}/coverage/${sample}.targetcoverage.cnn \
        ${OUTDIR}/coverage/${sample}.antitargetcoverage.cnn \
        ${OUTDIR}/reference.cnn \
        -o ${OUTDIR}/cnv/${sample}.cnr

    # Segment
    cnvkit.py segment ${OUTDIR}/cnv/${sample}.cnr \
        -o ${OUTDIR}/cnv/${sample}.cns

    # Call
    cnvkit.py call ${OUTDIR}/cnv/${sample}.cns \
        -o ${OUTDIR}/cnv/${sample}.call.cns
done

# === Step 5: Visualization ===
echo "=== Step 5: Visualization ==="

for bam in ${TUMOR_BAMS}; do
    sample=$(basename $bam .bam)

    # Scatter plot
    cnvkit.py scatter ${OUTDIR}/cnv/${sample}.cnr \
        -s ${OUTDIR}/cnv/${sample}.cns \
        -o ${OUTDIR}/plots/${sample}_scatter.pdf

    # Diagram
    cnvkit.py diagram ${OUTDIR}/cnv/${sample}.cnr \
        -s ${OUTDIR}/cnv/${sample}.cns \
        -o ${OUTDIR}/plots/${sample}_diagram.pdf
done

# Heatmap
cnvkit.py heatmap ${OUTDIR}/cnv/*.cns \
    -o ${OUTDIR}/plots/heatmap.pdf

# === Step 6: Export ===
echo "=== Step 6: Export ==="

# SEG format
cnvkit.py export seg ${OUTDIR}/cnv/*.cns \
    -o ${OUTDIR}/export/all_samples.seg

# Gene metrics
for bam in ${TUMOR_BAMS}; do
    sample=$(basename $bam .bam)
    cnvkit.py genemetrics ${OUTDIR}/cnv/${sample}.cnr \
        -s ${OUTDIR}/cnv/${sample}.cns \
        --threshold 0.2 \
        -o ${OUTDIR}/export/${sample}_genes.tsv
done

echo ""
echo "=== Pipeline Complete ==="
echo "Results in: ${OUTDIR}/"
echo "  - CNV calls: cnv/"
echo "  - Plots: plots/"
echo "  - Gene summaries: export/"
