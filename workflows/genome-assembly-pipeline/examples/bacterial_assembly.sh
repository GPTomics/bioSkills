#!/bin/bash
# Bacterial genome assembly from ONT + Illumina reads
set -e

THREADS=16
GENOME_SIZE="5m"
LONG_READS="nanopore.fastq.gz"
SHORT_R1="illumina_R1.fastq.gz"
SHORT_R2="illumina_R2.fastq.gz"
OUTDIR="bacterial_assembly"

mkdir -p ${OUTDIR}/{qc,flye,medaka,pilon,final,quast,busco}

echo "=== Bacterial Genome Assembly Pipeline ==="

# Step 1: QC
echo "=== Step 1: Quality Control ==="
NanoPlot --fastq ${LONG_READS} --outdir ${OUTDIR}/qc/nanoplot -t ${THREADS}
fastp -i ${SHORT_R1} -I ${SHORT_R2} \
    -o ${OUTDIR}/qc/trimmed_R1.fq.gz \
    -O ${OUTDIR}/qc/trimmed_R2.fq.gz \
    --detect_adapter_for_pe \
    --html ${OUTDIR}/qc/fastp.html -w ${THREADS}

# Step 2: Long-read assembly with Flye
echo "=== Step 2: Flye Assembly ==="
flye --nano-hq ${LONG_READS} \
    --out-dir ${OUTDIR}/flye \
    --threads ${THREADS} \
    --genome-size ${GENOME_SIZE}

echo "Flye assembly stats:"
head -20 ${OUTDIR}/flye/assembly_info.txt

# Step 3: Polish with medaka (long reads)
echo "=== Step 3: medaka Polishing ==="
medaka_consensus \
    -i ${LONG_READS} \
    -d ${OUTDIR}/flye/assembly.fasta \
    -o ${OUTDIR}/medaka \
    -t ${THREADS}

# Step 4: Polish with Pilon (short reads)
echo "=== Step 4: Pilon Polishing ==="
bwa index ${OUTDIR}/medaka/consensus.fasta
bwa mem -t ${THREADS} ${OUTDIR}/medaka/consensus.fasta \
    ${OUTDIR}/qc/trimmed_R1.fq.gz ${OUTDIR}/qc/trimmed_R2.fq.gz | \
    samtools sort -@ 4 -o ${OUTDIR}/pilon/aligned.bam
samtools index ${OUTDIR}/pilon/aligned.bam

pilon --genome ${OUTDIR}/medaka/consensus.fasta \
    --frags ${OUTDIR}/pilon/aligned.bam \
    --output ${OUTDIR}/pilon/polished \
    --threads ${THREADS} \
    --changes

# Copy final assembly
cp ${OUTDIR}/pilon/polished.fasta ${OUTDIR}/final/assembly.fasta

# Step 5: QUAST
echo "=== Step 5: QUAST ==="
quast.py ${OUTDIR}/final/assembly.fasta \
    -o ${OUTDIR}/quast \
    -t ${THREADS}

# Step 6: BUSCO
echo "=== Step 6: BUSCO ==="
busco -i ${OUTDIR}/final/assembly.fasta \
    -l bacteria_odb10 \
    -o busco \
    -m genome \
    -c ${THREADS} \
    --out_path ${OUTDIR}

# Summary
echo ""
echo "=== Assembly Complete ==="
echo "Final assembly: ${OUTDIR}/final/assembly.fasta"
echo ""
echo "QUAST Summary:"
cat ${OUTDIR}/quast/report.txt | grep -E "contigs|Total length|N50|L50"
echo ""
echo "BUSCO Summary:"
cat ${OUTDIR}/busco/short_summary*.txt | grep -E "Complete|Fragmented|Missing"
