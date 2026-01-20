#!/bin/bash
# HiFi genome assembly workflow with hifiasm

set -euo pipefail

SAMPLE="sample"
HIFI_READS="reads.hifi.fastq.gz"
HIC_R1="hic_R1.fastq.gz"
HIC_R2="hic_R2.fastq.gz"
THREADS=32
LINEAGE="mammalia_odb10"

echo "=== Step 1: Input QC ==="
seqkit stats ${HIFI_READS}

echo "=== Step 2: hifiasm assembly with Hi-C phasing ==="
hifiasm -o ${SAMPLE} -t ${THREADS} \
    --h1 ${HIC_R1} \
    --h2 ${HIC_R2} \
    ${HIFI_READS}

echo "=== Step 3: Convert GFA to FASTA ==="
for hap in hap1 hap2; do
    awk '/^S/{print ">"$2;print $3}' ${SAMPLE}.bp.${hap}.p_ctg.gfa > ${SAMPLE}.${hap}.fasta
done

awk '/^S/{print ">"$2;print $3}' ${SAMPLE}.bp.p_ctg.gfa > ${SAMPLE}.primary.fasta

echo "=== Step 4: Assembly statistics ==="
for asm in ${SAMPLE}.hap1.fasta ${SAMPLE}.hap2.fasta ${SAMPLE}.primary.fasta; do
    echo "Stats for ${asm}:"
    seqkit stats ${asm}
done

echo "=== Step 5: QUAST assessment ==="
quast.py -o quast_${SAMPLE} \
    -t ${THREADS} \
    ${SAMPLE}.hap1.fasta ${SAMPLE}.hap2.fasta ${SAMPLE}.primary.fasta

echo "=== Step 6: BUSCO completeness ==="
for asm in ${SAMPLE}.hap1.fasta ${SAMPLE}.hap2.fasta; do
    busco -i ${asm} -l ${LINEAGE} -o busco_$(basename ${asm} .fasta) -m genome -c ${THREADS}
done

echo "=== Assembly complete ==="
echo "Primary assembly: ${SAMPLE}.primary.fasta"
echo "Haplotype 1: ${SAMPLE}.hap1.fasta"
echo "Haplotype 2: ${SAMPLE}.hap2.fasta"
echo "QUAST report: quast_${SAMPLE}/report.html"
