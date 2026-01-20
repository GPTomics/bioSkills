#!/bin/bash
# Metagenome assembly and binning workflow

set -euo pipefail

SAMPLE="metagenome"
READS="ont_reads.fastq.gz"
THREADS=32

echo "=== Step 1: Input QC ==="
seqkit stats ${READS}

echo "=== Step 2: metaFlye assembly ==="
flye --nano-raw ${READS} \
    --meta \
    --out-dir ${SAMPLE}_flye \
    --threads ${THREADS}

ASSEMBLY="${SAMPLE}_flye/assembly.fasta"
echo "Assembly stats:"
seqkit stats ${ASSEMBLY}

echo "=== Step 3: Map reads for binning ==="
minimap2 -ax map-ont -t ${THREADS} ${ASSEMBLY} ${READS} | \
    samtools sort -@ ${THREADS} -o ${SAMPLE}.bam -
samtools index ${SAMPLE}.bam

echo "=== Step 4: Generate depth profile ==="
jgi_summarize_bam_contig_depths --outputDepth ${SAMPLE}_depth.txt ${SAMPLE}.bam

echo "=== Step 5: Bin with MetaBAT2 ==="
mkdir -p ${SAMPLE}_bins
metabat2 -i ${ASSEMBLY} -a ${SAMPLE}_depth.txt \
    -o ${SAMPLE}_bins/bin \
    -t ${THREADS} \
    --minContig 1500

echo "Number of bins: $(ls ${SAMPLE}_bins/*.fa 2>/dev/null | wc -l)"

echo "=== Step 6: Assess bin quality ==="
checkm2 predict \
    --input ${SAMPLE}_bins/ \
    --output-directory ${SAMPLE}_checkm2 \
    -x fa \
    --threads ${THREADS}

echo "=== Step 7: Identify circular genomes ==="
grep "Y" ${SAMPLE}_flye/assembly_info.txt | cut -f1 > circular_contigs.txt
echo "Circular contigs found: $(wc -l < circular_contigs.txt)"

if [ -s circular_contigs.txt ]; then
    seqkit grep -f circular_contigs.txt ${ASSEMBLY} > ${SAMPLE}_circular.fasta
fi

echo "=== Step 8: Taxonomic classification (optional) ==="
# Uncomment if GTDB-Tk database is available
# gtdbtk classify_wf --genome_dir ${SAMPLE}_bins/ \
#     --out_dir ${SAMPLE}_gtdbtk \
#     --cpus ${THREADS} \
#     -x fa

echo "=== Workflow complete ==="
echo "Assembly: ${ASSEMBLY}"
echo "Bins: ${SAMPLE}_bins/"
echo "CheckM2 report: ${SAMPLE}_checkm2/quality_report.tsv"
