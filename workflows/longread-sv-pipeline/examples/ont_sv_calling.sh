#!/bin/bash
# Reference: minimap2 2.28+, Sniffles 2.2+, bcftools 1.19+, samtools 1.19+ | Verify API if version differs
# ONT structural variant calling with Sniffles2. minimap2 >= 2.28 is required: 2.27 broke --MD.
set -e

THREADS=16
READS="nanopore_reads.fastq.gz"
REF="reference.fa"
TR_BED="human_GRCh38_TR.bed"          # reference-matched tandem-repeat BED: biggest FP lever in repeats
SAMPLE="sample1"
OUTDIR="sv_results"
PRESET="lr:hq"                        # accurate R10/Q20 ONT; use map-ont only for older noisy R9
MIN_COV=15                           # SV calling gets unreliable below ~10-15x (false negatives)
MIN_QUAL=20                          # Phred 20 = 1% error floor for a first-pass SV quality filter
MIN_SVLEN=50                         # GIAB >=50 bp SV convention (Sniffles2 default is 35)

mkdir -p ${OUTDIR}/{qc,aligned,sv}

echo "=== ONT SV Calling Pipeline ==="

# Step 1: QC
echo "=== Step 1: Quality Control ==="
NanoPlot \
    --fastq ${READS} \
    --outdir ${OUTDIR}/qc \
    --threads ${THREADS} \
    --plots dot    # 'hex' is deprecated in current NanoPlot and silently ignored unless --legacy hex
echo "QC complete. Check ${OUTDIR}/qc/NanoStats.txt (gate: read N50 >10 kb, mean quality >Q10)"

# Step 2: Alignment (-Y keeps breakpoint sequence on soft-clipped supplementary alignments)
echo "=== Step 2: Alignment ==="
minimap2 -ax ${PRESET} \
    -t ${THREADS} \
    --MD \
    -Y \
    ${REF} \
    ${READS} | \
samtools sort -@ 4 -o ${OUTDIR}/aligned/${SAMPLE}.bam
samtools index ${OUTDIR}/aligned/${SAMPLE}.bam

echo "Alignment stats:"
samtools flagstat ${OUTDIR}/aligned/${SAMPLE}.bam | head -5

avg_cov=$(samtools depth -a ${OUTDIR}/aligned/${SAMPLE}.bam | \
    awk '{sum+=$3; n++} END {printf "%.1f", sum/n}')
echo "Average coverage: ${avg_cov}x (gate: >=${MIN_COV}x for confident SV calling)"

# Step 3: SV Calling (--tandem-repeats clusters repeat-driven false positives)
echo "=== Step 3: SV Calling ==="
sniffles \
    --input ${OUTDIR}/aligned/${SAMPLE}.bam \
    --vcf ${OUTDIR}/sv/${SAMPLE}.raw.vcf.gz \
    --reference ${REF} \
    --tandem-repeats ${TR_BED} \
    --threads ${THREADS} \
    --minsvlen ${MIN_SVLEN} \
    --output-rnames

# Step 4: Filtering (ABS(SVLEN) is mandatory: DEL SVLEN is negative by convention)
# NOTE: an ABS(SVLEN) size gate silently EXCLUDES BND/translocation records (they carry no
# SVLEN), so the per-type counts below under-report translocations. Handle BNDs separately if
# translocations matter (e.g. keep SVTYPE="BND" records in a parallel pass).
echo "=== Step 4: Filtering ==="
bcftools view -i "QUAL>=${MIN_QUAL} && ABS(SVLEN)>=${MIN_SVLEN}" \
    ${OUTDIR}/sv/${SAMPLE}.raw.vcf.gz \
    -Oz -o ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz
bcftools index ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz

for svtype in DEL INS DUP INV; do
    bcftools view -i "SVTYPE=\"${svtype}\"" \
        ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz \
        -Oz -o ${OUTDIR}/sv/${SAMPLE}.${svtype}.vcf.gz
done

# Step 5: Statistics
echo "=== Step 5: Statistics ==="
bcftools stats ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz > ${OUTDIR}/sv/stats.txt

echo ""
echo "=== Summary ==="
echo "Total SVs: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz | wc -l)"
for svtype in DEL INS DUP INV; do
    echo "  ${svtype}: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.${svtype}.vcf.gz | wc -l)"
done
echo ""
echo "Results: ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz"
echo "Next: benchmark with 'truvari bench' vs GIAB HG002 Tier 1 AND CMRG; report refdist/pctsize/pctseq"
