#!/bin/bash
# Reference: bowtie2 2.5+, samtools 1.19+, bedtools 2.31+, SEACR 1.4+, MACS2 2.2.9+, cutadapt 4.6+ | Verify API if version differs
# CUT&RUN / CUT&Tag end-to-end pipeline: adapter trim, bowtie2 align to combined
# host + E. coli genome, fragment bedgraph generation, SEACR + MACS2 consensus
# peak calling, spike-in scaling factor calculation.
# Henikoff lab standard bowtie2 parameters; --keep-dup all in MACS2 for CUT&Tag.

set -euo pipefail

R1=$1               # paired-end R1.fq.gz
R2=$2               # paired-end R2.fq.gz
SAMPLE=$3           # output prefix
IGG_BG=$4           # IgG control bedgraph for SEACR (or "" if not available)
GENOME_INDEX=$5     # bowtie2 index combining hg38 + E. coli K12 chromosomes
HG38_SIZES=${6:-hg38.chrom.sizes}

OUT=cnr_${SAMPLE}
mkdir -p $OUT


# === 1. Adapter trimming ===
# CUT&Tag fragments are 25-75 bp; 100-150 bp reads readthrough adapters.
# Aggressive trimming with cutadapt; minimum length 25 to discard noise.
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -e 0.1 -O 5 --minimum-length 25 \
    -o $OUT/R1.trim.fq.gz -p $OUT/R2.trim.fq.gz \
    $R1 $R2


# === 2. Alignment with Henikoff lab parameters ===
# --local: allow soft-clipping for short fragments
# --very-sensitive: more thorough alignment
# --no-mixed: only paired alignments
# --no-discordant: require concordant pair
# -I 10 -X 700: insert size range (CUT&Tag fragments up to ~700 with di-nuc)
bowtie2 \
    --local --very-sensitive --no-mixed --no-discordant \
    --phred33 -I 10 -X 700 \
    -x $GENOME_INDEX \
    -1 $OUT/R1.trim.fq.gz -2 $OUT/R2.trim.fq.gz \
    -S $OUT/aln.sam

samtools view -bS $OUT/aln.sam | samtools sort -o $OUT/aln.bam
samtools index $OUT/aln.bam
rm $OUT/aln.sam $OUT/R1.trim.fq.gz $OUT/R2.trim.fq.gz


# === 3. Spike-in: count E. coli reads ===
# E. coli alignment fraction:
#   - Target ChIP: 0.5-2% (sample with target chromatin diluting carryover)
#   - IgG control: 2-5% (no target → carryover dominates)
#   - <0.1%: spike-in lost; cross-condition normalization unreliable
# Combined index should contain E. coli contig(s) like 'NC_000913' or 'ecoli_chr1';
# adjust ECOLI_CHROMS to match what is in the index header.
ECOLI_CHROMS="NC_000913"   # E. coli K12 MG1655 RefSeq accession; adjust to local index
TOTAL_READS=$(samtools view -c -F 4 $OUT/aln.bam)
ECOLI_READS=$(samtools view -c -F 4 $OUT/aln.bam $ECOLI_CHROMS)
ECOLI_FRAC=$(echo "scale=4; $ECOLI_READS / $TOTAL_READS" | bc)
echo "Total: $TOTAL_READS; E. coli: $ECOLI_READS ($ECOLI_FRAC); expected 0.005-0.05"

# Filter to host genome only for peak calling
samtools view -b $OUT/aln.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > $OUT/aln.hg38.bam
samtools index $OUT/aln.hg38.bam


# === 4. Fragment-size diagnostic ===
# Expect 25-75 bp peak for CUT&Tag (Tn5 staggered cuts) + ~150-200 bp secondary
# (mono-nucleosomal for H3K4me3). For CUT&RUN: mono- (~150) + di-nucleosomal (~300).
samtools view -f 0x2 $OUT/aln.hg38.bam | awk '{print $9}' | awk '$1>0' \
    | sort -n | uniq -c | awk '{print $2, $1}' > $OUT/fragment_sizes.tsv


# === 5. Generate fragment bedgraph for SEACR ===
# SEACR expects bedgraph from paired-end fragments, NOT bigWig.
samtools sort -n -o $OUT/aln.hg38.qsorted.bam $OUT/aln.hg38.bam
bedtools bamtobed -bedpe -i $OUT/aln.hg38.qsorted.bam > $OUT/aln.bedpe
awk '$1==$4 && $6-$2 < 1000 {print $0}' $OUT/aln.bedpe > $OUT/aln.clean.bedpe
cut -f 1,2,6 $OUT/aln.clean.bedpe | sort -k1,1 -k2,2n -k3,3n > $OUT/fragments.bed
bedtools genomecov -bg -i $OUT/fragments.bed -g $HG38_SIZES > $OUT/aln.bedgraph


# === 6. SEACR peak calling (with IgG if available, else top 1%) ===
# norm = scale target to IgG distribution; stringent = top-half of signal blocks.
# "non" requires upstream spike-in scaling; default to "norm" otherwise.
if [ -n "$IGG_BG" ] && [ -f "$IGG_BG" ]; then
    bash SEACR_1.4.sh $OUT/aln.bedgraph $IGG_BG norm stringent $OUT/${SAMPLE}_SEACR
    echo "SEACR with IgG control: $OUT/${SAMPLE}_SEACR.stringent.bed"
else
    bash SEACR_1.4.sh $OUT/aln.bedgraph 0.01 non stringent $OUT/${SAMPLE}_SEACR_top
    echo "SEACR top-1% (no IgG): $OUT/${SAMPLE}_SEACR_top.stringent.bed"
fi


# === 7. MACS2 peak calling (--keep-dup all for CUT&Tag) ===
# -f BAMPE: paired-end fragment span (modeling on 25-75 bp fragments fails with -f BAM)
# --keep-dup all: CUT&Tag duplicates contain biology (low PCR cycles)
# -p 0.01: tighter than ChIP -q 0.05 because background is so low; or use -q 0.01
macs2 callpeak \
    -t $OUT/aln.hg38.bam \
    -f BAMPE \
    -g hs \
    -n ${SAMPLE}_MACS2 \
    --outdir $OUT/ \
    --keep-dup all \
    -p 0.01 \
    --broad-cutoff 0.1


# === 8. Consensus: intersection of SEACR and MACS2 ===
# Per btaf375 2025 benchmark, MACS2 + SEACR consensus is recommended for publication.
if [ -f $OUT/${SAMPLE}_SEACR.stringent.bed ]; then
    SEACR_BED=$OUT/${SAMPLE}_SEACR.stringent.bed
else
    SEACR_BED=$OUT/${SAMPLE}_SEACR_top.stringent.bed
fi

bedtools intersect \
    -a $OUT/${SAMPLE}_MACS2_peaks.narrowPeak \
    -b $SEACR_BED \
    -u > $OUT/${SAMPLE}_consensus.bed


# === 9. Summary ===
echo ""
echo "=== Summary: $SAMPLE ==="
echo "Total reads: $TOTAL_READS"
echo "E. coli fraction: $ECOLI_FRAC (expect 0.005-0.05; ≥0.025 for IgG)"
echo "Fragment-size mode: $(awk '{print $1, $2}' $OUT/fragment_sizes.tsv | sort -k2,2nr | head -1)"
echo "SEACR peaks: $(wc -l < $SEACR_BED)"
echo "MACS2 peaks: $(wc -l < $OUT/${SAMPLE}_MACS2_peaks.narrowPeak)"
echo "Consensus: $(wc -l < $OUT/${SAMPLE}_consensus.bed)"
