#!/bin/bash
# Reference: pairtools 1.1+, bwa 0.7.17+ (or bwa-mem2 2.2+), samtools 1.19+, cooler 0.10+ | Verify API if version differs
# FASTQ -> deduplicated, filtered .pairs -> cooler, with the library-QC stats that decide whether the library worked.
# Usage: ./fastq_to_pairs.sh SAMPLE R1.fq.gz R2.fq.gz ref.fa chrom.sizes
set -euo pipefail

SAMPLE=$1
R1=$2
R2=$3
REF=$4          # bwa-indexed reference FASTA
CHROMSIZES=$5   # two-column chrom.sizes (also sets cooler bin order)

THREADS=16
MIN_MAPQ=1          # pairtools default: only MAPQ-0 treated as multi-mapping; raise to 30 for repeat-stringent loops
MAX_MISMATCH=3      # dedup wobble tolerance (bp); 0 over-splits, larger over-collapses complexity
MIN_DIST=1000       # Hi-C min-distance cut; DERIVE from the orientation-vs-distance plot (Micro-C signal is sub-1kb)
BIN=10000           # 10kb starting matrix resolution

# Align mates as INDEPENDENT single-end reads. -SP: no mate rescue/pairing (proper-pairing kills long-range/trans).
# -5: mark the 5'-most chimeric segment primary (anchors pairtools' 5' convention). -M: legacy compatibility only.
bwa mem -SP5M -t $THREADS $REF $R1 $R2 | \
    samtools view -b -@ 8 - > ${SAMPLE}.aligned.bam

# Parse to 5'-canonical pairs, sort (flips to upper-triangular), dedup (separate optical via by-tile), keep UU + UC.
pairtools parse -c $CHROMSIZES --walks-policy 5unique --min-mapq $MIN_MAPQ \
        --add-columns mapq --drop-sam ${SAMPLE}.aligned.bam | \
    pairtools sort --nproc 8 | \
    pairtools dedup --max-mismatch $MAX_MISMATCH --mark-dups \
        --output-stats ${SAMPLE}.dedup.stats \
        --output-bytile-stats ${SAMPLE}.bytile.stats | \
    pairtools select '(pair_type=="UU") or (pair_type=="UC")' \
        -o ${SAMPLE}.valid.pairs.gz

# Protocol-appropriate distance cut: keep trans plus cis beyond the orientation-equalization distance.
pairtools select "(chrom1!=chrom2) or (abs(pos2-pos1) > $MIN_DIST)" \
    -o ${SAMPLE}.filtered.pairs.gz ${SAMPLE}.valid.pairs.gz

# Library-QC readout: % long-range cis (signal), trans (noise floor), FF/FR/RF/RR orientation balance, frac_dups.
pairtools stats --bytile-dups -o ${SAMPLE}.stats.tsv ${SAMPLE}.filtered.pairs.gz

# Bin valid pairs into a cooler matrix (column order: chrom1=2 pos1=3 chrom2=4 pos2=5 in .pairs).
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
    ${CHROMSIZES}:${BIN} ${SAMPLE}.filtered.pairs.gz ${SAMPLE}.cool

echo "Done: ${SAMPLE}.cool built; QC in ${SAMPLE}.stats.tsv (read cis>=1kb, trans, orientation, frac_dups)"
