#!/bin/bash
# Reference: MACS2 2.2.9+, MACS3 3.0+, samtools 1.19+, idr 2.0.4+, bedtools 2.31+ | Verify API if version differs
# ENCODE 4 ATAC-seq peak-calling pipeline: per-replicate MACS2 + pseudoreplicate IDR self-consistency.
# `-f BAM --shift -75 --extsize 150` is the ENCODE pattern; `-f BAMPE` would silently ignore --shift/--extsize.

set -euo pipefail

REP1=${1:-rep1.filt.dedup.bam}
REP2=${2:-rep2.filt.dedup.bam}
GENOME=${3:-hs}                                  # Or e.g. 2.701e9 for hg38 100bp reads (deepTools)
BLACKLIST=${4:-hg38-blacklist.v2.bed}            # ENCODE blacklist (Amemiya 2019)
OUTDIR=${5:-peaks_out}

mkdir -p $OUTDIR/{rep1,rep2,pooled,psr1_1,psr1_2,psr2_1,psr2_2,idr}

call_peaks() {
    local bam=$1 name=$2 outdir=$3
    macs2 callpeak \
        -t $bam \
        -f BAM -g $GENOME \
        -n $name --outdir $outdir \
        --nomodel --shift -75 --extsize 150 \
        --keep-dup all \
        -B --SPMR \
        -p 0.01
}

# 1. Per-replicate
call_peaks $REP1 rep1 $OUTDIR/rep1
call_peaks $REP2 rep2 $OUTDIR/rep2

# 2. Pooled
macs2 callpeak \
    -t $REP1 $REP2 \
    -f BAM -g $GENOME -n pooled --outdir $OUTDIR/pooled \
    --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01

# 3. Pseudoreplicates (two independent 50% subsamples per rep, different seeds; approximate, not a disjoint partition)
samtools view -b -h -s 1.5 $REP1 > $OUTDIR/psr1_1/rep1.psr1.bam
samtools view -b -h -s 2.5 $REP1 > $OUTDIR/psr1_2/rep1.psr2.bam
samtools view -b -h -s 1.5 $REP2 > $OUTDIR/psr2_1/rep2.psr1.bam
samtools view -b -h -s 2.5 $REP2 > $OUTDIR/psr2_2/rep2.psr2.bam
for f in $OUTDIR/psr*/rep*.bam; do samtools index $f; done

call_peaks $OUTDIR/psr1_1/rep1.psr1.bam rep1_psr1 $OUTDIR/psr1_1
call_peaks $OUTDIR/psr1_2/rep1.psr2.bam rep1_psr2 $OUTDIR/psr1_2

# 4. IDR on true reps (threshold 0.05) and on pseudoreps (threshold 0.10)
sort_pk() { sort -k8,8nr $1 > ${1%.narrowPeak}.sorted.narrowPeak; }
sort_pk $OUTDIR/rep1/rep1_peaks.narrowPeak
sort_pk $OUTDIR/rep2/rep2_peaks.narrowPeak
sort_pk $OUTDIR/psr1_1/rep1_psr1_peaks.narrowPeak
sort_pk $OUTDIR/psr1_2/rep1_psr2_peaks.narrowPeak

idr --samples \
        $OUTDIR/rep1/rep1_peaks.sorted.narrowPeak \
        $OUTDIR/rep2/rep2_peaks.sorted.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file $OUTDIR/idr/true_reps.idr \
    --idr-threshold 0.05 --plot \
    --log-output-file $OUTDIR/idr/true_reps.log

idr --samples \
        $OUTDIR/psr1_1/rep1_psr1_peaks.sorted.narrowPeak \
        $OUTDIR/psr1_2/rep1_psr2_peaks.sorted.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file $OUTDIR/idr/rep1_pseudoreps.idr \
    --idr-threshold 0.10 --plot \
    --log-output-file $OUTDIR/idr/rep1_psr.log

# 5. Filter against ENCODE blacklist
bedtools intersect -v -a $OUTDIR/idr/true_reps.idr -b $BLACKLIST \
    > $OUTDIR/idr/true_reps.no_blacklist.bed

# 6. Self-consistency rule: max(Nt, Nself) / min(Nt, Nself) <= 2 to pass
NT=$(awk '$5 >= 540' $OUTDIR/idr/true_reps.idr | wc -l)
NSELF=$(awk '$5 >= 415' $OUTDIR/idr/rep1_pseudoreps.idr | wc -l)   # 0.10 IDR -> score >= 415
RATIO=$(awk -v a=$NT -v b=$NSELF 'BEGIN{m=(a>b)?a:b; n=(a<b)?a:b; print (n>0)? m/n : "inf"}')
echo "Nt=$NT  Nself=$NSELF  ratio=$RATIO  (ENCODE pass: ratio <= 2)"
