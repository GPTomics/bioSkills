#!/bin/bash
# Reference: Hostile 1.1+, Bowtie2 2.5+ | Verify API if version differs
# Remove host reads before profiling. Report reads removed - it is a QC metric. Prefer a T2T-CHM13
# index over GRCh38; high-sensitivity Bowtie2 drives removal more than the reference choice.
set -euo pipefail

R1=${1:-sample_R1.fastq.gz}
R2=${2:-sample_R2.fastq.gz}
OUTDIR=${3:-host_depleted}
mkdir -p "$OUTDIR"

# Hostile: removes >99.5% of human reads while discarding far fewer microbial reads than naive mapping.
hostile clean --fastq1 "$R1" --fastq2 "$R2" \
    --index human-t2t-hla --aligner bowtie2 \
    --output "$OUTDIR"
# For long reads use --aligner minimap2. Hostile prints the reads removed; capture it for QC.

# Transparent alignment alternative (keep the UNMAPPED reads). Mask host rDNA in the reference first
# or real microbial reads matching conserved rRNA will be deleted.
# bowtie2 -x t2t_chm13 -1 "$R1" -2 "$R2" --very-sensitive-local -p 8 | \
#     samtools view -b -f 12 -F 256 | samtools fastq -1 "$OUTDIR/clean_R1.fq.gz" -2 "$OUTDIR/clean_R2.fq.gz" -

echo "Host-depleted reads in ${OUTDIR}/ - record how many reads were removed and remaining depth."
