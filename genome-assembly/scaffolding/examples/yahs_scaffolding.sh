#!/bin/bash
# Reference: yahs 1.2+, bwa 0.7.17+, samtools 1.19+, seqkit 2.6+, tidk 0.2.3+, juicer_tools 1.22/2.0+ | Verify API if version differs
# Hi-C scaffolding from contigs to a CURATABLE draft. The output is a draft chromosome
# structure, not a finished genome - load the .hic + .assembly in Juicebox/PretextView and
# curate (break misjoins, flip inversions, remove false duplications) before trusting it.
set -euo pipefail

CONTIGS=$1                    # gap-free contigs (e.g. hifiasm/Flye output), haplotigs already purged
HIC_R1=$2
HIC_R2=$3
OUT=${4:-out}
THREADS=${5:-16}

ENZYME=${ENZYME:-GATC}        # GATC=DpnII/Arima; set ENZYME='' (empty) for Omni-C/enzyme-free Hi-C
MIN_MAPQ=10                   # drop repeat-ambiguous Hi-C reads that create spurious joins
TELO_MOTIF=${TELO_MOTIF:-TTAGGG}   # vertebrate telomere repeat; change per clade (tidk find lists clade motifs)
JUICER_MEM=48                 # GB heap for juicer_tools pre; scale to host RAM

echo "[1/6] Indexing contigs..."
samtools faidx "$CONTIGS"
bwa index "$CONTIGS"

# Hi-C reads are NOT a normal PE library: the two ends are different loci ligated together.
# Map each end SEPARATELY with no mate rescue. -5 = 5' (junction) portion primary; -S/-P skip rescue/pairing.
echo "[2/6] Mapping Hi-C (each end separately, no mate rescue)..."
bwa mem -5SP -T0 -t "$THREADS" "$CONTIGS" "$HIC_R1" "$HIC_R2" | \
    samtools view -@ 8 -b - > aligned.bam

# Hi-C is duplicate-rich; mark/remove PCR + optical duplicates before scaffolding (mandatory).
echo "[3/6] Marking duplicates..."
samtools sort -@ "$THREADS" -n aligned.bam | \
    samtools fixmate -@ "$THREADS" -m - - | \
    samtools sort -@ "$THREADS" - | \
    samtools markdup -@ "$THREADS" - "${OUT}_hic.bam"
samtools index "${OUT}_hic.bam"

echo "[4/6] Running YaHS (contig error-correction ON by default - breaks chimeras)..."
if [ -n "$ENZYME" ]; then
    yahs -e "$ENZYME" -q "$MIN_MAPQ" --telo-motif "$TELO_MOTIF" -o "$OUT" "$CONTIGS" "${OUT}_hic.bam"
else
    yahs -q "$MIN_MAPQ" --telo-motif "$TELO_MOTIF" -o "$OUT" "$CONTIGS" "${OUT}_hic.bam"
fi
SCAFFOLDS="${OUT}_scaffolds_final.fa"

# Prepare .hic + .assembly for Juicebox Assembly Tools (JBAT) - this is the curation handoff.
echo "[5/6] Preparing contact map for Juicebox curation..."
juicer pre -a -o "${OUT}_JBAT" "${OUT}.bin" "${OUT}_scaffolds_final.agp" "${CONTIGS}.fai" 2>"${OUT}_jbat.log"
java "-Xmx${JUICER_MEM}G" -jar juicer_tools.jar pre "${OUT}_JBAT.txt" "${OUT}_JBAT.hic" \
    <(grep PRE_C_SIZE "${OUT}_jbat.log" | awk '{print $2" "$3}')
echo "    Load ${OUT}_JBAT.hic + ${OUT}_JBAT.assembly in Juicebox, curate, export ${OUT}_JBAT.review.assembly,"
echo "    then: juicer post -o ${OUT}_curated ${OUT}_JBAT.review.assembly ${OUT}_JBAT.liftover.agp ${CONTIGS}"

# QC: report contig N50 and scaffold N50 SEPARATELY (scaffolding adds N, not sequence), plus telomeres.
echo "[6/6] QC: contig vs scaffold contiguity + telomeres..."
{ echo "## contigs:"; seqkit stats -a "$CONTIGS"; echo "## scaffolds:"; seqkit stats -a "$SCAFFOLDS"; } > "${OUT}_stats.txt"
tidk search -s "$TELO_MOTIF" -o "${OUT}_telo" -d . "$SCAFFOLDS" || true

echo "Done. Scaffolds=$SCAFFOLDS (DRAFT - curate before trusting chromosome structure)"
echo "  Contact map: ${OUT}_JBAT.hic | Stats (contig vs scaffold N50): ${OUT}_stats.txt"
