#!/bin/bash
# Reference: kmc 3.2+, genomescope2 2.0+, merqury 1.3+ (best_k.sh), smudgeplot 0.2.5+ | Verify API if version differs
# Pre-assembly genome profiling: count k-mers from accurate reads, model the spectrum
# with GenomeScope2 (size, heterozygosity, repeat content, ploidy), and infer ploidy
# independently with Smudgeplot. Run on Illumina or PacBio HiFi reads only -- never noisy ONT.
set -euo pipefail

READS_LIST=fastq_list.txt   # one fastq(.gz) path per line; accurate reads only
OUTPREFIX=sample
THREADS=16
MEM_GB=64
EXPECTED_SIZE=1000000000    # ~1 Gb; only used to derive k via best_k.sh
PLOIDY=2                    # diploid default; raise if Smudgeplot/biology says polyploid

# --- choose k from the expected genome size ---
# best_k.sh solves k = log4(G(1-p)/p) with tolerable collision rate p=0.001 (Illumina WGS default);
# ~1 Gb -> k=20 (~3 Gb -> k=21, the common GenomeScope default). The same k is used for counting and modeling.
K=$(sh "${MERQURY}/best_k.sh" "${EXPECTED_SIZE}" | tail -1 | awk '{print int($1+0.5)}')
echo "best_k.sh chose k=${K}"

# --- count k-mers (KMC) ---
# -ci1 keeps singleton k-mers (the error shoulder GenomeScope models); dropping them breaks the fit.
# -cs10000 / -cx10000 cap the histogram so a single organelle/high-copy-repeat spike at very high
# multiplicity does not dominate the file. KMC counts canonical k-mers by default (WGS is unstranded).
mkdir -p kmc_tmp
kmc -k"${K}" -t"${THREADS}" -m"${MEM_GB}" -ci1 -cs10000 "@${READS_LIST}" kmc_db kmc_tmp/
kmc_tools transform kmc_db histogram "${OUTPREFIX}.histo" -cx10000

# --- model the spectrum (GenomeScope2) ---
# Reads out: haploid genome length, heterozygosity rate, % unique (inverse of repeat content),
# and the homozygous-peak k-mer coverage (kmercov/lambda) -- which is per-k-mer, NOT per-base depth.
genomescope2 -i "${OUTPREFIX}.histo" -o "${OUTPREFIX}_genomescope" -k "${K}" -p "${PLOIDY}"
echo "GenomeScope2 summary:"
cat "${OUTPREFIX}_genomescope/summary.txt"

# --- infer ploidy independently (Smudgeplot, classic KMC backend) ---
# L ~0.5x and U ~8.5x the haploid k-mer peak: below L is error, above U is high-copy repeat.
# Smudgeplot classifies het k-mer PAIRS by CovB/(CovA+CovB): ~0.5 diploid AB, ~0.33 triploid AAB,
# ~0.25 tetraploid AABB. Disagreement with GenomeScope's -p is itself informative.
L=$(smudgeplot.py cutoff "${OUTPREFIX}.histo" L)
U=$(smudgeplot.py cutoff "${OUTPREFIX}.histo" U)
echo "Smudgeplot coverage cutoffs: L=${L} U=${U}"
kmc_tools transform kmc_db -ci"${L}" -cx"${U}" dump -s "kmc_L${L}_U${U}.dump"
smudgeplot.py hetkmers -o "${OUTPREFIX}_pairs" < "kmc_L${L}_U${U}.dump"
smudgeplot.py plot "${OUTPREFIX}_pairs_coverages.tsv" -o "${OUTPREFIX}"

echo "Profiling done. Use the GenomeScope haploid size as the NG50 denominator (assembly-qc),"
echo "Flye -g, and hifiasm --hom-cov. An assembly 1.5-2x this size = uncollapsed haplotigs."
