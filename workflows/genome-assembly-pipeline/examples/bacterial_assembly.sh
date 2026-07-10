#!/bin/bash
# Reference: genomescope2 2.0+, meryl 1.4+, merqury 1.3+, fastp 0.23+, flye 2.9+, medaka 2.0+, fcs-gx 0.5+, quast 5.2+, busco 5.5+ | Verify API if version differs
# End-to-end single-organism (bacterial isolate) assembly: profile -> QC -> assemble -> polish -> decontaminate -> three-axis QC
set -e

THREADS=16
LONG_READS="nanopore.fastq.gz"           # ONT R10 / Dorado-SUP basecalled
SHORT_R1="illumina_R1.fastq.gz"
SHORT_R2="illumina_R2.fastq.gz"
# shellcheck disable=SC2034  # consumed only when swapping medaka's --bacteria for -m (see Step 3)
BASECALLER_MODEL="r1041_e82_400bps_sup_v5.0.0"   # generic fallback for `-m` if medaka <2.0 lacks --bacteria; MUST match the model the reads were basecalled with
TAXID=562                                 # NCBI tax-id of the organism, for FCS-GX foreign screen
GXDB="/path/to/gxdb"                       # dir holding the ~470 GiB FCS-GX database (gxdb); resident in RAM
KMER=21                                   # common GenomeScope2 default; best_k.sh derives a smaller k for Mb genomes (harmless here, the QV step recomputes)
BUSCO_LINEAGE="bacteria_odb10"
OUTDIR="bacterial_assembly"

# NOTE: 'busco' is intentionally NOT pre-created -- BUSCO aborts if its --out_path/-o dir already exists (needs -f).
mkdir -p ${OUTDIR}/{profile,qc,flye,medaka,fcsgx,final,quast,merqury}

echo "=== Step 0: Profile the genome (size, heterozygosity, ploidy) ==="
# Count from ACCURATE Illumina reads, never noisy ONT (error k-mers swamp the spectrum and break the fit).
meryl count k=${KMER} output ${OUTDIR}/profile/reads.meryl ${SHORT_R1} ${SHORT_R2}
meryl histogram ${OUTDIR}/profile/reads.meryl > ${OUTDIR}/profile/reads.hist
genomescope2 -i ${OUTDIR}/profile/reads.hist -o ${OUTDIR}/profile/gscope -k ${KMER}
echo "Read the estimated genome size from ${OUTDIR}/profile/gscope/ -- it sets the NG50 denominator below."
GENOME_SIZE="5m"      # Flye-style suffix size for --genome-size; replace with the GenomeScope2 estimate
GENOME_BP=5000000     # same size as INTEGER bp for QUAST --est-ref-size and best_k.sh

echo "=== Step 1: Read QC ==="
NanoPlot --fastq ${LONG_READS} --outdir ${OUTDIR}/qc/nanoplot -t ${THREADS}
fastp -i ${SHORT_R1} -I ${SHORT_R2} \
    -o ${OUTDIR}/qc/trimmed_R1.fq.gz -O ${OUTDIR}/qc/trimmed_R2.fq.gz \
    --detect_adapter_for_pe --html ${OUTDIR}/qc/fastp.html -w ${THREADS}

echo "=== Step 2: Assembly with Flye (--nano-hq is the R10 default; --nano-raw on R10 silently collapses repeats) ==="
flye --nano-hq ${LONG_READS} --out-dir ${OUTDIR}/flye --threads ${THREADS} --genome-size ${GENOME_SIZE}
head -20 ${OUTDIR}/flye/assembly_info.txt

echo "=== Step 3: Polish with medaka (model-matched); stop at the Merqury QV plateau, do not over-iterate ==="
medaka_consensus -i ${LONG_READS} -d ${OUTDIR}/flye/assembly.fasta \
    -o ${OUTDIR}/medaka -t ${THREADS} --bacteria    # methylation-aware bacterial model (medaka 2.0+); beats the generic -m ${BASECALLER_MODEL} on isolates
cp ${OUTDIR}/medaka/consensus.fasta ${OUTDIR}/final/assembly.fasta

echo "=== Step 4: Decontaminate -- FCS-GX foreign screen (single-organism path; GenBank-mandatory) ==="
# Cross-kingdom foreign screen. Keep host-integrated foreign sequence; act on near-complete foreign contigs.
# fcs.py screen genome writes an action report; verify flags against the installed FCS-GX version.
python3 ./fcs.py screen genome --fasta ${OUTDIR}/final/assembly.fasta --out-dir ${OUTDIR}/fcsgx \
    --gx-db "${GXDB}/gxdb" --tax-id ${TAXID} || echo "FCS-GX needs the ~470 GiB GX DB resident in RAM; run on a high-memory or cloud host."

echo "=== Step 5: (scaffolding skipped -- no Hi-C for a single bacterial chromosome) ==="

echo "=== Step 6: Three-axis QC -- contiguity + completeness + correctness ==="
# Axis 1: contiguity against the profiled size (NG50/auN, not bare N50)
quast.py ${OUTDIR}/final/assembly.fasta -o ${OUTDIR}/quast -t ${THREADS} --est-ref-size ${GENOME_BP}

# Axis 2: completeness on the deepest applicable clade
busco -i ${OUTDIR}/final/assembly.fasta -l ${BUSCO_LINEAGE} \
    -o busco -m genome -c ${THREADS} --out_path ${OUTDIR}

# Axis 3: correctness -- Merqury QV from ACCURATE reads (k from best_k.sh, not hardcoded)
QV_KMER=$(sh $MERQURY/best_k.sh ${GENOME_BP} | tail -n1 | awk '{print int($1+0.5)}')   # round float->int
meryl count k=${QV_KMER} output ${OUTDIR}/merqury/illumina.meryl \
    ${OUTDIR}/qc/trimmed_R1.fq.gz ${OUTDIR}/qc/trimmed_R2.fq.gz
( cd ${OUTDIR}/merqury && merqury.sh illumina.meryl ../final/assembly.fasta merqury_out )

echo ""
echo "=== Assembly Complete -- report the TRIAD, never N50 alone ==="
echo "Final assembly: ${OUTDIR}/final/assembly.fasta"
grep -E "contigs|Total length|N50|auN|NG50|L50" ${OUTDIR}/quast/report.txt || true
grep -E "Complete|Fragmented|Missing" ${OUTDIR}/busco/short_summary*.txt || true
echo "Merqury QV (correctness): ${OUTDIR}/merqury/merqury_out.qv"
