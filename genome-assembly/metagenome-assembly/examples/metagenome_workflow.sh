#!/bin/bash
# Reference: flye 2.9+, metabat2 2.15+, semibin 2.0+, das_tool 1.1.6+, checkm2 1.0+, gunc 1.0+, gtdbtk 2.3+, minimap2 2.26+, samtools 1.19+ | Verify API if version differs
# Metagenome -> MAGs: assemble (--meta), multi-binner -> DAS_Tool consolidation, CheckM2 + GUNC QC, GTDB-Tk taxonomy.
# The deliverable is a set of MIMAG-gated MAGs, not a single assembly N50.

set -euo pipefail

SAMPLE="community"
THREADS=32
# Multiple per-sample BAMs drive DIFFERENTIAL-COVERAGE binning, the strongest signal.
# A single sample collapses binning to weak composition-only TNF; list all samples here.
READS_LIST=("sample1.fastq.gz" "sample2.fastq.gz" "sample3.fastq.gz")
ASSEMBLY_READS="${READS_LIST[0]}"   # assemble from one (or co-assemble) then map ALL back
MIN_CONTIG=1500   # MetaBAT2 default; contigs below ~1000 bp have unreliable TNF/coverage -> chimeras

echo "=== Step 1: metaFlye assembly (--meta mandatory for communities) ==="
flye --meta --nano-hq "${ASSEMBLY_READS}" --out-dir "${SAMPLE}_flye" -t "${THREADS}"
ASSEMBLY="${SAMPLE}_flye/assembly.fasta"
seqkit stats "${ASSEMBLY}"

echo "=== Step 2: Map every sample back for differential coverage ==="
BAMS=()
for reads in "${READS_LIST[@]}"; do
    bam="${SAMPLE}_$(basename "${reads}" .fastq.gz).sorted.bam"
    minimap2 -ax map-ont -t "${THREADS}" "${ASSEMBLY}" "${reads}" | \
        samtools sort -@ "${THREADS}" -o "${bam}" -
    samtools index "${bam}"
    BAMS+=("${bam}")
done
jgi_summarize_bam_contig_depths --outputDepth "${SAMPLE}_depth.txt" "${BAMS[@]}"

echo "=== Step 3: Run MULTIPLE binners (never just one; each recovers a different genome set) ==="
mkdir -p metabat_bins semibin_out
metabat2 -i "${ASSEMBLY}" -a "${SAMPLE}_depth.txt" -o metabat_bins/bin -m "${MIN_CONTIG}" -t "${THREADS}"
SemiBin2 single_easy_bin -i "${ASSEMBLY}" -b "${BAMS[@]}" -o semibin_out

echo "=== Step 4: Consolidate with DAS_Tool (ensemble beats any single binner) ==="
# DAS_Tool assumes the SAME contig set across binners; do not mix assemblies.
Fasta_to_Contig2Bin.sh -i metabat_bins/ -e fa > metabat.tsv
Fasta_to_Contig2Bin.sh -i semibin_out/output_bins/ -e fa.gz > semibin.tsv
DAS_Tool -i metabat.tsv,semibin.tsv -l metabat,semibin \
    -c "${ASSEMBLY}" -o dastool/DAS --write_bins -t "${THREADS}"
BINS="dastool/DAS_DASTool_bins"

echo "=== Step 5: QC every MAG with CheckM2 AND GUNC (completeness lies about chimeras) ==="
checkm2 predict --input "${BINS}" -x fa --output-directory "${SAMPLE}_checkm2" --threads "${THREADS}"
gunc run --input_dir "${BINS}" --file_suffix .fa --out_dir "${SAMPLE}_gunc" --threads "${THREADS}"

echo "=== Step 6: Identify complete circular genomes (Flye flags them) ==="
# assembly_info.txt has a circularity column ('circ.'); HiFi recovers circular MAGs short reads cannot.
awk 'NR>1 && $4=="Y"{print $1}' "${SAMPLE}_flye/assembly_info.txt" > circular_contigs.txt
echo "Circular contigs: $(wc -l < circular_contigs.txt)"

echo "=== Step 7: Taxonomy with GTDB-Tk (DB release must match the binary) ==="
# gtdbtk classify_wf --genome_dir "${BINS}" -x fa --out_dir "${SAMPLE}_gtdbtk" --cpus "${THREADS}"

echo "=== Done. Report MAG count by MIMAG tier, not assembly N50 ==="
echo "MAGs: ${BINS}/  | CheckM2: ${SAMPLE}_checkm2/quality_report.tsv  | GUNC: ${SAMPLE}_gunc/"
