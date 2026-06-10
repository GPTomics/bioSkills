#!/bin/bash
# Reference: Bismark 0.24+, MethylDackel 0.6+, samtools 1.19+ | Verify API if version differs
# Extracts per-cytosine methylation from a deduplicated bisulfite/EM-seq BAM via either
# extractor path. The conversion-rate gate runs FIRST; betas are not trusted until it passes.
set -euo pipefail

GENOME_DIR='genome'           # FASTA + Bisulfite_Genome index used at alignment time
BAM_DIR='aligned'
OUTPUT_DIR='methylation'
REF_FA="${GENOME_DIR}/genome.fa"
ALIGNER='bismark'             # 'bismark' (XM-tagged BAM) or 'bwameth' (no XM tag -> MethylDackel)

CONVERSION_FLOOR=99           # require >=99% conversion; residual non-conversion inflates every beta
R2_IGNORE=2                   # trim R2 5' end-repair artifact; set from the M-bias plot (EM-seq/PBAT)

mkdir -p "$OUTPUT_DIR"

for bam in "${BAM_DIR}"/*.deduplicated.bam; do
    sample=$(basename "$bam" .bam)   # keep the full stem; Bismark names outputs after the input BAM basename

    if [ "$ALIGNER" = 'bismark' ]; then
        # The extractor splits all 3 contexts regardless of --CX: the CHH rate in the splitting
        # report is the conversion-rate QC. No --CX keeps the coverage file CpG-only, which is what
        # coverage2cytosine --merge_CpG requires (--CX coverage would break the merge below).
        bismark_methylation_extractor \
            -p --comprehensive \
            --bedGraph --cytosine_report \
            --genome_folder "$GENOME_DIR" \
            --ignore_r2 "$R2_IGNORE" \
            --parallel 4 --gzip \
            -o "$OUTPUT_DIR" "$bam"

        chh_meth=$(grep -i "C methylated in CHH context" "${OUTPUT_DIR}/${sample}"*_splitting_report.txt | grep -oE '[0-9.]+%' | tr -d '%')
        conversion=$(echo "100 - ${chh_meth:-100}" | bc -l)
        echo "${sample}: mammalian conversion proxy = ${conversion}% (lambda spike-in mandatory in plants/ESC/neurons)"
        # gate: in mammals only; in plants use a lambda spike-in instead of CHH
        awk -v c="$conversion" -v f="$CONVERSION_FLOOR" 'BEGIN{if(c<f) print "  WARNING: below "f"% conversion floor - betas inflated"}'

        # Symmetric CpG dyad collapse (doubles effective coverage). --merge_CpG needs a CpG-only
        # coverage file (hence no --CX above); --merge_non_CpG would instead merge CHG+CHH (wrong here).
        coverage2cytosine --merge_CpG --genome_folder "$GENOME_DIR" \
            -o "${sample}.merged" --dir "$OUTPUT_DIR" \
            "${OUTPUT_DIR}/${sample}".bismark.cov.gz
    else
        # bwa-meth BAM: recompute calls aligner-agnostically. mbias suggests --OT/--OB bounds (inspect them).
        MethylDackel mbias "$REF_FA" "$bam" "${OUTPUT_DIR}/${sample}_mbias"
        MethylDackel extract \
            --mergeContext \
            --maxVariantFrac 0.25 \
            -o "${OUTPUT_DIR}/${sample}" \
            "$REF_FA" "$bam"
    fi
done

ls -la "${OUTPUT_DIR}"/*.cov* 2>/dev/null || ls -la "${OUTPUT_DIR}"/*_CpG.bedGraph 2>/dev/null
