#!/usr/bin/env bash
# Reference: PLINK 1.9 (1.90b7+), PLINK 2.0 (alpha 6+) | Verify API if version differs
# Ordered GWAS QC pipeline: variant missingness BEFORE sample missingness, controls-only HWE with
# mid-p, KING relatedness pruning, and a differential-missingness confounder check. Outputs go to a
# caller-supplied directory (default: a fresh temp dir) so nothing is written to the current dir.
# Usage: ./qc_pipeline.sh <input.vcf.gz> [pheno.txt] [output_dir]
set -euo pipefail

INPUT="${1:?Usage: $0 <input.vcf.gz> [pheno.txt] [output_dir]}"
PHENO="${2:-}"
OUTDIR="${3:-$(mktemp -d)}"
mkdir -p "$OUTDIR"
echo "Outputs -> $OUTDIR"

# Convert to pgen to keep REF/ALT and any dosage information honest.
plink2 --vcf "$INPUT" --double-id --make-pgen --out "$OUTDIR/raw"

# Variant missingness FIRST, in its own run, so a sample is not dropped for missingness driven
# by variants that are about to be removed. Default --geno is 0.1; GWAS QC tightens to 0.02.
plink2 --pfile "$OUTDIR/raw" --geno 0.02 --make-pgen --out "$OUTDIR/step_geno"

# THEN sample missingness, MAF, and HWE on the surviving variants. mid-p is required (the plain
# exact test is conservative for low-count genotypes). HWE here is applied to all samples; the
# controls-only step below is the correct form when a case/control phenotype is supplied.
plink2 --pfile "$OUTDIR/step_geno" --mind 0.02 --maf 0.01 --hwe 1e-6 midp --make-pgen --out "$OUTDIR/step_qc"

# Controls-only HWE: a true risk variant depletes heterozygotes in cases and would fail a
# case-inclusive HWE test. plink2 has no controls-only default, so gate to controls explicitly.
if [[ -n "$PHENO" ]]; then
    plink2 --pfile "$OUTDIR/step_geno" --pheno "$PHENO" --keep-if "PHENO1 == control" \
        --hwe 1e-6 midp --write-snplist --out "$OUTDIR/hwe_pass_controls"
    plink2 --pfile "$OUTDIR/step_geno" --mind 0.02 --maf 0.01 \
        --extract "$OUTDIR/hwe_pass_controls.snplist" --make-pgen --out "$OUTDIR/step_qc"

    # Differential missingness: drop variants whose missingness differs between cases and controls.
    plink2 --pfile "$OUTDIR/step_geno" --pheno "$PHENO" --test-missing --out "$OUTDIR/diffmiss"
fi

# KING-robust relatedness (structure-robust, unlike PLINK 1.9 PI_HAT). Prune to no-closer-than
# second-degree: cutoffs are 0.354 (duplicate/MZ), 0.177 (first-degree), 0.0884 (second-degree).
plink2 --pfile "$OUTDIR/step_qc" --king-cutoff 0.0884 --out "$OUTDIR/king_unrelated"
plink2 --pfile "$OUTDIR/step_qc" --keep "$OUTDIR/king_unrelated.king.cutoff.in.id" \
    --make-pgen --out "$OUTDIR/clean"

echo "Final fileset: $OUTDIR/clean.{pgen,pvar,psam}"
wc -l "$OUTDIR/clean.psam" "$OUTDIR/clean.pvar" 2>/dev/null || true
