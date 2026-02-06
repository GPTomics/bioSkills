#!/usr/bin/env bash
## FINEMAP pipeline
##
## Prepares input files and runs FINEMAP for fine-mapping a GWAS locus.
## Requires: FINEMAP binary, plink, GWAS summary stats, reference panel.

set -euo pipefail

# --- Configuration ---
GWAS_FILE="gwas_sumstats.txt"    # Columns: SNP CHR POS A1 A2 MAF BETA SE
REF_PANEL="1000G_EUR"            # Plink binary prefix for LD reference
CHR=6
START=30000000
END=31000000
N_SAMPLES=50000
N_CAUSAL=5                       # Max number of causal SNPs to model

# --- Step 1: Extract locus SNPs ---
awk -v chr="$CHR" -v start="$START" -v end="$END" \
  '$2 == chr && $3 >= start && $3 <= end' "$GWAS_FILE" > locus_gwas.txt

# Extract SNP IDs for plink
awk '{print $1}' locus_gwas.txt > locus_snps.txt

# --- Step 2: Generate LD matrix from reference panel ---
plink --bfile "$REF_PANEL" \
  --chr "$CHR" --from-bp "$START" --to-bp "$END" \
  --extract locus_snps.txt \
  --r square \
  --out locus_ld

# --- Step 3: Create .z file for FINEMAP ---
# Format: rsid chromosome position allele1 allele2 maf beta se
echo "rsid chromosome position allele1 allele2 maf beta se" > locus.z
awk -v OFS=" " '{print $1, $2, $3, $4, $5, $6, $7, $8}' locus_gwas.txt >> locus.z

# Rename LD file
mv locus_ld.ld locus.ld

# --- Step 4: Create master file ---
cat > master.txt << EOF
z;ld;snp;config;cred;log;n_samples
locus.z;locus.ld;locus.snp;locus.config;locus.cred;locus.log;${N_SAMPLES}
EOF

# --- Step 5: Run FINEMAP ---
# --sss: Shotgun stochastic search (default, recommended)
# --n-causal-snps: Max causal SNPs (5 is reasonable for most loci)
finemap --sss --in-files master.txt --n-causal-snps "$N_CAUSAL"

# --- Step 6: Report results ---
echo ""
echo "=== FINEMAP Results ==="
echo "Top 10 variants by PIP:"
head -1 locus.snp
sort -k11 -nr locus.snp | head -10

echo ""
echo "Credible sets:"
cat locus.cred

echo ""
echo "Model configurations:"
head -5 locus.config
