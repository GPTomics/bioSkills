#!/usr/bin/env bash
# Reference: PLINK 2.0 (alpha 6+), numpy 1.26+, pandas 2.2+, scipy 1.12+ | Verify API if version differs
# Single-variant GWAS with PC adjustment and Firth fallback, then a lambda diagnostic that frames
# inflation as polygenicity vs confounding (NOT a license to genomic-control). Outputs (PCA,
# association, plots) go to a caller-supplied directory (default: a fresh temp dir) so nothing is
# written to the current dir. Use a mixed model (SAIGE/regenie/GEMMA/BOLT-LMM) instead of this GLM
# when the sample has relatedness, or case:control imbalance with low MAC.
# Usage: ./gwas_pipeline.sh <plink_prefix> <pheno_file> [output_dir]
set -euo pipefail

BFILE="${1:?Usage: $0 <plink_prefix> <pheno_file> [output_dir]}"
PHENO="${2:?Usage: $0 <plink_prefix> <pheno_file> [output_dir]}"
OUTDIR="${3:-$(mktemp -d)}"
mkdir -p "$OUTDIR"
echo "Outputs -> $OUTDIR"

# PCs to absorb continuous ancestry. Compute on LD-pruned, MAF-filtered genotypes upstream
# (see population-structure); too few PCs leave stratification, too many absorb real signal.
plink2 --bfile "$BFILE" --pca 10 --out "$OUTDIR/pca"

# Single-variant association. firth-fallback is the binary-trait default and writes
# .glm.logistic.hybrid (mixed logistic and Firth rows); a quantitative phenotype writes .glm.linear.
# hide-covar drops per-covariate rows from the output.
plink2 --bfile "$BFILE" \
    --pheno "$PHENO" \
    --covar "$OUTDIR/pca.eigenvec" --covar-name PC1-PC10 \
    --glm firth-fallback hide-covar cols=+a1freq \
    --out "$OUTDIR/gwas"

RESULT_FILE=""
for f in "$OUTDIR"/gwas.*.glm.*; do
    [[ -f "$f" && "$f" != *.log ]] && { RESULT_FILE="$f"; break; }
done
if [[ -z "${RESULT_FILE:-}" || ! -f "$RESULT_FILE" ]]; then
    echo "No association output produced"
    exit 1
fi
echo "Results: $RESULT_FILE"

# Summarize and diagnose. Lambda > 1 under a polygenic trait is EXPECTED and mostly true signal:
# this prints lambda for diagnosis only, it does NOT genomic-control. The confounding diagnostic
# is the LDSC intercept (run ldsc on the sumstats), not lambda.
python3 - "$RESULT_FILE" "$OUTDIR" <<'PY'
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

result_file, outdir = sys.argv[1], sys.argv[2]
GENOME_WIDE = 5e-8   # European common-variant convention; stricter for African ancestry / WGS
SUGGESTIVE = 1e-5

df = pd.read_csv(result_file, sep='\t')
df = df[df['TEST'] == 'ADD'].copy()
df = df.dropna(subset=['P'])
df['neglog10p'] = -np.log10(df['P'])

n_gw = int((df['P'] < GENOME_WIDE).sum())
n_sug = int((df['P'] < SUGGESTIVE).sum())
chisq = stats.chi2.ppf(1 - df['P'].clip(upper=1 - 1e-300), 1)
lam = np.median(chisq) / stats.chi2.ppf(0.5, 1)
print(f'tests={len(df)} genome_wide<{GENOME_WIDE}={n_gw} suggestive<{SUGGESTIVE}={n_sug} lambda={lam:.3f}')
print('lambda > 1 is expected under polygenicity; use the LDSC intercept to judge confounding')

chrom_col = '#CHROM' if '#CHROM' in df.columns else 'CHROM'
df = df.sort_values([chrom_col, 'POS'])
offset, cum, ticks, labels = 0, [], [], []
for chrom, block in df.groupby(chrom_col, sort=False):
    pos = block['POS'].to_numpy() + offset
    cum.extend(pos)
    ticks.append((pos.min() + pos.max()) / 2)
    labels.append(str(chrom))
    offset = pos.max()
df['cum_pos'] = cum

plt.figure(figsize=(14, 5))
colors = np.where(df.groupby(chrom_col).ngroup() % 2 == 0, '#1f77b4', '#ff7f0e')
plt.scatter(df['cum_pos'], df['neglog10p'], c=colors, s=2)
plt.axhline(-np.log10(GENOME_WIDE), color='red', ls='--', lw=0.8)
plt.axhline(-np.log10(SUGGESTIVE), color='blue', ls='--', lw=0.8)
plt.xticks(ticks, labels, fontsize=7)
plt.xlabel('Chromosome'); plt.ylabel('-log10(P)')
plt.tight_layout(); plt.savefig(f'{outdir}/manhattan.png', dpi=150); plt.close()

obs = np.sort(df['P'].to_numpy())
exp = (np.arange(1, len(obs) + 1) - 0.5) / len(obs)
plt.figure(figsize=(5, 5))
plt.scatter(-np.log10(exp), -np.log10(obs), s=2)
hi = max(-np.log10(exp).max(), -np.log10(obs).max())
plt.plot([0, hi], [0, hi], 'r--', lw=0.8)
plt.xlabel('Expected -log10(P)'); plt.ylabel('Observed -log10(P)')
plt.title(f'lambda = {lam:.3f}')
plt.tight_layout(); plt.savefig(f'{outdir}/qqplot.png', dpi=150); plt.close()
print(f'plots -> {outdir}/manhattan.png, {outdir}/qqplot.png')
PY

echo "Done. All outputs in $OUTDIR"
