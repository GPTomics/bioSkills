#!/usr/bin/env bash
# Reference: Panaroo 1.5.1+, PPanGGOLiN 2.2+, Bakta 1.10+, BUSCO 5.7+, Python 3.10+, pandas 2.2+ | Verify API if version differs
# Bacterial pangenome construction with Panaroo + PPanGGOLiN cross-validation.

set -euo pipefail

GENOMES_DIR=${1:?usage: $0 GENOMES_DIR OUTPUT_DIR}
OUTPUT_DIR=${2:?missing output dir}
THREADS=${THREADS:-16}
GENUS=${GENUS:-Escherichia}
SPECIES=${SPECIES:-coli}

mkdir -p "$OUTPUT_DIR"/{annotated,panaroo,ppanggolin}

echo "[1/5] Annotate all genomes with Bakta"
for fa in "$GENOMES_DIR"/*.fa; do
    name=$(basename "$fa" .fa)
    if [ ! -f "$OUTPUT_DIR/annotated/${name}/${name}.gff3" ]; then
        bakta --db "${BAKTA_DB:-/path/to/bakta-db}" --threads "$THREADS" \
            --output "$OUTPUT_DIR/annotated/${name}" --prefix "$name" \
            --genus "$GENUS" --species "$SPECIES" \
            "$fa"
    fi
done

echo "[2/5] Run Panaroo in strict mode"
panaroo -i "$OUTPUT_DIR/annotated"/*/*.gff3 \
    -o "$OUTPUT_DIR/panaroo" -t "$THREADS" \
    --clean-mode strict --remove-invalid-genes

echo "[3/5] Extract Tettelin partition (core/shell/cloud)"
python3 - <<PY
import pandas as pd
import numpy as np

presence = pd.read_csv('$OUTPUT_DIR/panaroo/gene_presence_absence_roary.csv', sep=',')
# Long-format strain columns are at column index 14 onward in Panaroo output
strain_cols = [c for c in presence.columns if c not in ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']]
n_strains = len(strain_cols)
print(f'Strains: {n_strains}')

# Compute presence fraction per gene
mask = presence[strain_cols].notna() & (presence[strain_cols] != '')
fraction = mask.sum(axis=1) / n_strains

# Tettelin classification at 95%/15% thresholds
def classify(f):
    if f >= 0.95: return 'core'
    if f >= 0.15: return 'shell'
    return 'cloud'

presence['fraction'] = fraction
presence['class'] = fraction.apply(classify)
counts = presence['class'].value_counts()
print('\nTettelin partition:')
print(counts)
print(f'\nTotal pangene families: {len(presence)}')
print(f'Core (>=95%):  {counts.get("core", 0)}')
print(f'Shell (15-95%): {counts.get("shell", 0)}')
print(f'Cloud (<15%):   {counts.get("cloud", 0)}')
PY

echo "[4/5] PPanGGOLiN HMM partition (cross-validation)"
# PPanGGOLiN --anno expects a TSV index (name<TAB>path), not a glob.
for f in "$OUTPUT_DIR/annotated"/*/*.gff3; do
    printf "%s\t%s\n" "$(basename "$f" .gff3)" "$(realpath "$f")"
done > "$OUTPUT_DIR/ppanggolin/gff_index.tsv"

ppanggolin all --anno "$OUTPUT_DIR/ppanggolin/gff_index.tsv" \
    -o "$OUTPUT_DIR/ppanggolin" -c "$THREADS"

echo "[5/5] Heaps law fitting"
python3 - <<PY
import pandas as pd
import numpy as np

presence = pd.read_csv('$OUTPUT_DIR/panaroo/gene_presence_absence_roary.csv', sep=',')
strain_cols = [c for c in presence.columns if c not in ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']]
mat = presence[strain_cols].notna().astype(int).values
n_strains = mat.shape[1]

# Resample 100 orderings; report pangenome growth
alphas = []
np.random.seed(42)
for _ in range(100):
    perm = np.random.permutation(n_strains)
    pan = set()
    sizes = []
    for i in perm:
        present = np.where(mat[:, i] == 1)[0]
        pan.update(present.tolist())
        sizes.append(len(pan))
    log_n = np.log(np.arange(1, n_strains + 1))
    log_pan = np.log(sizes)
    alpha = np.polyfit(log_n, log_pan, 1)[0]
    alphas.append(alpha)

print(f'\nHeaps law alpha: {np.mean(alphas):.3f} (95% CI [{np.percentile(alphas, 2.5):.3f}, {np.percentile(alphas, 97.5):.3f}])')
print(f'Pangenome is {"OPEN" if np.mean(alphas) < 1 else "CLOSED"}')
PY

echo "Done. Panaroo: $OUTPUT_DIR/panaroo; PPanGGOLiN: $OUTPUT_DIR/ppanggolin"
