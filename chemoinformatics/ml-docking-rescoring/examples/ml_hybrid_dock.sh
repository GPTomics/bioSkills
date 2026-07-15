#!/usr/bin/env bash
# Reference: DiffDock-L, GNINA 1.1+, PoseBusters 0.6+ | Verify API if version differs
# Hybrid VS: DiffDock-L pose sampling + GNINA CNN rescore + PoseBusters QC

set -euo pipefail

RECEPTOR_PDB="${1:-receptor.pdb}"
COMPLEX_CSV="${2:-complexes.csv}"
DIFFDOCK_DIR="${3:-/path/to/DiffDock}"
OUT_DIR="${4:-out}"

mkdir -p "${OUT_DIR}/diffdock" "${OUT_DIR}/rescored" "${OUT_DIR}/validated"

# This reference script assumes every row in complexes.csv uses RECEPTOR_PDB.
# For multi-receptor batches, resolve each complex's protein_path separately.

# DiffDock is run from its checkout because the official configuration and model
# paths are relative to that directory. Resolve caller paths before changing cwd.
COMPLEX_CSV="$(cd "$(dirname "${COMPLEX_CSV}")" && pwd -P)/$(basename "${COMPLEX_CSV}")"
OUT_DIR="$(cd "${OUT_DIR}" && pwd -P)"

# Inference runs from DIFFDOCK_DIR, so file paths stored inside the CSV must be
# absolute. Ligand descriptions that are SMILES are left unchanged; ligand file
# paths, like protein paths, must be absolute to avoid cwd-dependent failures.
COMPLEX_CSV="${COMPLEX_CSV}" python - <<'PY'
import csv
import os
from pathlib import Path

csv_path = Path(os.environ['COMPLEX_CSV'])
file_endings = ('.sdf', '.sdf.gz', '.mol', '.mol2', '.pdb', '.pdbqt')
with csv_path.open(newline='') as handle:
    for row_number, row in enumerate(csv.DictReader(handle), start=2):
        protein_path = row.get('protein_path', '').strip()
        if protein_path and not Path(protein_path).is_absolute():
            raise SystemExit(
                f'{csv_path}:{row_number}: protein_path must be absolute: '
                f'{protein_path!r}'
            )
        ligand = row.get('ligand_description', '').strip()
        if (ligand and ligand.lower().endswith(file_endings)
                and not Path(ligand).is_absolute()):
            raise SystemExit(
                f'{csv_path}:{row_number}: ligand file path must be absolute: '
                f'{ligand!r}'
            )
PY

# Step 1: DiffDock-L pose sampling. complexes.csv follows the official
# protein_ligand_csv schema, including complex_name, protein_path, and
# ligand_description columns. DiffDock does not install diffdock_inference.
# Forty samples and 20 inference steps are repository starting defaults; tune
# them against pose quality, runtime, and the selected DiffDock checkpoint.
(
    cd "${DIFFDOCK_DIR}"
    python -m inference \
        --config default_inference_args.yaml \
        --protein_ligand_csv "${COMPLEX_CSV}" \
        --out_dir "${OUT_DIR}/diffdock" \
        --samples_per_complex 40 \
        --inference_steps 20
)

# Step 2: GNINA CNN rescoring on all nested DiffDock poses. `find -print0`
# supports spaces and works with the macOS Bash 3 shipped on many systems.
find "${OUT_DIR}/diffdock" -type f -name '*.sdf' -print0 |
while IFS= read -r -d '' sdf; do
    rel="${sdf#"${OUT_DIR}/diffdock/"}"
    stem=$(basename "${rel%.sdf}")
    digest=$(printf '%s' "${rel}" | shasum -a 256 | awk '{print $1}')
    name="${stem}__${digest}"
    gnina \
        -r "${RECEPTOR_PDB}" \
        -l "${sdf}" \
        --cnn_scoring rescore \
        -o "${OUT_DIR}/rescored/${name}_rescored.sdf" \
        --score_only
done

# Step 3: PoseBusters physical validation. Hashing the relative DiffDock path
# above prevents basename collisions when nested complex directories repeat.
find "${OUT_DIR}/rescored" -type f -name '*.sdf' -print0 |
while IFS= read -r -d '' sdf; do
    name=$(basename "${sdf}" .sdf)
    bust "${sdf}" -p "${RECEPTOR_PDB}" --outfmt=csv \
        > "${OUT_DIR}/validated/${name}_pb.csv"
done

# Step 4: Combine and filter to PB-valid + ranked
OUT_DIR="${OUT_DIR}" python3 <<'PY'
import os
import pandas as pd

out_dir = os.environ['OUT_DIR']
rows = []
validated_dir = os.path.join(out_dir, 'validated')
for f in os.listdir(validated_dir):
    if not f.endswith('_pb.csv'):
        continue
    df = pd.read_csv(os.path.join(validated_dir, f))
    bool_cols = df.select_dtypes(include='bool').columns
    df['pb_valid'] = df[bool_cols].all(axis=1)
    df['ligand'] = f.replace('_rescored_pb.csv', '')
    rows.append(df)

combined = pd.concat(rows) if rows else pd.DataFrame()
valid = combined[combined['pb_valid']] if not combined.empty else combined
print(f'Total poses: {len(combined)}; PB-valid: {len(valid)}')
PY

echo "Output in ${OUT_DIR}/validated/"
