#!/usr/bin/env bash
# Reference: DiffDock-L, GNINA 1.1+, PoseBusters 0.6+ | Verify API if version differs
# Hybrid VS: DiffDock-L pose sampling + GNINA CNN rescore + PoseBusters QC

set -euo pipefail

RECEPTOR_PDB="${1:-receptor.pdb}"
LIGANDS_SMI="${2:-ligands.smi}"
OUT_DIR="${3:-out}"

mkdir -p "${OUT_DIR}/diffdock" "${OUT_DIR}/rescored" "${OUT_DIR}/validated"

# Step 1: DiffDock-L pose sampling
# samples_per_complex=40 is standard for production; 20 for screening
diffdock_inference \
    --protein_path "${RECEPTOR_PDB}" \
    --ligand "${LIGANDS_SMI}" \
    --out_dir "${OUT_DIR}/diffdock" \
    --samples_per_complex 40 \
    --inference_steps 20

# Step 2: GNINA CNN rescoring on all DiffDock poses
# CNN rescore mode: keep DiffDock pose, replace ranking score with GNINA CNN
for sdf in "${OUT_DIR}"/diffdock/*.sdf; do
    name=$(basename "${sdf}" .sdf)
    gnina \
        -r "${RECEPTOR_PDB}" \
        -l "${sdf}" \
        --cnn_scoring rescore \
        -o "${OUT_DIR}/rescored/${name}_rescored.sdf.gz" \
        --score_only
done

# Step 3: PoseBusters validation (mandatory for AI poses; ~50% fail)
for sdf in "${OUT_DIR}"/rescored/*.sdf.gz; do
    name=$(basename "${sdf}" .sdf.gz)
    posebusters bust \
        --mol_pred "${sdf}" \
        --mol_cond "${RECEPTOR_PDB}" \
        --config dock \
        --output "${OUT_DIR}/validated/${name}_pb.csv"
done

# Step 4: Combine and filter to PB-valid + ranked
python3 <<'PY'
import os
import pandas as pd

rows = []
for f in os.listdir('out/validated'):
    if not f.endswith('_pb.csv'):
        continue
    df = pd.read_csv(os.path.join('out/validated', f))
    bool_cols = df.select_dtypes(include='bool').columns
    df['pb_valid'] = df[bool_cols].all(axis=1)
    df['ligand'] = f.replace('_rescored_pb.csv', '')
    rows.append(df)

combined = pd.concat(rows) if rows else pd.DataFrame()
valid = combined[combined['pb_valid']] if not combined.empty else combined
print(f'Total poses: {len(combined)}; PB-valid: {len(valid)}')
PY

echo "Output in ${OUT_DIR}/validated/"
