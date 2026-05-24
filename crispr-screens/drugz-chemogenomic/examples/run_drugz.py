# Reference: drugZ Aug-2019+ (hart-lab/drugz), pandas 2.2+, numpy 1.26+ | Verify API if version differs
#
# drugZ chemogenomic screen analysis end-to-end.
# Compares drug-treated vs vehicle-treated arms at matched timepoint.

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

# === INPUTS ===
counts_file = 'counts.txt'                          # sgRNA, GENE, sample columns
vehicle_samples = 'Veh_r1,Veh_r2,Veh_r3'
drug_samples = 'Drug_r1,Drug_r2,Drug_r3'
output_dir = Path('drugz_output')
output_dir.mkdir(exist_ok=True)

# === DOWNLOAD CEGv2 (for excluding essentials from null) ===
ceg_file = Path('CEGv2.txt')
if not ceg_file.exists():
    subprocess.run(['curl', '-L', '-o', str(ceg_file),
                    'https://raw.githubusercontent.com/hart-lab/bagel/master/CEGv2.txt'],
                   check=True)

# === STEP 1: STANDARD drugZ RUN ===
output_file = output_dir / 'drugz_standard.txt'
subprocess.run(['python', 'drugz.py',
                '-i', counts_file,
                '-o', str(output_file),
                '-c', vehicle_samples,
                '-x', drug_samples,
                '-p', '5'], check=True)

# === STEP 2: drugZ EXCLUDING ESSENTIALS FROM NULL ===
output_file_clean = output_dir / 'drugz_ceg_excluded.txt'
subprocess.run(['python', 'drugz.py',
                '-i', counts_file,
                '-o', str(output_file_clean),
                '-c', vehicle_samples,
                '-x', drug_samples,
                '-r', str(ceg_file),
                '-p', '5'], check=True)

# === STEP 3: INTERPRET ===
df = pd.read_csv(output_file, sep='\t')

# fdr_synth = sensitizer (negative Z; gene KO sensitizes to drug)
# fdr_supp = suppressor (positive Z; gene KO confers resistance)
sensitizers = df[df['fdr_synth'] < 0.05].sort_values('normZ').head(50)
suppressors = df[df['fdr_supp'] < 0.05].sort_values('normZ', ascending=False).head(50)

# Save tier-stratified output
sensitizers.to_csv(output_dir / 'top50_sensitizers.tsv', sep='\t', index=False)
suppressors.to_csv(output_dir / 'top50_suppressors.tsv', sep='\t', index=False)

# === STEP 4: SUMMARY ===
print(f'Sensitizers (drug + KO synergistic):  {len(sensitizers)}')
print(f'Suppressors (KO confers resistance):  {len(suppressors)}')
print(f'Top 5 sensitizers: {", ".join(sensitizers.head(5)["GENE"].astype(str))}')
print(f'Top 5 suppressors: {", ".join(suppressors.head(5)["GENE"].astype(str))}')

# === STEP 5: HIGH-CONFIDENCE HITS ===
# Apply stringent threshold: fdr_synth < 0.01 AND normZ < -3
strong_sens = df[(df['fdr_synth'] < 0.01) & (df['normZ'] < -3)]
print(f'Strong sensitizers (fdr<0.01, normZ<-3): {len(strong_sens)}')

# === STEP 6: DOSE CONSISTENCY (if multi-dose) ===
# For multi-dose screens, run drugZ at each dose then check consistency
# def per_dose_drugz(doses_dict):
#     # doses_dict: {'low': 'Drug_low_r1,Drug_low_r2', 'mid': 'Drug_mid_r1,Drug_mid_r2', ...}
#     for dose, samples in doses_dict.items():
#         subprocess.run(['python', 'drugz.py', ...])
#     # Then aggregate to find dose-consistent hits
