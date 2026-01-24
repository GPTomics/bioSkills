---
name: bio-proteomics-ptm-analysis
description: Post-translational modification analysis including phosphorylation, acetylation, and ubiquitination. Covers site localization, motif analysis, and quantitative PTM analysis. Use when analyzing phosphoproteomic data or other modification-enriched samples.
tool_type: mixed
primary_tool: pyOpenMS
---

# Post-Translational Modification Analysis

## Common PTMs and Mass Shifts

```python
PTM_MASSES = {
    'Phosphorylation': 79.966331,      # STY
    'Oxidation': 15.994915,             # M
    'Acetylation': 42.010565,           # K, N-term
    'Methylation': 14.015650,           # KR
    'Dimethylation': 28.031300,         # KR
    'Trimethylation': 42.046950,        # K
    'Ubiquitination': 114.042927,       # K (GlyGly remnant)
    'Deamidation': 0.984016,            # NQ
    'Carbamidomethyl': 57.021464,       # C (fixed mod from IAA)
}
```

## Processing MaxQuant PTM Output

```python
import pandas as pd
import numpy as np

# Phospho(STY)Sites.txt from MaxQuant
phospho = pd.read_csv('Phospho (STY)Sites.txt', sep='\t', low_memory=False)

# Filter valid sites
phospho = phospho[
    (phospho['Reverse'] != '+') &
    (phospho['Potential contaminant'] != '+')
]

# Filter by localization probability
phospho_confident = phospho[phospho['Localization prob'] >= 0.75]
print(f'Confident sites (prob >= 0.75): {len(phospho_confident)}')

# Extract site information
phospho_confident['site'] = phospho_confident.apply(
    lambda r: f"{r['Gene names']}_{r['Amino acid']}{r['Position']}", axis=1
)
```

## Site Localization Scoring

```python
def calculate_ascore_simple(peak_matches_with_ptm, peak_matches_without_ptm, total_peaks):
    '''Simplified A-score calculation'''
    if peak_matches_without_ptm >= peak_matches_with_ptm:
        return 0
    p = peak_matches_with_ptm / total_peaks if total_peaks > 0 else 0
    if p <= 0 or p >= 1:
        return 0

    from scipy.stats import binom
    p_value = 1 - binom.cdf(peak_matches_with_ptm - 1, total_peaks, 0.5)
    return -10 * np.log10(p_value) if p_value > 0 else 100
```

## Motif Analysis

```python
from collections import Counter

def extract_motifs(sites_df, sequence_col, position_col, window=7):
    '''Extract sequence windows around modification sites'''
    motifs = []
    for _, row in sites_df.iterrows():
        seq = row[sequence_col]
        pos = row[position_col] - 1  # 0-indexed
        start = max(0, pos - window)
        end = min(len(seq), pos + window + 1)

        # Pad if at sequence boundary
        motif = '_' * (window - (pos - start)) + seq[start:end] + '_' * (window - (end - pos - 1))
        motifs.append(motif)

    return motifs

def count_amino_acids_by_position(motifs, center=7):
    '''Count amino acid frequencies by position'''
    position_counts = {i: Counter() for i in range(-center, center + 1)}
    for motif in motifs:
        for i, aa in enumerate(motif):
            position_counts[i - center][aa] += 1
    return position_counts
```

## R: Site-Level Quantification with MSstatsPTM

```r
library(MSstatsPTM)

# Prepare input from MaxQuant
ptm_input <- MaxQtoMSstatsPTMFormat(
    evidence = read.table('evidence.txt', sep = '\t', header = TRUE),
    annotation = read.csv('annotation.csv'),
    fasta = 'uniprot_human.fasta',
    mod_type = 'Phospho'
)

# Process data
processed_ptm <- dataSummarizationPTM(ptm_input, method = 'msstats')

# Differential PTM analysis (adjusting for protein-level changes)
ptm_results <- groupComparisonPTM(processed_ptm, contrast.matrix = comparison_matrix)
```

## Related Skills

- peptide-identification - Identify modified peptides
- quantification - Quantify PTM sites
- pathway-analysis/go-enrichment - Enrichment of modified proteins
