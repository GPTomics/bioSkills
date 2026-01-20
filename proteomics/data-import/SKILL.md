---
name: bio-proteomics-data-import
description: Load and parse mass spectrometry data formats including mzML, mzXML, and quantification tool outputs like MaxQuant proteinGroups.txt. Handles common preprocessing steps like contaminant filtering and missing value assessment.
tool_type: mixed
primary_tool: pyOpenMS
---

# Mass Spectrometry Data Import

## Loading mzML/mzXML Files with pyOpenMS

```python
from pyopenms import MSExperiment, MzMLFile, MzXMLFile

exp = MSExperiment()
MzMLFile().load('sample.mzML', exp)

for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        mz, intensity = spectrum.get_peaks()
    elif spectrum.getMSLevel() == 2:
        precursor = spectrum.getPrecursors()[0]
        precursor_mz = precursor.getMZ()
```

## Loading MaxQuant Output

```python
import pandas as pd

protein_groups = pd.read_csv('proteinGroups.txt', sep='\t', low_memory=False)

# Filter contaminants and reverse hits
protein_groups = protein_groups[
    (protein_groups['Potential contaminant'] != '+') &
    (protein_groups['Reverse'] != '+') &
    (protein_groups['Only identified by site'] != '+')
]

# Extract intensity columns
intensity_cols = [c for c in protein_groups.columns if c.startswith('LFQ intensity')]
intensities = protein_groups[['Protein IDs', 'Gene names'] + intensity_cols]
```

## Loading Spectronaut/DIA-NN Output

```python
diann_report = pd.read_csv('report.tsv', sep='\t')

# Pivot to protein-level matrix
protein_matrix = diann_report.pivot_table(
    index='Protein.Group', columns='Run', values='PG.MaxLFQ', aggfunc='first'
)
```

## R: Loading with MSnbase

```r
library(MSnbase)

raw_data <- readMSData('sample.mzML', mode = 'onDisk')
spectra <- spectra(raw_data)
header_info <- fData(raw_data)
```

## Missing Value Assessment

```python
def assess_missing_values(df, intensity_cols):
    missing_per_protein = df[intensity_cols].isna().sum(axis=1)
    missing_per_sample = df[intensity_cols].isna().sum(axis=0)

    total_missing = df[intensity_cols].isna().sum().sum()
    total_values = df[intensity_cols].size
    missing_pct = 100 * total_missing / total_values

    return {'per_protein': missing_per_protein, 'per_sample': missing_per_sample, 'total_pct': missing_pct}
```

## Related Skills

- **quantification** - Process imported data for quantification
- **peptide-identification** - Identify peptides from raw spectra
- **expression-matrix/counts-ingest** - Similar data loading patterns
