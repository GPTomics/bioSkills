# Mass Spectrometry Data Import

## Overview

This skill covers loading and parsing mass spectrometry data from various formats used in proteomics workflows.

## Supported Formats

| Format | Description | Tool |
|--------|-------------|------|
| mzML | Open standard for MS data | pyOpenMS, MSnbase |
| mzXML | Legacy open format | pyOpenMS |
| proteinGroups.txt | MaxQuant protein output | pandas |
| evidence.txt | MaxQuant peptide output | pandas |
| report.tsv | DIA-NN output | pandas |

## Common Workflows

### From Raw to Analysis-Ready

1. Load raw mzML files with pyOpenMS
2. Or load preprocessed MaxQuant/DIA-NN output with pandas
3. Filter contaminants and decoys
4. Assess missing values
5. Log-transform intensities

### MaxQuant Filtering

Standard filtering removes:
- **Potential contaminant**: Common contaminants (keratins, trypsin)
- **Reverse**: Decoy sequences for FDR
- **Only identified by site**: Proteins identified only by modification site

### Missing Value Patterns

- **MCAR**: Missing completely at random (rare in proteomics)
- **MAR**: Missing at random (can impute)
- **MNAR**: Missing not at random (low abundance) - use left-censored imputation

## Tips

- Use `low_memory=False` when loading large MaxQuant files
- Check for batch effects in missing value patterns
- Log2-transform intensities before statistical analysis
