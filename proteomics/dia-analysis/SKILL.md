---
name: bio-proteomics-dia-analysis
description: Data-independent acquisition (DIA) proteomics analysis with DIA-NN and other tools. Use when analyzing DIA mass spectrometry data with library-free or library-based workflows for deep proteome profiling.
tool_type: cli
primary_tool: diann
---

# DIA Proteomics Analysis

## DIA-NN Library-Free Analysis

```bash
# Library-free mode (generates library from data)
diann \
    --f sample1.mzML \
    --f sample2.mzML \
    --lib "" \
    --threads 8 \
    --verbose 1 \
    --out report.tsv \
    --qvalue 0.01 \
    --matrices \
    --out-lib generated_lib.tsv \
    --gen-spec-lib \
    --predictor \
    --fasta uniprot_human.fasta \
    --fasta-search \
    --min-fr-mz 200 \
    --max-fr-mz 1800 \
    --met-excision \
    --cut K*,R* \
    --missed-cleavages 1 \
    --min-pep-len 7 \
    --max-pep-len 30 \
    --min-pr-mz 300 \
    --max-pr-mz 1800 \
    --min-pr-charge 1 \
    --max-pr-charge 4 \
    --unimod4 \
    --var-mods 1 \
    --var-mod UniMod:35,15.994915,M \
    --reanalyse \
    --smart-profiling
```

## DIA-NN with Spectral Library

```bash
# Use pre-built or predicted library
diann \
    --f sample1.mzML \
    --f sample2.mzML \
    --lib spectral_library.tsv \
    --threads 8 \
    --verbose 1 \
    --out report.tsv \
    --qvalue 0.01 \
    --matrices \
    --reanalyse \
    --smart-profiling
```

## DIA-NN Output Files

```
report.tsv                    # Main quantification report (long format)
report.stats.tsv              # Run statistics
report.pg_matrix.tsv          # Protein group quantities (wide format)
report.pr.matrix.tsv          # Precursor quantities (wide format)
report.gg_matrix.tsv          # Gene group quantities (wide format)
generated_lib.tsv             # Generated spectral library (if requested)
```

## Load DIA-NN Results in R

```r
library(tidyverse)

# Load main report
report <- read_tsv('report.tsv')

# Load protein matrix (already wide format)
proteins <- read_tsv('report.pg_matrix.tsv')

# Filter and reshape for analysis
protein_matrix <- proteins %>%
    column_to_rownames('Protein.Group') %>%
    select(starts_with('sample')) %>%
    as.matrix()

# Log2 transform (DIA-NN outputs raw intensities)
log2_matrix <- log2(protein_matrix)
log2_matrix[is.infinite(log2_matrix)] <- NA
```

## Load DIA-NN Results in Python

```python
import pandas as pd
import numpy as np

# Load main report
report = pd.read_csv('report.tsv', sep='\t')

# Load protein matrix
proteins = pd.read_csv('report.pg_matrix.tsv', sep='\t')
proteins = proteins.set_index('Protein.Group')

# Log2 transform
log2_proteins = np.log2(proteins.replace(0, np.nan))
```

## MSFragger-DIA Analysis

```bash
# MSFragger for DIA (alternative to DIA-NN)
# Requires FragPipe GUI or command-line workflow

# Generate predicted library with EasyPQP
easypqp library \
    --in psm_results.tsv \
    --out library.pqp \
    --psmtsv \
    --rt_reference irt.tsv

# Convert to DIA-NN format
easypqp convert \
    --in library.pqp \
    --out library.tsv \
    --format diann
```

## Spectronaut Export Processing

```r
# Load Spectronaut report
spectronaut <- read_tsv('spectronaut_report.tsv')

# Pivot to protein matrix
protein_matrix <- spectronaut %>%
    select(PG.ProteinGroups, R.FileName, PG.Quantity) %>%
    pivot_wider(names_from = R.FileName, values_from = PG.Quantity) %>%
    column_to_rownames('PG.ProteinGroups')
```

## DIA Quality Metrics

```r
library(tidyverse)

report <- read_tsv('report.tsv')

# Identifications per run
ids_per_run <- report %>%
    group_by(Run) %>%
    summarise(
        precursors = n_distinct(Precursor.Id),
        proteins = n_distinct(Protein.Group),
        genes = n_distinct(Genes)
    )

# Missing value analysis
proteins <- read_tsv('report.pg_matrix.tsv')
protein_values <- proteins %>% select(-Protein.Group)
missing_pct <- colSums(protein_values == 0 | is.na(protein_values)) / nrow(protein_values) * 100
```

## Match Between Runs

```bash
# DIA-NN MBR is automatic with --reanalyse flag
# First pass: identifies peptides per run
# Second pass: transfers IDs between runs

diann \
    --f *.mzML \
    --lib library.tsv \
    --reanalyse \
    --out report_mbr.tsv
```

## DIA vs DDA Comparison

| Feature | DIA | DDA |
|---------|-----|-----|
| Acquisition | All precursors fragmented | Top-N precursors selected |
| Missing values | Lower (5-20%) | Higher (30-50%) |
| Dynamic range | Better for low-abundance | Better for high-abundance |
| Library required | Optional (library-free) | Not applicable |
| Quantification | More reproducible | More variable |
| Analysis tools | DIA-NN, Spectronaut | MaxQuant, MSFragger |

## Related Skills

- data-import - Load raw MS data
- spectral-libraries - Build and use spectral libraries
- quantification - Normalization methods
- differential-abundance - Statistical testing
