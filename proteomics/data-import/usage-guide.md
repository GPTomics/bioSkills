# Data Import - Usage Guide

## Overview
Load mass-spectrometry data (mzML/mzXML raw spectra, MaxQuant proteinGroups.txt, DIA-NN report.parquet) into Python or R, and -- in the same step -- enforce the two contracts that decide whether downstream numbers mean anything: strip the search engine's bookkeeping (decoys, contaminants, site-only groups, semicolon razor-ID ambiguity) and pick the quant column that matches the question (Intensity vs LFQ intensity vs iBAQ). Import is also where the acquisition mode's missingness structure (DDA MNAR vs DIA MCAR) is inherited and diagnosed.

## Prerequisites
```bash
pip install pyopenms pandas numpy pyarrow
# R alternative (current Bioconductor): BiocManager::install(c("Spectra", "QFeatures"))
# Legacy R (maintenance mode): BiocManager::install("MSnbase")
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my MaxQuant proteinGroups.txt, strip decoys and contaminants, and give me a log2 LFQ matrix"
- "Import the DIA-NN report.parquet and build a protein-by-run matrix at 1% FDR"
- "Read the MS1 and MS2 spectra from my mzML and report precursor and isolation-window info"
- "Tell me whether my missing values look MNAR or MCAR before I impute"

## Example Prompts

### Loading Search Engine Output
> "Load MaxQuant proteinGroups.txt, remove Reverse/contaminant/site-only rows, take the leading protein and gene from the semicolon lists, set zeros to NaN, and log2-transform the LFQ intensities"

> "Import DIA-NN report.parquet, filter Q.Value and PG.Q.Value to 1%, and pivot PG.MaxLFQ into a protein-by-run matrix"

> "I have a proteinGroups.txt with Intensity, LFQ intensity, and iBAQ columns -- which one should I use for comparing two conditions, and why?"

### Loading Raw MS Data
> "Read my mzML with pyOpenMS and iterate spectra by MS level"

> "Extract precursor m/z and isolation-window width from the MS2 spectra in my mzML"

### Diagnosing Missingness
> "Check whether missingness in my MaxQuant matrix correlates with abundance so I know if it is MNAR"

> "My data is DDA -- which imputation methods are legitimate and which will bias my results?"

## What the Agent Will Do
1. Detect the source format (mzML/mzXML, MaxQuant txt, DIA-NN parquet) and load with the matching reader.
2. Strip search-engine bookkeeping: Reverse/REV__ decoys, Potential contaminant/CON__, and Only-identified-by-site rows.
3. Resolve semicolon protein-ID and gene-name lists to the leading (razor) entry, guarding blank gene names.
4. Select the correct quant column for the question (LFQ intensity for between-sample, iBAQ within-sample, Intensity raw).
5. Set MaxQuant zeros to NaN, then log2-transform.
6. For DIA-NN, filter precursor- and protein-group q-values to 1% before pivoting PG.MaxLFQ.
7. Diagnose the missingness pattern (MNAR vs MCAR) and flag which imputation class is legitimate downstream.

## Supported Formats

| Format | Description | Reader |
|--------|-------------|--------|
| mzML / mzXML | Open standards for raw MS data | pyOpenMS, Spectra (R) |
| proteinGroups.txt | MaxQuant protein-group output (LFQ/iBAQ/Intensity; flag columns) | pandas |
| evidence.txt | MaxQuant per-PSM output (MSstats input) | pandas |
| report.parquet | DIA-NN main report (default 1.9+, only default 2.0) | pandas.read_parquet |
| report.tsv | DIA-NN legacy report (pre-2.0) | pandas |

## Tips
- Use `low_memory=False` when reading MaxQuant TSVs to avoid mixed-type warnings.
- `Protein IDs`, `Majority protein IDs`, and `Gene names` are semicolon lists -- take the first entry as the leading/razor identifier.
- `Only identified by site` exists ONLY in proteinGroups.txt; guard the lookup when parsing other tables.
- A MaxQuant zero means "not quantified", not "zero abundance" -- convert to NaN before log2.
- Use `LFQ intensity` for between-sample comparison, `iBAQ` for within-sample molar abundance, raw `Intensity` only for custom normalization.
- DIA-NN 2.0 dropped the TSV default; read `report.parquet` and q-filter before pivoting.
- Diagnose missingness here: a negative correlation between abundance and missingness is the MNAR signature that forbids mean/KNN imputation.

## Related Skills
- peptide-identification - search raw spectra and convert vendor RAW to mzML
- quantification - compute MaxLFQ and TMT reporter-ion quantities from imported data
- protein-inference - resolve protein-group parsimony and razor assignment
- differential-abundance - normalize, impute (per the missingness diagnosis), and test
- proteomics-qc - assess run-level identification and quant quality
- dia-analysis - run DIA-NN to produce the report this skill imports
- expression-matrix/normalization - general intensity-matrix normalization patterns
- workflows/proteomics-pipeline - end-to-end pipeline that begins with this import step
