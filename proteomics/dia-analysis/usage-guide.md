# DIA Analysis Usage Guide

## Overview

Data-independent acquisition (DIA) is a mass spectrometry method where all precursors within isolation windows are fragmented simultaneously. This provides comprehensive coverage with fewer missing values than DDA.

## Key Tools

- **DIA-NN** - Fast, accurate, open-source (recommended)
- **Spectronaut** - Commercial, user-friendly
- **MSFragger-DIA** - Part of FragPipe suite
- **OpenSWATH** - OpenMS-based workflow

## Library-Free vs Library-Based

### Library-Free (Recommended for most cases)
- No prior DDA experiments needed
- Uses deep learning (DIA-NN predictor)
- Generates library from data automatically
- Slightly lower sensitivity than library-based

### Library-Based
- Requires pre-built spectral library
- Higher sensitivity for known targets
- Better for targeted panels
- Library from DDA or predicted (Prosit, DeepLC)

## Key Parameters

### Mass Ranges
- `--min-fr-mz 200 --max-fr-mz 1800` - Fragment ion range
- `--min-pr-mz 300 --max-pr-mz 1800` - Precursor range

### Digestion
- `--cut K*,R*` - Trypsin specificity
- `--missed-cleavages 1` - Allow 1 missed cleavage

### Modifications
- `--unimod4` - Carbamidomethyl C (fixed)
- `--var-mod UniMod:35,15.994915,M` - Oxidation M (variable)

### Quality Control
- `--qvalue 0.01` - 1% FDR at precursor and protein level
- `--reanalyse` - Two-pass analysis for MBR
- `--smart-profiling` - Improved quantification

## Typical Workflow

1. Convert raw files to mzML (ProteoWizard)
2. Run DIA-NN (library-free or library-based)
3. Load protein matrix in R/Python
4. Normalize (median centering)
5. Impute missing values
6. Statistical testing (limma)

## References

- DIA-NN: https://github.com/vdemichev/DiaNN
- DIA-NN paper: doi:10.1038/s41592-019-0638-x
