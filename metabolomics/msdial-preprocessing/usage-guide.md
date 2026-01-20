# MS-DIAL Preprocessing - Usage Guide

## Overview

MS-DIAL is an open-source software for mass spectrometry data processing. It provides a user-friendly GUI alternative to XCMS with built-in annotation capabilities and excellent lipidomics support.

## Prerequisites

```bash
# Download MS-DIAL
# https://systemsomicslab.github.io/compms/msdial/main.html

# For data conversion
# ProteoWizard: https://proteowizard.sourceforge.io/

# Python (for exported data)
pip install pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Process my LC-MS data with MS-DIAL and export the feature table"
- "Load MS-DIAL alignment results for downstream analysis"

## Example Prompts

### Data Preparation
> "Convert my vendor raw files to mzML format for MS-DIAL"
> "Help me set up MS-DIAL parameters for positive mode Orbitrap data"

### Processing
> "What peak detection parameters should I use for 10 ppm mass accuracy?"
> "Set up alignment using my QC sample as reference"
> "Configure gap filling to recover missing peak values"

### Annotation
> "Set up LipidBlast library for lipid annotation in MS-DIAL"
> "Configure MS/MS matching against MassBank spectra"

### Export and Analysis
> "Load my MS-DIAL alignment_result.csv into Python for analysis"
> "Parse MS-DIAL output and create a clean feature matrix"
> "Merge MS-DIAL peak table with sample metadata"

## What the Agent Will Do

1. Guide MS-DIAL parameter selection
2. Help interpret processing results
3. Load exported alignment results
4. Parse feature annotations
5. Create analysis-ready data matrices
6. Integrate with downstream statistical analysis

## Tips

- Convert vendor formats to mzML using ProteoWizard before MS-DIAL
- Use a QC sample as alignment reference for best results
- Lower minimum peak height if few features detected
- Check centroid vs profile mode matches your data type
- LipidBlast library is excellent for lipidomics annotation

## MS-DIAL vs XCMS

| Feature | MS-DIAL | XCMS |
|---------|---------|------|
| Interface | GUI + Console | R package |
| Annotation | Built-in | Separate |
| Lipidomics | Excellent | Manual |
| Learning curve | Lower | Higher |
| Batch processing | Console mode | R scripts |

## Key Parameters

| Parameter | Typical Value | Description |
|-----------|---------------|-------------|
| Minimum peak height | 1000-10000 | Intensity threshold |
| Mass accuracy | 5-20 ppm | Instrument dependent |
| RT tolerance | 0.1-0.5 min | Alignment window |

## Output Files

- `alignment_result.csv` - Main feature table with annotations
- `peak_area_matrix.csv` - Peak areas only
- `msp_output.msp` - MS/MS spectra for library

## Spectral Libraries

- **LipidBlast** - In-silico lipid library (built-in)
- **MassBank** - Experimental spectra
- **GNPS** - Community library
- **HMDB** - Human metabolome

## References

- MS-DIAL: doi:10.1038/nmeth.4512
- MS-DIAL 4: doi:10.1038/s41592-023-01888-3
- LipidBlast: doi:10.1038/nmeth.2442
