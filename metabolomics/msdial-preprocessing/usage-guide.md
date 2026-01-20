# MS-DIAL Preprocessing Usage Guide

## Overview

MS-DIAL is an open-source software for mass spectrometry data processing. It provides a user-friendly GUI alternative to XCMS with built-in annotation capabilities.

## When to Use MS-DIAL

- Prefer GUI over command-line
- Need integrated spectral library matching
- Working with lipidomics data (excellent lipid support)
- Want all-in-one preprocessing + annotation
- Processing data from multiple vendors

## MS-DIAL vs XCMS

| Feature | MS-DIAL | XCMS |
|---------|---------|------|
| Interface | GUI + Console | R package |
| Peak detection | Good | Excellent |
| Alignment | Good | Excellent |
| Annotation | Built-in | Separate |
| Lipidomics | Excellent | Manual |
| Batch processing | Console mode | R scripts |
| Learning curve | Lower | Higher |

## Supported Data Types

- LC-MS DDA (data-dependent acquisition)
- LC-MS DIA (SWATH, MSE, AIF)
- GC-MS (EI)
- LC-IM-MS (ion mobility)
- Imaging MS

## Workflow Steps

### 1. Project Setup
- Choose analysis method (LC-MS, GC-MS, etc.)
- Select positive or negative ion mode
- Define data type (centroid/profile)

### 2. Data Import
- Support for mzML, ABF, mzXML, CDF
- Convert vendor formats first (ProteoWizard)

### 3. Peak Detection
- Automatic parameter optimization available
- Key parameters: minimum peak height, peak width

### 4. Alignment
- Reference file selection (typically QC sample)
- RT tolerance, m/z tolerance settings

### 5. Gap Filling
- Fill missing values from raw data
- Important for complete data matrix

### 6. Identification
- MSP spectral library matching
- LipidBlast for lipidomics
- GNPS, MassBank libraries

### 7. Export
- Alignment result table
- Peak area/height matrix
- Identification results

## Key Parameters

### Peak Detection
| Parameter | Typical Value | Description |
|-----------|---------------|-------------|
| Minimum peak height | 1000-10000 | Intensity threshold |
| Minimum peak width | 3-10 scans | Filter noise |
| Mass accuracy | 5-20 ppm | Instrument dependent |

### Alignment
| Parameter | Typical Value | Description |
|-----------|---------------|-------------|
| RT tolerance | 0.1-0.5 min | Retention time window |
| MS1 tolerance | 0.01-0.05 Da | Mass tolerance |

## Output Files

- `alignment_result.csv` - Main feature table
- `peak_area_matrix.csv` - Areas only
- `msp_output.msp` - Spectra for library

## Common Issues

### Few features detected
- Lower minimum peak height
- Check data quality
- Verify centroid vs profile mode

### Poor alignment
- Use better reference file
- Increase RT tolerance
- Check for large RT shifts

### Missing annotations
- Add more spectral libraries
- Lower identification score cutoff
- Check adduct settings

## Spectral Libraries

- **LipidBlast** - In-silico lipid library
- **MassBank** - Experimental spectra
- **GNPS** - Community library
- **HMDB** - Human metabolome
- **Custom MSP** - User libraries

## References

- MS-DIAL: doi:10.1038/nmeth.4512
- MS-DIAL 4: doi:10.1038/s41592-023-01888-3
- LipidBlast: doi:10.1038/nmeth.2442
