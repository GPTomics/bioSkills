# Metabolite Annotation - Usage Guide

## Overview

Metabolite annotation assigns chemical identities to detected features based on m/z, retention time, and MS/MS spectra. Confidence levels range from Level 1 (authenticated standard) to Level 4 (mass match only).

## Prerequisites

```bash
# Python
pip install matchms rdkit-pypi

# R/Bioconductor
BiocManager::install("CompoundDb")

# External tools
# SIRIUS: https://bio.informatik.uni-jena.de/software/sirius/
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my metabolomics features against HMDB"
- "Match MS/MS spectra to MassBank database"

## Example Prompts

### Database Matching
> "Search my features against HMDB with 10 ppm mass tolerance in positive mode"
> "Match m/z values to KEGG compounds considering [M+H]+ and [M+Na]+ adducts"

### MS/MS Annotation
> "Compare my MS/MS spectra to MassBank using cosine similarity"
> "Run SIRIUS for molecular formula prediction on features with MS/MS"

### Adduct Consideration
> "Annotate features considering common positive mode adducts: [M+H]+, [M+Na]+, [M+NH4]+"
> "Check for in-source fragments and adduct clusters"

### Confidence Assignment
> "Assign MSI confidence levels to my annotations based on evidence"
> "Filter annotations to keep only Level 2 or better matches"

## What the Agent Will Do

1. Load feature table with m/z, RT, and optional MS/MS spectra
2. Calculate expected masses for different adducts
3. Search against metabolite databases within mass tolerance
4. Score MS/MS matches if spectra available
5. Assign confidence levels
6. Export annotations with evidence

## Tips

- Always specify ion mode (positive/negative) and expected adducts
- MS/MS matching greatly improves annotation confidence
- Report confidence levels per MSI guidelines
- Consider in-source fragments that may appear as separate features
- Use multiple databases for better coverage

## Annotation Confidence Levels

| Level | Evidence Required |
|-------|-------------------|
| 1 | Authentic standard (m/z, RT, MS/MS) |
| 2 | MS/MS match to database spectra |
| 3 | Molecular formula from accurate mass |
| 4 | Mass match only (multiple candidates) |

## Key Databases

| Database | Content | Access |
|----------|---------|--------|
| HMDB | Human metabolites | hmdb.ca |
| KEGG | Pathway metabolites | kegg.jp |
| MassBank | MS/MS spectra | massbank.eu |
| GNPS | MS/MS spectra | gnps.ucsd.edu |

## References

- MSI Reporting: doi:10.1007/s11306-007-0082-2
- SIRIUS: doi:10.1038/s41592-019-0344-8
- matchms: doi:10.21105/joss.02411
