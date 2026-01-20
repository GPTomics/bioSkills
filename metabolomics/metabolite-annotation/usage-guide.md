# Metabolite Annotation Usage Guide

## Overview

Metabolite annotation (identification) assigns chemical identities to detected features based on m/z, retention time, and MS/MS spectra.

## Annotation Approaches

### Level 1: Authentic Standard
- Match m/z, RT, and MS/MS to run standard
- Highest confidence, limited coverage

### Level 2: Database MS/MS Match
- Match MS/MS to reference spectra (HMDB, GNPS, MassBank)
- Good confidence for known compounds

### Level 3: Formula Assignment
- Accurate mass + isotope pattern → molecular formula
- SIRIUS/CSI:FingerID for structure prediction

### Level 4: Mass Match
- m/z matches database entry within tolerance
- Multiple candidates usually

## Key Databases

| Database | Content | Access |
|----------|---------|--------|
| HMDB | Human metabolites | hmdb.ca |
| KEGG | Pathway metabolites | kegg.jp |
| PubChem | All chemicals | pubchem.ncbi.nlm.nih.gov |
| MassBank | MS/MS spectra | massbank.eu |
| GNPS | MS/MS spectra | gnps.ucsd.edu |

## Adduct Considerations

### Positive Mode
- [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+

### Negative Mode
- [M-H]-, [M+Cl]-, [M+FA-H]-

## Best Practices

1. Always specify polarity and adducts
2. Use MS/MS when available
3. Report confidence levels
4. Validate with RT if possible
5. Consider in-source fragments

## References

- MSI Reporting: doi:10.1007/s11306-007-0082-2
- SIRIUS: doi:10.1038/s41592-019-0344-8
- matchms: doi:10.21105/joss.02411
