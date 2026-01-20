# Peptide Identification

## Overview

Peptide identification matches MS/MS spectra to peptide sequences from a protein database or spectral library.

## Search Strategies

### Database Search
1. Digest protein database in-silico
2. Generate theoretical spectra for peptides
3. Match experimental spectra to theoretical
4. Score matches and estimate FDR

### Spectral Library Search
1. Use previously identified spectra as reference
2. Match by spectral similarity
3. Faster but limited to known peptides

## Key Parameters

| Parameter | Typical Value | Effect |
|-----------|---------------|--------|
| Precursor tolerance | 10-20 ppm | Match window for precursor m/z |
| Fragment tolerance | 0.02-0.05 Da | Match window for fragment ions |
| Missed cleavages | 2 | Allow incomplete digestion |
| Fixed modifications | Carbamidomethyl (C) | Always present |
| Variable modifications | Oxidation (M), Phospho (STY) | May be present |

## FDR Control

The target-decoy approach:
1. Search against forward + reversed sequences
2. Decoys estimate false positives
3. FDR = #decoys / #targets at score threshold
4. Typically filter to 1% FDR

## Common Enzymes

| Enzyme | Cleavage Rule |
|--------|---------------|
| Trypsin | After K, R (not before P) |
| Lys-C | After K |
| Chymotrypsin | After F, W, Y, L |
| Asp-N | Before D |

## Output Formats

- **mzIdentML (.mzid)**: Standard exchange format
- **idXML**: OpenMS internal format
- **pepXML**: Trans-Proteomic Pipeline format
